#!/usr/bin/env bash
# =============================================================================
# GVCF Scalability Benchmark — Main Orchestrator
# Tests GLnexus, NVIDIA CLARA Parabricks, and GATK at dataset sizes
# 10, 20, 30 … 100 GVCFs.
#
# Usage:
#   bash run_benchmark.sh [--skip-prep] [--only <software>] [--force]
#
# Options:
#   --skip-prep       Skip dataset preparation (assumes datasets already exist)
#   --only <sw>       Run only the specified software (glnexus/parabricks/gatk)
#   --force           Ignore .done flags and rerun all steps
#   --report-only     Skip all benchmarks; only regenerate reports
#
# Prerequisites:
#   1. Populate ./input_gvcfs/ with ≥100 GVCF files (.g.vcf.gz or .gvcf.gz)
#   2. Edit config.sh to set REF_GENOME and any tool-specific options
#   3. Ensure the relevant bioinformatics tools are in PATH
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Source configuration (sets BENCHMARK_DIR, INPUT_DIR, etc.)
source "${SCRIPT_DIR}/config.sh"
source "${SCRIPT_DIR}/scripts/common.sh"

# ─── Parse arguments ─────────────────────────────────────────────────────────
SKIP_PREP=0
REPORT_ONLY=0
FORCE=0
ONLY_SW=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --skip-prep)   SKIP_PREP=1;    shift ;;
        --report-only) REPORT_ONLY=1;  shift ;;
        --force)       FORCE=1;        shift ;;
        --only)        ONLY_SW="$2";   shift 2 ;;
        -h|--help)
            sed -n '3,20p' "$0"
            exit 0 ;;
        *) echo "Unknown option: $1" >&2; exit 1 ;;
    esac
done

# ─── Bootstrap logging ────────────────────────────────────────────────────────
mkdir -p "${BENCHMARK_DIR}/04_reports"
LOG_FILE="${BENCHMARK_DIR}/04_reports/execution.log"
touch "${LOG_FILE}"

# ─── Setup ────────────────────────────────────────────────────────────────────
setup_directories() {
    log "SETUP" "Creating benchmark directory structure..."
    mkdir -p \
        "${BENCHMARK_DIR}/01_prep" \
        "${BENCHMARK_DIR}/02_execution/glnexus" \
        "${BENCHMARK_DIR}/02_execution/parabricks" \
        "${BENCHMARK_DIR}/02_execution/parabricks_glnexus" \
        "${BENCHMARK_DIR}/02_execution/gatk" \
        "${BENCHMARK_DIR}/02_execution/gatk_genomicsdb" \
        "${BENCHMARK_DIR}/03_metrics" \
        "${BENCHMARK_DIR}/04_reports"
    log "SETUP" "Directory structure ready."
}

# ─── Environment validation ───────────────────────────────────────────────────
validate_environment() {
    log "VALIDATE" "Validating environment..."

    # Input GVCFs
    if [[ ! -d "${INPUT_DIR}" ]]; then
        mkdir -p "${INPUT_DIR}"
        warn "VALIDATE" "Created ${INPUT_DIR} — please populate with GVCF files before running"
    fi

    local gvcf_count
    gvcf_count=$(find "${INPUT_DIR}" -maxdepth 1 \
        \( -name "*.g.vcf" -o -name "*.g.vcf.gz" \
           -o -name "*.gvcf" -o -name "*.gvcf.gz" \) 2>/dev/null | wc -l)
    log "VALIDATE" "GVCF pool: ${gvcf_count} files in ${INPUT_DIR}"

    local max_size="${DATASET_SIZES[-1]}"
    if [[ "${gvcf_count}" -lt "${max_size}" ]] && [[ "${REPORT_ONLY}" -eq 0 ]]; then
        warn "VALIDATE" "Pool has ${gvcf_count} GVCFs; max dataset size is ${max_size}. Larger datasets will be skipped."
    fi

    # Software availability
    for software in "${BENCHMARK_SOFTWARES[@]}"; do
        case "${software}" in
            glnexus)
                command -v "${GLNEXUS_BIN}"    &>/dev/null \
                    && log "VALIDATE" "GLnexus   : found ($(command -v "${GLNEXUS_BIN}"))" \
                    || warn "VALIDATE" "GLnexus   : NOT FOUND — ${GLNEXUS_BIN} (will skip)"
                ;;
            parabricks)
                if command -v docker >/dev/null 2>&1 && \
                    docker run --rm --gpus "device=${PARABRICKS_GPU:-0}" "${PARABRICKS_DOCKER_IMAGE:-nvcr.io/nvidia/clara/clara-parabricks:4.5.1-1}" \
                    which pbrun >/dev/null 2>&1; then
                    log "VALIDATE" "Parabricks: found (docker image ${PARABRICKS_DOCKER_IMAGE})"
                else
                    warn "VALIDATE" "Parabricks: NOT FOUND — docker/pbrun unavailable (will skip)"
                fi
                ;;
            parabricks_glnexus)
                if command -v docker >/dev/null 2>&1 && \
                    docker image inspect "${PARABRICKS_GLNEXUS_DOCKER_IMAGE:-nvcr.io/nvidia/clara/clara-parabricks:3.6.1-1}" >/dev/null 2>&1; then
                    log "VALIDATE" "Parabricks GLnexus: found (docker image ${PARABRICKS_GLNEXUS_DOCKER_IMAGE})"
                else
                    warn "VALIDATE" "Parabricks GLnexus: NOT FOUND — docker image ${PARABRICKS_GLNEXUS_DOCKER_IMAGE} not available (will skip)"
                fi
                ;;
            gatk)
                command -v "${GATK_BIN}"       &>/dev/null \
                    && log "VALIDATE" "GATK      : found ($(command -v "${GATK_BIN}"))" \
                    || warn "VALIDATE" "GATK      : NOT FOUND — ${GATK_BIN} (will skip)"
                ;;
            gatk_genomicsdb)
                command -v "${GATK_BIN}"       &>/dev/null \
                    && log "VALIDATE" "GATK GDB  : found ($(command -v "${GATK_BIN}"))" \
                    || warn "VALIDATE" "GATK GDB  : NOT FOUND — ${GATK_BIN} (will skip)"
                ;;
        esac
    done

    # GPU
    if command -v nvidia-smi &>/dev/null; then
        local gpu_info
        gpu_info=$(nvidia-smi --query-gpu=name,memory.total --format=csv,noheader 2>/dev/null | head -1 || echo "Unknown")
        log "VALIDATE" "GPU: ${gpu_info}"
    else
        warn "VALIDATE" "nvidia-smi not found — GPU metrics will not be collected"
    fi

    # Helper tools
    for tool in bcftools bgzip tabix; do
        command -v "${tool}" &>/dev/null \
            && log "VALIDATE" "${tool}: found" \
            || warn "VALIDATE" "${tool}: NOT FOUND — some steps may fail"
    done

    # Python
    python3 --version &>/dev/null \
        && log "VALIDATE" "Python3: $(python3 --version)" \
        || { err "VALIDATE" "python3 not found — required for preparation and reports"; exit 1; }

    log "VALIDATE" "Validation complete."
}

# ─── Preparation ──────────────────────────────────────────────────────────────
run_preparation() {
    log "PREP" "=== Step 1: Dataset Preparation ==="
    local prep_log="${BENCHMARK_DIR}/preparation.log"

    local force_flag=""
    [[ "${FORCE}" -eq 1 ]] && force_flag="--force"

    python3 "${SCRIPTS_DIR}/01_prepare_datasets.py" \
        --input-dir  "${INPUT_DIR}" \
        --output-dir "${BENCHMARK_DIR}/01_prep" \
        --sizes      "${DATASET_SIZES[@]}" \
        --replicates "${BENCHMARK_REPLICATES}" \
        --seed       "${RANDOM_SEED}" \
        --log        "${prep_log}" \
        ${force_flag}

    log "PREP" "Preparation log: ${prep_log}"
}

# ─── Per-software benchmark ───────────────────────────────────────────────────
run_software() {
    local software="$1"
    log "BENCH" "=== Step 2: ${software} Benchmark ==="

    local runner="${SCRIPTS_DIR}/02_run_${software}.sh"
    if [[ ! -f "${runner}" ]]; then
        warn "BENCH" "Runner not found: ${runner}"
        return 1
    fi

    if bash "${runner}"; then
        log "BENCH" "${software} benchmark completed."
    else
        warn "BENCH" "${software} benchmark finished with errors (check logs)."
    fi
}

# ─── Reports ──────────────────────────────────────────────────────────────────
run_reports() {
    log "REPORT" "=== Step 3a: Variant Quality Analysis ==="

    # env -u LD_PRELOAD: libjemalloc (set for GLnexus) conflicts with bcftools
    # C extensions and openpyxl/lxml.  Always unset before Python invocations.
    env -u LD_PRELOAD python3 "${SCRIPTS_DIR}/03_analyze_variants.py" \
        --benchmark-dir "${BENCHMARK_DIR}" \
    || warn "REPORT" "Variant analysis failed or incomplete (continuing to report)"

    log "REPORT" "=== Step 3b: Generating Performance Report ==="

    env -u LD_PRELOAD python3 "${SCRIPTS_DIR}/04a_generate_reports.py" \
        --benchmark-dir "${BENCHMARK_DIR}" \
        --output-dir    "${BENCHMARK_DIR}/04_reports"

    log "REPORT" "=== Step 3c: Generating Concordance Report ==="

    env -u LD_PRELOAD python3 "${SCRIPTS_DIR}/04b_concordance_report.py" \
        --benchmark-dir "${BENCHMARK_DIR}" \
        --output-dir    "${BENCHMARK_DIR}/04_reports" \
    || warn "REPORT" "Concordance report failed or incomplete (continuing)"

    log "REPORT" "HTML         : ${BENCHMARK_DIR}/04_reports/benchmark_report.html"
    log "REPORT" "Excel        : ${BENCHMARK_DIR}/04_reports/benchmark_data.xlsx"
    log "REPORT" "Concordância : ${BENCHMARK_DIR}/04_reports/concordance_report.html"
}

# ─── Main ─────────────────────────────────────────────────────────────────────
main() {
    local ts_start; ts_start=$(date '+%Y-%m-%d %H:%M:%S')

    log "START" "================================================================="
    log "START" "GVCF Scalability Benchmark"
    log "START" "Started : ${ts_start}"
    log "START" "Sizes      : ${DATASET_SIZES[*]}"
    log "START" "Replicates : ${BENCHMARK_REPLICATES}"
    log "START" "Software   : ${BENCHMARK_SOFTWARES[*]}"
    [[ -n "${ONLY_SW}" ]] && log "START" "Filter  : --only ${ONLY_SW}"
    log "START" "================================================================="

    setup_directories
    validate_environment

    if [[ "${REPORT_ONLY}" -eq 0 ]]; then
        if [[ "${SKIP_PREP}" -eq 0 ]]; then
            run_preparation
        else
            log "START" "Skipping preparation (--skip-prep)"
        fi

        local softwares_to_run=("${BENCHMARK_SOFTWARES[@]}")
        if [[ -n "${ONLY_SW}" ]]; then
            softwares_to_run=("${ONLY_SW}")
        fi

        for sw in "${softwares_to_run[@]}"; do
            run_software "${sw}"
        done
    else
        log "START" "Skipping benchmarks (--report-only)"
    fi

    run_reports

    local ts_end; ts_end=$(date '+%Y-%m-%d %H:%M:%S')
    log "DONE" "================================================================="
    log "DONE" "Benchmark complete!"
    log "DONE" "Started : ${ts_start}"
    log "DONE" "Finished: ${ts_end}"
    log "DONE" "================================================================="

    echo ""
    echo "╔══════════════════════════════════════════════════════════╗"
    echo "║              BENCHMARK COMPLETE                          ║"
    echo "╠══════════════════════════════════════════════════════════╣"
    echo "║  Performance : benchmarks/04_reports/benchmark_report.html"
    echo "║  Concordância: benchmarks/04_reports/concordance_report.html"
    echo "║  Excel       : benchmarks/04_reports/benchmark_data.xlsx"
    echo "║  Full log    : benchmarks/04_reports/execution.log"
    echo "╚══════════════════════════════════════════════════════════╝"
}

main "$@"
