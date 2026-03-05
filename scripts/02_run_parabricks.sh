#!/usr/bin/env bash
# =============================================================================
# NVIDIA CLARA Parabricks Benchmark Runner
# Runs pbrun joint_genotyping for each dataset size on the A100 GPU.
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../config.sh"
source "${SCRIPT_DIR}/common.sh"
: "${PARABRICKS_GPU_DEVICES:=${PARABRICKS_GPU:-0}}"

SOFTWARE="parabricks"

# ─── Prerequisites check ──────────────────────────────────────────────────────
check_prerequisites() {
    # if ! command -v "${PARABRICKS_BIN}" &>/dev/null; then
    #     warn "PARABRICKS" "pbrun not found (${PARABRICKS_BIN}) — skipping Parabricks benchmark"
    #     exit 0
    # fi
    if ! command -v docker &>/dev/null; then
        warn "PARABRICKS" "docker not found — skipping Parabricks benchmark"
        exit 0
    fi

    # opcional: checar se a imagem existe localmente (não puxa nada)
    if [[ -n "${PARABRICKS_DOCKER_IMAGE:-}" ]]; then
        if ! docker image inspect "${PARABRICKS_DOCKER_IMAGE}" >/dev/null 2>&1; then
            warn "PARABRICKS" "Docker image not found locally: ${PARABRICKS_DOCKER_IMAGE}"
            exit 0
        fi
    fi

    if [[ -z "${REF_GENOME}" ]]; then
        warn "PARABRICKS" "REF_GENOME not set in config.sh — Parabricks requires a reference"
        exit 0
    fi

    if [[ ! -f "${REF_GENOME}" ]]; then
        warn "PARABRICKS" "REF_GENOME file not found: ${REF_GENOME}"
        exit 0
    fi

    log "PARABRICKS" "Using docker image: ${PARABRICKS_DOCKER_IMAGE}"
    log "PARABRICKS" "Reference: ${REF_GENOME}"
    log "PARABRICKS" "GPU device: ${PARABRICKS_GPU_DEVICES}"

    if command -v nvidia-smi &>/dev/null; then
        local gpu_name
        gpu_name=$(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null \
                   | sed -n "$((PARABRICKS_GPU + 1))p" || echo "Unknown")
        log "PARABRICKS" "GPU: ${gpu_name}"
    fi
}

# ─── Per-size runner ──────────────────────────────────────────────────────────
run_size() {
    local size="$1"
    local rep="$2"

    local dataset_dir="${BENCHMARK_DIR}/01_prep/dataset_${size}_rep${rep}"
    local manifest="${dataset_dir}/manifest.txt"
    local output_dir="${BENCHMARK_DIR}/02_execution/${SOFTWARE}/dataset_${size}_rep${rep}"
    local output_vcf="${output_dir}/output.vcf"
    local combined_gvcf="${output_dir}/combined.g.vcf.gz"
    local stderr_log="${output_dir}/stderr.log"
    local monitor_json="${output_dir}/monitor.json"
    local metrics_json="${BENCHMARK_DIR}/03_metrics/metrics_${SOFTWARE}_${size}_rep${rep}.json"
    local tmp_dir="${output_dir}/tmp"

    mkdir -p "${output_dir}" "${tmp_dir}"

    if is_done "${output_dir}"; then
        log "PARABRICKS" "[dataset_${size}_rep${rep}] Already completed — skipping"
        return 0
    fi

    if [[ ! -f "${manifest}" ]]; then
        warn "PARABRICKS" "[dataset_${size}_rep${rep}] manifest.txt not found — run preparation first"
        return 1
    fi

    local gvcf_count
    gvcf_count=$(wc -l < "${manifest}")
    log "PARABRICKS" "[dataset_${size}_rep${rep}] Starting — ${gvcf_count} GVCFs"

    # Build -V arguments for CombineGVCFs (Parabricks genotypegvcf requires a
    # pre-combined gVCF — passing individual per-sample GVCFs produces a 1-sample output)
    local v_args=()
    while IFS= read -r gvcf; do
        [[ -z "${gvcf}" ]] && continue
        if command -v bgzip &>/dev/null && command -v tabix &>/dev/null; then
            ready_gvcf=$(ensure_bgzipped "${gvcf}")
        else
            ready_gvcf="${gvcf}"
        fi
        v_args+=(-V "${ready_gvcf}")
    done < "${manifest}"

    # Clean up any leftover intermediate files from a previous failed run
    rm -f "${combined_gvcf}" "${combined_gvcf}.tbi"

    local exit_code=0

    start_monitor "${monitor_json}"

    local start_epoch; start_epoch=$(date +%s)
    local start_iso;   start_iso=$(date -u +%Y-%m-%dT%H:%M:%S)

    # ── Step A: CombineGVCFs ─────────────────────────────────────────────────
    # NOTE: LD_PRELOAD (libjemalloc) must be unset before invoking the JVM —
    # libjemalloc + JVM = SIGSEGV (exit 245).
    log "PARABRICKS" "[dataset_${size}_rep${rep}] Step A: CombineGVCFs..."
    env -u LD_PRELOAD "${GATK_BIN}" --java-options "${GATK_JAVA_OPTS}" \
        CombineGVCFs \
        -R "${REF_GENOME}" \
        "${v_args[@]}" \
        -O "${combined_gvcf}" \
        --tmp-dir "${tmp_dir}" \
        >> "${stderr_log}" 2>&1 \
    || exit_code=$?

    local combine_end_epoch; combine_end_epoch=$(date +%s)
    local wall_time_combine=$(( combine_end_epoch - start_epoch ))

    if [[ "${exit_code}" -ne 0 ]]; then
        warn "PARABRICKS" "[dataset_${size}_rep${rep}] CombineGVCFs failed (exit_code=${exit_code}) — see ${stderr_log}"
    else
        log "PARABRICKS" "[dataset_${size}_rep${rep}] Step A done — combine_time=${wall_time_combine}s"
    fi

    # ── Step B: pbrun genotypegvcf ────────────────────────────────────────────
    local genotype_start_epoch; genotype_start_epoch=$(date +%s)

    if [[ "${exit_code}" -eq 0 ]]; then
        log "PARABRICKS" "[dataset_${size}_rep${rep}] Step B: pbrun genotypegvcf..."
        docker run --rm \
            --gpus "\"device=${PARABRICKS_GPU_DEVICES}\"" \
            -v /nfs:/nfs -v /home:/home \
            -v /home/alisongonpereira/raid:/home/alisongonpereira/raid \
            -w "${PWD}" \
            "${PARABRICKS_DOCKER_IMAGE}" \
            pbrun genotypegvcf \
                --ref "${REF_GENOME}" \
                --in-gvcf "${combined_gvcf}" \
                --out-vcf "${output_vcf}" \
                2>>"${stderr_log}" \
        || exit_code=$?
    fi

    local end_epoch; end_epoch=$(date +%s)
    local end_iso;   end_iso=$(date -u +%Y-%m-%dT%H:%M:%S)
    local wall_time=$(( end_epoch - start_epoch ))
    local wall_time_genotype=$(( end_epoch - genotype_start_epoch ))

    stop_monitor

    # Clean up intermediate combined GVCF and tmp (keep final VCF)
    rm -f "${combined_gvcf}" "${combined_gvcf}.tbi"
    rm -rf "${tmp_dir}"

    local status="success"
    local output_valid="false"
    local variant_count=0
    VARIANT_COUNT=0

    if [[ "${exit_code}" -ne 0 ]]; then
        status="failed"
        warn "PARABRICKS" "[dataset_${size}_rep${rep}] exit_code=${exit_code} — see ${stderr_log}"
    else
        if validate_vcf "${output_vcf}"; then
            output_valid="true"
            variant_count=${VARIANT_COUNT}
        else
            status="invalid_output"
            warn "PARABRICKS" "[dataset_${size}_rep${rep}] Output VCF validation failed"
        fi
    fi

    write_metrics_json \
        "${metrics_json}" "${SOFTWARE}" "${size}" "${status}" "${exit_code}" \
        "${wall_time}"    "${start_iso}" "${end_iso}" \
        "${output_vcf}"   "${output_valid}" "${variant_count}" \
        "${monitor_json}"

    # Inject per-step breakdown into the metrics JSON
    python3 - "${metrics_json}" "${wall_time_combine}" "${wall_time_genotype}" <<'PYEOF'
import json, sys
path, combine, genotype = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
with open(path) as f:
    d = json.load(f)
d["wall_time_combine_s"]  = combine
d["wall_time_genotype_s"] = genotype
with open(path, "w") as f:
    json.dump(d, f, indent=2)
PYEOF

    log "PARABRICKS" "[dataset_${size}_rep${rep}] Done — status=${status} wall_time=${wall_time}s (combine=${wall_time_combine}s + genotype=${wall_time_genotype}s) variants=${variant_count}"

    if [[ "${status}" == "success" ]]; then
        mark_done "${output_dir}"
    fi
}

# ─── Main ─────────────────────────────────────────────────────────────────────
main() {
    log "PARABRICKS" "============================================================"
    log "PARABRICKS" "Parabricks Benchmark — sizes: ${DATASET_SIZES[*]}  replicates: ${BENCHMARK_REPLICATES}"
    log "PARABRICKS" "============================================================"

    check_prerequisites

    local any_failed=0
    for size in "${DATASET_SIZES[@]}"; do
        for rep in $(seq 1 "${BENCHMARK_REPLICATES}"); do
            run_size "${size}" "${rep}" \
                || { warn "PARABRICKS" "dataset_${size}_rep${rep} failed (continuing)"; any_failed=1; }
        done
    done

    log "PARABRICKS" "Parabricks benchmark complete."
    return "${any_failed}"
}

main "$@"
