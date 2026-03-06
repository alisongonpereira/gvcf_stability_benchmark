#!/usr/bin/env bash
# =============================================================================
# NVIDIA Clara Parabricks Benchmark Runner — split-by-chromosome strategy
#
# Pipeline per dataset:
#   [Shared preprocessing]  CombineGVCFs per chr, in parallel  (02_preprocess_by_chr.sh)
#   [Concat combined GVCFs] bcftools concat per-chr combined GVCFs → combined.g.vcf.gz
#   [Parabricks genotyping] pbrun genotypegvcf (GPU) → output.vcf
#
# Metrics captured:
#   wall_time_preprocess_s  – CombineGVCFs phase (shared with GATK)
#   wall_time_concat_s      – bcftools concat of per-chr combined GVCFs
#   wall_time_genotype_s    – pbrun genotypegvcf (GPU)
#   wall_time_s             – total (preprocess + concat + genotype)
#   monitor.json            – CPU/RAM/GPU during concat + genotyping
#   preprocessing/monitor.json – CPU/RAM during preprocessing
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../config.sh"
source "${SCRIPT_DIR}/common.sh"
source "${SCRIPT_DIR}/02_preprocess_by_chr.sh"

: "${PARABRICKS_GPU_DEVICES:=${PARABRICKS_GPU:-0}}"
SOFTWARE="parabricks"

# ─── Prerequisites check ──────────────────────────────────────────────────────
check_prerequisites() {
    if ! command -v docker &>/dev/null; then
        warn "PARABRICKS" "docker not found — skipping Parabricks benchmark"
        exit 0
    fi

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

    if ! command -v bcftools &>/dev/null; then
        warn "PARABRICKS" "bcftools not found — concat step will fail"
    fi

    log "PARABRICKS" "Docker image: ${PARABRICKS_DOCKER_IMAGE}"
    log "PARABRICKS" "Reference:    ${REF_GENOME}"
    log "PARABRICKS" "GPU device:   ${PARABRICKS_GPU_DEVICES}"

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
    local label="dataset_${size}_rep${rep}"

    local dataset_dir="${BENCHMARK_DIR}/01_prep/${label}"
    local manifest="${dataset_dir}/manifest.txt"
    local preproc_dir="${BENCHMARK_DIR}/02_execution/preprocessing/${label}"
    local output_dir="${BENCHMARK_DIR}/02_execution/${SOFTWARE}/${label}"
    local combined_gvcf="${output_dir}/combined.g.vcf.gz"
    local output_vcf="${output_dir}/output.vcf"
    local stderr_log="${output_dir}/stderr.log"
    local monitor_json="${output_dir}/monitor.json"
    local metrics_json="${BENCHMARK_DIR}/03_metrics/metrics_${SOFTWARE}_${size}_rep${rep}.json"
    local tmp_dir="${output_dir}/tmp"

    mkdir -p "${output_dir}" "${tmp_dir}"

    if is_done "${output_dir}"; then
        log "PARABRICKS" "[${label}] Already completed — skipping"
        return 0
    fi

    if [[ ! -f "${manifest}" ]]; then
        warn "PARABRICKS" "[${label}] manifest.txt not found — run preparation first"
        return 1
    fi

    local gvcf_count
    gvcf_count=$(wc -l < "${manifest}")
    log "PARABRICKS" "[${label}] Starting — ${gvcf_count} GVCFs"

    local overall_start; overall_start=$(date +%s)
    local overall_start_iso; overall_start_iso=$(date -u +%Y-%m-%dT%H:%M:%S)
    local exit_code=0

    # ── Step A: Shared preprocessing (CombineGVCFs per chr, parallel) ────────
    log "PARABRICKS" "[${label}] Step A: preprocessing (CombineGVCFs per chr)..."
    run_preprocessing "${size}" "${rep}" "${manifest}" "${preproc_dir}" \
        || { exit_code=$?; warn "PARABRICKS" "[${label}] Preprocessing failed"; }

    local wall_time_preprocess="${PREPROC_WALL_TIME}"
    local n_chrs="${#PREPROC_CHROMOSOMES[@]}"

    if [[ "${exit_code}" -ne 0 || "${n_chrs}" -eq 0 ]]; then
        warn "PARABRICKS" "[${label}] Aborting — preprocessing did not complete"
        local overall_end; overall_end=$(date +%s)
        local overall_end_iso; overall_end_iso=$(date -u +%Y-%m-%dT%H:%M:%S)
        write_metrics_json \
            "${metrics_json}" "${SOFTWARE}" "${size}" "failed" "${exit_code}" \
            "$(( overall_end - overall_start ))" "${overall_start_iso}" "${overall_end_iso}" \
            "${output_vcf}" "false" "0" ""
        return 1
    fi

    log "PARABRICKS" "[${label}] Step A done — ${n_chrs} chrs combined in ${wall_time_preprocess}s"

    # Start resource monitor — covers concat + pbrun steps
    start_monitor "${monitor_json}"
    rm -f "${combined_gvcf}" "${combined_gvcf}.tbi"

    # ── Step B: bcftools concat per-chr combined GVCFs → combined.g.vcf.gz ───
    log "PARABRICKS" "[${label}] Step B: concat per-chr combined GVCFs..."
    local concat_start; concat_start=$(date +%s)

    local chr_combined_gvcfs=()
    for chr in "${PREPROC_CHROMOSOMES[@]}"; do
        chr_combined_gvcfs+=("${preproc_dir}/chr_${chr}.combined.g.vcf.gz")
    done

    bcftools concat -a -D -O z \
        --threads "${BCFTOOLS_THREADS}" \
        -o "${combined_gvcf}" \
        "${chr_combined_gvcfs[@]}" \
        >> "${stderr_log}" 2>&1 \
    && bcftools index -t --threads "${BCFTOOLS_THREADS}" "${combined_gvcf}" \
        >> "${stderr_log}" 2>&1 \
    || { exit_code=$?; warn "PARABRICKS" "[${label}] bcftools concat failed"; }

    local concat_end; concat_end=$(date +%s)
    local wall_time_concat=$(( concat_end - concat_start ))

    # ── Step C: pbrun genotypegvcf (GPU) ──────────────────────────────────────
    local geno_start; geno_start=$(date +%s)

    if [[ "${exit_code}" -eq 0 ]]; then
        log "PARABRICKS" "[${label}] Step C: pbrun genotypegvcf..."
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
        || { exit_code=$?; warn "PARABRICKS" "[${label}] pbrun genotypegvcf failed"; }
    fi

    local geno_end; geno_end=$(date +%s)
    local geno_end_iso; geno_end_iso=$(date -u +%Y-%m-%dT%H:%M:%S)
    local wall_time_genotype=$(( geno_end - geno_start ))
    local wall_time_total=$(( wall_time_preprocess + wall_time_concat + wall_time_genotype ))

    stop_monitor

    # Remove intermediate combined GVCF and tmp
    rm -f "${combined_gvcf}" "${combined_gvcf}.tbi"
    rm -rf "${tmp_dir}"

    # ── Validate + write metrics ──────────────────────────────────────────────
    local status="success"
    local output_valid="false"
    local variant_count=0
    VARIANT_COUNT=0

    if [[ "${exit_code}" -ne 0 ]]; then
        status="failed"
        warn "PARABRICKS" "[${label}] exit_code=${exit_code} — see ${stderr_log}"
    else
        if validate_vcf "${output_vcf}"; then
            output_valid="true"
            variant_count=${VARIANT_COUNT}
        else
            status="invalid_output"
            warn "PARABRICKS" "[${label}] Output VCF validation failed"
        fi
    fi

    write_metrics_json \
        "${metrics_json}" "${SOFTWARE}" "${size}" "${status}" "${exit_code}" \
        "${wall_time_total}" "${overall_start_iso}" "${geno_end_iso}" \
        "${output_vcf}" "${output_valid}" "${variant_count}" \
        "${monitor_json}"

    # Inject per-step timings and preprocessing resource reference
    python3 - "${metrics_json}" \
        "${wall_time_preprocess}" "${wall_time_concat}" "${wall_time_genotype}" \
        "${PREPROC_MONITOR_JSON}" \
        <<'PYEOF'
import json, os, sys
path, t_pre, t_concat, t_geno, pre_monitor = sys.argv[1:]
with open(path) as f:
    d = json.load(f)
d["wall_time_preprocess_s"] = int(t_pre)
d["wall_time_concat_s"]     = int(t_concat)
d["wall_time_genotype_s"]   = int(t_geno)
if pre_monitor and os.path.isfile(pre_monitor):
    with open(pre_monitor) as f:
        mon = json.load(f)
    d["preprocess_resources"] = mon.get("aggregates", {})
with open(path, "w") as f:
    json.dump(d, f, indent=2)
PYEOF

    log "PARABRICKS" "[${label}] Done — status=${status} total=${wall_time_total}s" \
        "(preprocess=${wall_time_preprocess}s + concat=${wall_time_concat}s + genotype=${wall_time_genotype}s)" \
        "variants=${variant_count}"

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
