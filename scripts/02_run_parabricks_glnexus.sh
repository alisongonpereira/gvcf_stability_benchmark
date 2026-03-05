#!/usr/bin/env bash
# =============================================================================
# NVIDIA Clara Parabricks — GLnexus GPU Benchmark Runner
# Uses pbrun glnexus (Parabricks 3.6) to joint-genotype GVCFs directly on GPU.
# No CombineGVCFs / GenomicsDBImport pre-step required.
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../config.sh"
source "${SCRIPT_DIR}/common.sh"
: "${PARABRICKS_GPU_DEVICES:=${PARABRICKS_GPU:-0}}"

SOFTWARE="parabricks_glnexus"

# ─── Prerequisites check ──────────────────────────────────────────────────────
check_prerequisites() {
    if ! command -v docker &>/dev/null; then
        warn "PB_GLNEXUS" "docker not found — skipping Parabricks GLnexus benchmark"
        exit 0
    fi

    if [[ -n "${PARABRICKS_GLNEXUS_DOCKER_IMAGE:-}" ]]; then
        if ! docker image inspect "${PARABRICKS_GLNEXUS_DOCKER_IMAGE}" >/dev/null 2>&1; then
            warn "PB_GLNEXUS" "Docker image not found locally: ${PARABRICKS_GLNEXUS_DOCKER_IMAGE}"
            exit 0
        fi
    fi

    if ! command -v bcftools &>/dev/null; then
        warn "PB_GLNEXUS" "bcftools not found — required for BCF→VCF conversion; skipping"
        exit 0
    fi

    log "PB_GLNEXUS" "Using docker image : ${PARABRICKS_GLNEXUS_DOCKER_IMAGE}"
    log "PB_GLNEXUS" "GLnexus config     : ${PARABRICKS_GLNEXUS_CONFIG}"
    log "PB_GLNEXUS" "GPU device         : ${PARABRICKS_GPU_DEVICES}"

    if command -v nvidia-smi &>/dev/null; then
        local gpu_name
        gpu_name=$(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null \
                   | sed -n "$((PARABRICKS_GPU + 1))p" || echo "Unknown")
        log "PB_GLNEXUS" "GPU: ${gpu_name}"
    fi
}

# ─── Per-size runner ──────────────────────────────────────────────────────────
run_size() {
    local size="$1"
    local rep="$2"

    local dataset_dir="${BENCHMARK_DIR}/01_prep/dataset_${size}_rep${rep}"
    local manifest="${dataset_dir}/manifest.txt"
    local output_dir="${BENCHMARK_DIR}/02_execution/${SOFTWARE}/dataset_${size}_rep${rep}"
    local output_bcf="${output_dir}/output.bcf"
    local output_vcf="${output_dir}/output.vcf"
    local stderr_log="${output_dir}/stderr.log"
    local monitor_json="${output_dir}/monitor.json"
    local metrics_json="${BENCHMARK_DIR}/03_metrics/metrics_${SOFTWARE}_${size}_rep${rep}.json"

    mkdir -p "${output_dir}"

    if is_done "${output_dir}"; then
        log "PB_GLNEXUS" "[dataset_${size}_rep${rep}] Already completed — skipping"
        return 0
    fi

    if [[ ! -f "${manifest}" ]]; then
        warn "PB_GLNEXUS" "[dataset_${size}_rep${rep}] manifest.txt not found — run preparation first"
        return 1
    fi

    local gvcf_count
    gvcf_count=$(wc -l < "${manifest}")
    log "PB_GLNEXUS" "[dataset_${size}_rep${rep}] Starting — ${gvcf_count} GVCFs, config=${PARABRICKS_GLNEXUS_CONFIG}"

    # Build --in-gvcf arguments; ensure each GVCF is bgzipped + indexed
    local in_gvcf_args=()
    while IFS= read -r gvcf; do
        [[ -z "${gvcf}" ]] && continue
        if command -v bgzip &>/dev/null && command -v tabix &>/dev/null; then
            ready_gvcf=$(ensure_bgzipped "${gvcf}")
        else
            ready_gvcf="${gvcf}"
        fi
        in_gvcf_args+=(--in-gvcf "${ready_gvcf}")
    done < "${manifest}"

    # Clean up any leftovers from a previous failed run
    rm -f "${output_bcf}" "${output_vcf}"

    local exit_code=0

    start_monitor "${monitor_json}"

    local start_epoch; start_epoch=$(date +%s)
    local start_iso;   start_iso=$(date -u +%Y-%m-%dT%H:%M:%S)

    # ── Step A: pbrun glnexus (GPU) ──────────────────────────────────────────
    log "PB_GLNEXUS" "[dataset_${size}_rep${rep}] Step A: pbrun glnexus..."
    docker run --rm \
        --gpus "\"device=${PARABRICKS_GPU_DEVICES}\"" \
        -v /nfs:/nfs -v /home:/home \
        -v /home/alisongonpereira/raid:/home/alisongonpereira/raid \
        -w "${PWD}" \
        "${PARABRICKS_GLNEXUS_DOCKER_IMAGE}" \
        pbrun glnexus \
            "${in_gvcf_args[@]}" \
            --out-bcf "${output_bcf}" \
            --config  "${PARABRICKS_GLNEXUS_CONFIG}" \
            2>>"${stderr_log}" \
    || exit_code=$?

    # ── Step B: BCF → VCF conversion (CPU) ──────────────────────────────────
    if [[ "${exit_code}" -eq 0 ]]; then
        log "PB_GLNEXUS" "[dataset_${size}_rep${rep}] Step B: BCF→VCF..."
        bcftools view "${output_bcf}" \
            -O v \
            --threads "${BCFTOOLS_THREADS}" \
            -o "${output_vcf}" \
            2>>"${stderr_log}" \
        || exit_code=$?
    fi

    local end_epoch; end_epoch=$(date +%s)
    local end_iso;   end_iso=$(date -u +%Y-%m-%dT%H:%M:%S)
    local wall_time=$(( end_epoch - start_epoch ))

    stop_monitor

    # Keep final VCF, remove intermediate BCF
    rm -f "${output_bcf}"

    local status="success"
    local output_valid="false"
    local variant_count=0
    VARIANT_COUNT=0

    if [[ "${exit_code}" -ne 0 ]]; then
        status="failed"
        warn "PB_GLNEXUS" "[dataset_${size}_rep${rep}] exit_code=${exit_code} — see ${stderr_log}"
    else
        if validate_vcf "${output_vcf}"; then
            output_valid="true"
            variant_count=${VARIANT_COUNT}
        else
            status="invalid_output"
            warn "PB_GLNEXUS" "[dataset_${size}_rep${rep}] Output VCF validation failed"
        fi
    fi

    write_metrics_json \
        "${metrics_json}" "${SOFTWARE}" "${size}" "${status}" "${exit_code}" \
        "${wall_time}"    "${start_iso}" "${end_iso}" \
        "${output_vcf}"   "${output_valid}" "${variant_count}" \
        "${monitor_json}"

    log "PB_GLNEXUS" "[dataset_${size}_rep${rep}] Done — status=${status} wall_time=${wall_time}s variants=${variant_count}"

    if [[ "${status}" == "success" ]]; then
        mark_done "${output_dir}"
    fi
}

# ─── Main ─────────────────────────────────────────────────────────────────────
main() {
    log "PB_GLNEXUS" "============================================================"
    log "PB_GLNEXUS" "Parabricks GLnexus Benchmark — sizes: ${DATASET_SIZES[*]}  replicates: ${BENCHMARK_REPLICATES}"
    log "PB_GLNEXUS" "============================================================"

    check_prerequisites

    local any_failed=0
    for size in "${DATASET_SIZES[@]}"; do
        for rep in $(seq 1 "${BENCHMARK_REPLICATES}"); do
            run_size "${size}" "${rep}" \
                || { warn "PB_GLNEXUS" "dataset_${size}_rep${rep} failed (continuing)"; any_failed=1; }
        done
    done

    log "PB_GLNEXUS" "Parabricks GLnexus benchmark complete."
    return "${any_failed}"
}

main "$@"
