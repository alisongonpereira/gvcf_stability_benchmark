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

    local dataset_dir="${BENCHMARK_DIR}/01_prep/dataset_${size}"
    local manifest="${dataset_dir}/manifest.txt"
    local output_dir="${BENCHMARK_DIR}/02_execution/${SOFTWARE}/dataset_${size}"
    local output_vcf="${output_dir}/output.vcf"
    local stderr_log="${output_dir}/stderr.log"
    local monitor_json="${output_dir}/monitor.json"
    local metrics_json="${BENCHMARK_DIR}/03_metrics/metrics_${SOFTWARE}_${size}.json"

    mkdir -p "${output_dir}"

    if is_done "${output_dir}"; then
        log "PARABRICKS" "[dataset_${size}] Already completed — skipping"
        return 0
    fi

    if [[ ! -f "${manifest}" ]]; then
        warn "PARABRICKS" "[dataset_${size}] manifest.txt not found — run preparation first"
        return 1
    fi

    local gvcf_count
    gvcf_count=$(wc -l < "${manifest}")
    log "PARABRICKS" "[dataset_${size}] Starting — ${gvcf_count} GVCFs"

    # Build --in-gvcf arguments
    local ingvcf_args=()
    while IFS= read -r gvcf; do
        [[ -z "${gvcf}" ]] && continue
        # Parabricks works with bgzipped GVCFs
        if command -v bgzip &>/dev/null && command -v tabix &>/dev/null; then
            ready_gvcf=$(ensure_bgzipped "${gvcf}")
        else
            ready_gvcf="${gvcf}"
        fi
        ingvcf_args+=(--in-gvcf "${ready_gvcf}")
    done < "${manifest}"

    # Start resource monitor
    start_monitor "${monitor_json}"

    local start_epoch; start_epoch=$(date +%s)
    local start_iso;   start_iso=$(date -u +%Y-%m-%dT%H:%M:%S)
    local exit_code=0

    log "PARABRICKS" "[dataset_${size}] Executing pbrun joint_genotyping..."
    # "${PARABRICKS_BIN}" joint_genotyping \
    #     --ref     "${REF_GENOME}" \
    #     "${ingvcf_args[@]}" \
    #     --out-vcf "${output_vcf}" \
    #     --gpus '"device=0,1,2,4"' \
    #     2>"${stderr_log}" \
    # || exit_code=$?
    docker run --rm \
    --gpus "\"device=${PARABRICKS_GPU_DEVICES}\"" \
    -v /nfs:/nfs -v /home:/home \
    -v /home/alisongonpereira/raid:/home/alisongonpereira/raid \
    -w "${PWD}" \
    "${PARABRICKS_DOCKER_IMAGE}" \
    pbrun genotypegvcf \
        --ref "${REF_GENOME}" \
        "${ingvcf_args[@]}" \
        --out-vcf "${output_vcf}" \
        2>"${stderr_log}" \
    || exit_code=$?

    local end_epoch; end_epoch=$(date +%s)
    local end_iso;   end_iso=$(date -u +%Y-%m-%dT%H:%M:%S)
    local wall_time=$(( end_epoch - start_epoch ))

    stop_monitor

    local status="success"
    local output_valid="false"
    local variant_count=0
    VARIANT_COUNT=0

    if [[ "${exit_code}" -ne 0 ]]; then
        status="failed"
        warn "PARABRICKS" "[dataset_${size}] exit_code=${exit_code} — see ${stderr_log}"
    else
        if validate_vcf "${output_vcf}"; then
            output_valid="true"
            variant_count=${VARIANT_COUNT}
        else
            status="invalid_output"
            warn "PARABRICKS" "[dataset_${size}] Output VCF validation failed"
        fi
    fi

    write_metrics_json \
        "${metrics_json}" "${SOFTWARE}" "${size}" "${status}" "${exit_code}" \
        "${wall_time}"    "${start_iso}" "${end_iso}" \
        "${output_vcf}"   "${output_valid}" "${variant_count}" \
        "${monitor_json}"

    log "PARABRICKS" "[dataset_${size}] Done — status=${status} wall_time=${wall_time}s variants=${variant_count}"

    if [[ "${status}" == "success" ]]; then
        mark_done "${output_dir}"
    fi
}

# ─── Main ─────────────────────────────────────────────────────────────────────
main() {
    log "PARABRICKS" "============================================================"
    log "PARABRICKS" "Parabricks Benchmark — sizes: ${DATASET_SIZES[*]}"
    log "PARABRICKS" "============================================================"

    check_prerequisites

    local any_failed=0
    for size in "${DATASET_SIZES[@]}"; do
        run_size "${size}" || { warn "PARABRICKS" "dataset_${size} failed (continuing)"; any_failed=1; }
    done

    log "PARABRICKS" "Parabricks benchmark complete."
    return "${any_failed}"
}

main "$@"
