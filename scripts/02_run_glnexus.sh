#!/usr/bin/env bash
# =============================================================================
# GLnexus Benchmark Runner
# Runs glnexus_cli for each dataset size, collects resource metrics, and
# validates the output VCF.
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../config.sh"
source "${SCRIPT_DIR}/common.sh"

SOFTWARE="glnexus"

# ─── Prerequisites check ──────────────────────────────────────────────────────
check_prerequisites() {
    if ! command -v "${GLNEXUS_BIN}" &>/dev/null; then
        warn "GLNEXUS" "glnexus_cli not found (${GLNEXUS_BIN}) — skipping GLnexus benchmark"
        exit 0
    fi
    if ! command -v bcftools &>/dev/null; then
        warn "GLNEXUS" "bcftools not found — required to convert BCF→VCF; skipping"
        exit 0
    fi
    log "GLNEXUS" "Using: $(command -v ${GLNEXUS_BIN})"
    log "GLNEXUS" "Config preset: ${GLNEXUS_CONFIG}"
    log "GLNEXUS" "Threads: ${GLNEXUS_THREADS}  Mem: ${GLNEXUS_MEM_GB}g  bcftools threads: ${BCFTOOLS_THREADS}"
}

# ─── Per-size runner ──────────────────────────────────────────────────────────
run_size() {
    local size="$1"
    local rep="$2"

    local dataset_dir="${BENCHMARK_DIR}/01_prep/dataset_${size}_rep${rep}"
    local manifest="${dataset_dir}/manifest.txt"
    local output_dir="${BENCHMARK_DIR}/02_execution/${SOFTWARE}/dataset_${size}_rep${rep}"
    local output_vcf="${output_dir}/output.vcf"
    local stderr_log="${output_dir}/stderr.log"
    local monitor_json="${output_dir}/monitor.json"
    local metrics_json="${BENCHMARK_DIR}/03_metrics/metrics_${SOFTWARE}_${size}_rep${rep}.json"
    local work_dir="/tmp/glnexus_bench_${size}_rep${rep}_$$"

    mkdir -p "${output_dir}"

    if is_done "${output_dir}"; then
        log "GLNEXUS" "[dataset_${size}_rep${rep}] Already completed — skipping"
        return 0
    fi

    if [[ ! -f "${manifest}" ]]; then
        warn "GLNEXUS" "[dataset_${size}_rep${rep}] manifest.txt not found — run preparation first"
        return 1
    fi

    local gvcf_count
    gvcf_count=$(wc -l < "${manifest}")
    log "GLNEXUS" "[dataset_${size}_rep${rep}] Starting — ${gvcf_count} GVCFs, config=${GLNEXUS_CONFIG}"

    # Ensure GVCFs are bgzipped + indexed (GLnexus requirement)
    local ready_manifest="${output_dir}/manifest_ready.txt"
    : > "${ready_manifest}"
    while IFS= read -r gvcf; do
        [[ -z "${gvcf}" ]] && continue
        if command -v bgzip &>/dev/null && command -v tabix &>/dev/null; then
            ready_gvcf=$(ensure_bgzipped "${gvcf}")
        else
            ready_gvcf="${gvcf}"
            warn "GLNEXUS" "bgzip/tabix not available — using GVCFs as-is (may fail)"
        fi
        echo "${ready_gvcf}" >> "${ready_manifest}"
    done < "${manifest}"

    # Start resource monitor
    start_monitor "${monitor_json}"

    local start_epoch; start_epoch=$(date +%s)
    local start_iso;   start_iso=$(date -u +%Y-%m-%dT%H:%M:%S)
    local exit_code=0

    rm -rf "${work_dir}"

    log "GLNEXUS" "[dataset_${size}_rep${rep}] Executing glnexus_cli..."
    "${GLNEXUS_BIN}" \
        --config    "${GLNEXUS_CONFIG}" \
        --dir       "${work_dir}" \
        --list      "${ready_manifest}" \
        --threads   "${GLNEXUS_THREADS}" \
        --mem-gbytes "${GLNEXUS_MEM_GB}" \
        2>"${stderr_log}" \
    | bcftools view - -O v --threads "${BCFTOOLS_THREADS}" -o "${output_vcf}" \
        2>>"${stderr_log}" \
    || exit_code=$?

    local end_epoch; end_epoch=$(date +%s)
    local end_iso;   end_iso=$(date -u +%Y-%m-%dT%H:%M:%S)
    local wall_time=$(( end_epoch - start_epoch ))

    stop_monitor

    # Cleanup GLnexus work directory
    rm -rf "${work_dir}"

    local status="success"
    local output_valid="false"
    local variant_count=0
    VARIANT_COUNT=0

    if [[ "${exit_code}" -ne 0 ]]; then
        status="failed"
        warn "GLNEXUS" "[dataset_${size}_rep${rep}] exit_code=${exit_code} — see ${stderr_log}"
    else
        if validate_vcf "${output_vcf}"; then
            output_valid="true"
            variant_count=${VARIANT_COUNT}
        else
            status="invalid_output"
            warn "GLNEXUS" "[dataset_${size}_rep${rep}] Output VCF validation failed"
        fi
    fi

    write_metrics_json \
        "${metrics_json}" "${SOFTWARE}" "${size}" "${status}" "${exit_code}" \
        "${wall_time}"    "${start_iso}" "${end_iso}" \
        "${output_vcf}"   "${output_valid}" "${variant_count}" \
        "${monitor_json}"

    log "GLNEXUS" "[dataset_${size}_rep${rep}] Done — status=${status} wall_time=${wall_time}s variants=${variant_count}"

    if [[ "${status}" == "success" ]]; then
        mark_done "${output_dir}"
    fi
}

# ─── Main ─────────────────────────────────────────────────────────────────────
main() {
    log "GLNEXUS" "============================================================"
    log "GLNEXUS" "GLnexus Benchmark — sizes: ${DATASET_SIZES[*]}  replicates: ${BENCHMARK_REPLICATES}"
    log "GLNEXUS" "============================================================"

    check_prerequisites

    local any_failed=0
    for size in "${DATASET_SIZES[@]}"; do
        for rep in $(seq 1 "${BENCHMARK_REPLICATES}"); do
            run_size "${size}" "${rep}" \
                || { warn "GLNEXUS" "dataset_${size}_rep${rep} failed (continuing)"; any_failed=1; }
        done
    done

    log "GLNEXUS" "GLnexus benchmark complete."
    return "${any_failed}"
}

main "$@"
