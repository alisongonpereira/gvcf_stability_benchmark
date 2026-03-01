#!/usr/bin/env bash
# =============================================================================
# GATK Benchmark Runner
# For each dataset size: CombineGVCFs → GenotypeGVCFs
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../config.sh"
source "${SCRIPT_DIR}/common.sh"

SOFTWARE="gatk"

# ─── Prerequisites check ──────────────────────────────────────────────────────
check_prerequisites() {
    if ! command -v "${GATK_BIN}" &>/dev/null; then
        warn "GATK" "gatk not found (${GATK_BIN}) — skipping GATK benchmark"
        exit 0
    fi

    if [[ -z "${REF_GENOME}" ]]; then
        warn "GATK" "REF_GENOME not set in config.sh — GATK requires a reference"
        exit 0
    fi

    if [[ ! -f "${REF_GENOME}" ]]; then
        warn "GATK" "REF_GENOME file not found: ${REF_GENOME}"
        exit 0
    fi

    # GATK needs .fai and .dict
    if [[ ! -f "${REF_GENOME}.fai" ]]; then
        warn "GATK" "Reference index (.fai) not found: ${REF_GENOME}.fai"
        warn "GATK" "Run: samtools faidx ${REF_GENOME}"
    fi

    local dict_file="${REF_GENOME%.fa*}.dict"
    if [[ ! -f "${dict_file}" ]]; then
        warn "GATK" "Sequence dictionary (.dict) not found — GATK may fail"
        warn "GATK" "Run: gatk CreateSequenceDictionary -R ${REF_GENOME}"
    fi

    log "GATK" "Using: $(command -v ${GATK_BIN})"
    log "GATK" "Reference: ${REF_GENOME}"
    log "GATK" "Java opts: ${GATK_JAVA_OPTS}"
}

# ─── Per-size runner ──────────────────────────────────────────────────────────
run_size() {
    local size="$1"

    local dataset_dir="${BENCHMARK_DIR}/01_prep/dataset_${size}"
    local manifest="${dataset_dir}/manifest.txt"
    local output_dir="${BENCHMARK_DIR}/02_execution/${SOFTWARE}/dataset_${size}"
    local output_vcf="${output_dir}/output.vcf.gz"
    local combined_gvcf="${output_dir}/combined.g.vcf.gz"
    local stderr_log="${output_dir}/stderr.log"
    local monitor_json="${output_dir}/monitor.json"
    local metrics_json="${BENCHMARK_DIR}/03_metrics/metrics_${SOFTWARE}_${size}.json"
    local tmp_dir="${output_dir}/tmp"

    mkdir -p "${output_dir}" "${tmp_dir}"

    if is_done "${output_dir}"; then
        log "GATK" "[dataset_${size}] Already completed — skipping"
        return 0
    fi

    if [[ ! -f "${manifest}" ]]; then
        warn "GATK" "[dataset_${size}] manifest.txt not found — run preparation first"
        return 1
    fi

    local gvcf_count
    gvcf_count=$(wc -l < "${manifest}")
    log "GATK" "[dataset_${size}] Starting — ${gvcf_count} GVCFs"

    # Ensure bgzipped + indexed (GATK requires .tbi)
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

    # Start resource monitor (before both GATK steps)
    start_monitor "${monitor_json}"

    local start_epoch; start_epoch=$(date +%s)
    local start_iso;   start_iso=$(date -u +%Y-%m-%dT%H:%M:%S)
    local exit_code=0

    # ── Step A: CombineGVCFs ─────────────────────────────────────────────────
    log "GATK" "[dataset_${size}] Step A: CombineGVCFs..."
    "${GATK_BIN}" --java-options "${GATK_JAVA_OPTS}" \
        CombineGVCFs \
        -R "${REF_GENOME}" \
        "${v_args[@]}" \
        -O "${combined_gvcf}" \
        --tmp-dir "${tmp_dir}" \
        2>>"${stderr_log}" \
    || exit_code=$?

    if [[ "${exit_code}" -ne 0 ]]; then
        warn "GATK" "[dataset_${size}] CombineGVCFs failed (exit_code=${exit_code})"
    fi

    # ── Step B: GenotypeGVCFs ────────────────────────────────────────────────
    if [[ "${exit_code}" -eq 0 ]]; then
        log "GATK" "[dataset_${size}] Step B: GenotypeGVCFs..."
        "${GATK_BIN}" --java-options "${GATK_JAVA_OPTS}" \
            GenotypeGVCFs \
            -R "${REF_GENOME}" \
            -V "${combined_gvcf}" \
            -O "${output_vcf}" \
            --tmp-dir "${tmp_dir}" \
            2>>"${stderr_log}" \
        || exit_code=$?
    fi

    local end_epoch; end_epoch=$(date +%s)
    local end_iso;   end_iso=$(date -u +%Y-%m-%dT%H:%M:%S)
    local wall_time=$(( end_epoch - start_epoch ))

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
        warn "GATK" "[dataset_${size}] exit_code=${exit_code} — see ${stderr_log}"
    else
        # output is .vcf.gz — validate via bcftools
        if validate_vcf "${output_vcf}"; then
            output_valid="true"
            variant_count=${VARIANT_COUNT}
        else
            status="invalid_output"
            warn "GATK" "[dataset_${size}] Output VCF validation failed"
        fi
    fi

    write_metrics_json \
        "${metrics_json}" "${SOFTWARE}" "${size}" "${status}" "${exit_code}" \
        "${wall_time}"    "${start_iso}" "${end_iso}" \
        "${output_vcf}"   "${output_valid}" "${variant_count}" \
        "${monitor_json}"

    log "GATK" "[dataset_${size}] Done — status=${status} wall_time=${wall_time}s variants=${variant_count}"

    if [[ "${status}" == "success" ]]; then
        mark_done "${output_dir}"
    fi
}

# ─── Main ─────────────────────────────────────────────────────────────────────
main() {
    log "GATK" "============================================================"
    log "GATK" "GATK Benchmark — sizes: ${DATASET_SIZES[*]}"
    log "GATK" "============================================================"

    check_prerequisites

    local any_failed=0
    for size in "${DATASET_SIZES[@]}"; do
        run_size "${size}" || { warn "GATK" "dataset_${size} failed (continuing)"; any_failed=1; }
    done

    log "GATK" "GATK benchmark complete."
    return "${any_failed}"
}

main "$@"
