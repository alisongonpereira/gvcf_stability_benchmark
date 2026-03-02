#!/usr/bin/env bash
# =============================================================================
# GATK GenomicsDBImport Benchmark Runner
# For each dataset size: GenomicsDBImport → GenotypeGVCFs
#
# GenomicsDBImport stores variants in a columnar DB rather than a merged GVCF,
# avoiding the I/O-heavy BCF merge step of CombineGVCFs and scaling better
# for large cohorts.
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../config.sh"
source "${SCRIPT_DIR}/common.sh"

SOFTWARE="gatk_genomicsdb"

# ─── Prerequisites check ──────────────────────────────────────────────────────
check_prerequisites() {
    if ! command -v "${GATK_BIN}" &>/dev/null; then
        warn "GATK_GDB" "gatk not found (${GATK_BIN}) — skipping GATK GenomicsDB benchmark"
        exit 0
    fi

    if [[ -z "${REF_GENOME}" ]]; then
        warn "GATK_GDB" "REF_GENOME not set in config.sh — GATK requires a reference"
        exit 0
    fi

    if [[ ! -f "${REF_GENOME}" ]]; then
        warn "GATK_GDB" "REF_GENOME file not found: ${REF_GENOME}"
        exit 0
    fi

    if [[ ! -f "${REF_GENOME}.fai" ]]; then
        warn "GATK_GDB" "Reference index (.fai) not found: ${REF_GENOME}.fai"
        warn "GATK_GDB" "Run: samtools faidx ${REF_GENOME}"
    fi

    local dict_file="${REF_GENOME%.fa*}.dict"
    if [[ ! -f "${dict_file}" ]]; then
        warn "GATK_GDB" "Sequence dictionary (.dict) not found — GATK may fail"
        warn "GATK_GDB" "Run: gatk CreateSequenceDictionary -R ${REF_GENOME}"
    fi

    log "GATK_GDB" "Using: $(command -v ${GATK_BIN})"
    log "GATK_GDB" "Reference: ${REF_GENOME}"
    log "GATK_GDB" "Java opts: ${GATK_JAVA_OPTS}"

    # Log interval source
    if [[ -n "${GENOMICSDB_INTERVALS:-}" ]] && [[ -f "${GENOMICSDB_INTERVALS}" ]]; then
        log "GATK_GDB" "Intervals: ${GENOMICSDB_INTERVALS} (from config)"
    elif [[ -f "${REF_GENOME}.fai" ]]; then
        local n_chrs
        n_chrs=$(awk '{print $1}' "${REF_GENOME}.fai" \
                 | grep -cE '^(chr[0-9]+|chrX|chrY|chrM)$' || true)
        log "GATK_GDB" "Intervals: auto from .fai (${n_chrs} main chromosomes)"
    else
        warn "GATK_GDB" "No .fai found — interval generation may fail"
    fi
}

# ─── Build interval -L arguments ──────────────────────────────────────────────
build_interval_args() {
    # Outputs the -L args array into the caller's l_args variable.
    # Priority:
    #   1. GENOMICSDB_INTERVALS file from config.sh (BED or list)
    #   2. Auto-generate main chromosomes from .fai
    local -n _l_args=$1       # nameref
    local intervals_file="$2" # temp file to write auto-generated list

    _l_args=()

    if [[ -n "${GENOMICSDB_INTERVALS:-}" ]] && [[ -f "${GENOMICSDB_INTERVALS}" ]]; then
        _l_args=(-L "${GENOMICSDB_INTERVALS}")
        return
    fi

    if [[ ! -f "${REF_GENOME}.fai" ]]; then
        warn "GATK_GDB" "Cannot generate intervals: ${REF_GENOME}.fai not found"
        return 1
    fi

    # Extract main chromosomes (handles both "chr1" and "1" naming conventions)
    awk '{print $1}' "${REF_GENOME}.fai" \
        | grep -E '^(chr[0-9]+|chrX|chrY|chrM|[0-9]+|X|Y|MT)$' \
        > "${intervals_file}"

    local n_chrs
    n_chrs=$(wc -l < "${intervals_file}")
    if [[ "${n_chrs}" -eq 0 ]]; then
        warn "GATK_GDB" "No chromosomes matched from .fai — using all contigs"
        awk '{print $1}' "${REF_GENOME}.fai" > "${intervals_file}"
    fi

    while IFS= read -r chr; do
        [[ -z "${chr}" ]] && continue
        _l_args+=(-L "${chr}")
    done < "${intervals_file}"
}

# ─── Per-size runner ──────────────────────────────────────────────────────────
run_size() {
    local size="$1"

    local dataset_dir="${BENCHMARK_DIR}/01_prep/dataset_${size}"
    local manifest="${dataset_dir}/manifest.txt"
    local output_dir="${BENCHMARK_DIR}/02_execution/${SOFTWARE}/dataset_${size}"
    local output_vcf="${output_dir}/output.vcf.gz"
    local workspace="${output_dir}/genomicsdb_workspace"
    local intervals_file="${output_dir}/intervals.list"
    local stderr_log="${output_dir}/stderr.log"
    local monitor_json="${output_dir}/monitor.json"
    local metrics_json="${BENCHMARK_DIR}/03_metrics/metrics_${SOFTWARE}_${size}.json"
    local tmp_dir="${output_dir}/tmp"

    mkdir -p "${output_dir}" "${tmp_dir}"

    if is_done "${output_dir}"; then
        log "GATK_GDB" "[dataset_${size}] Already completed — skipping"
        return 0
    fi

    if [[ ! -f "${manifest}" ]]; then
        warn "GATK_GDB" "[dataset_${size}] manifest.txt not found — run preparation first"
        return 1
    fi

    local gvcf_count
    gvcf_count=$(wc -l < "${manifest}")
    log "GATK_GDB" "[dataset_${size}] Starting — ${gvcf_count} GVCFs"

    # Clean up any leftovers from a previous failed run
    rm -rf "${workspace}" "${output_vcf}" "${output_vcf}.tbi"

    # Ensure bgzipped + tabix-indexed (GenomicsDBImport requirement)
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

    # Build interval arguments
    local l_args=()
    build_interval_args l_args "${intervals_file}"
    if [[ "${#l_args[@]}" -eq 0 ]]; then
        warn "GATK_GDB" "[dataset_${size}] No intervals — aborting"
        return 1
    fi
    log "GATK_GDB" "[dataset_${size}] Using ${#l_args[@]} interval argument(s)"

    # Start resource monitor (covers both GATK steps)
    start_monitor "${monitor_json}"

    local start_epoch; start_epoch=$(date +%s)
    local start_iso;   start_iso=$(date -u +%Y-%m-%dT%H:%M:%S)
    local exit_code=0

    # ── Step A: GenomicsDBImport ──────────────────────────────────────────────
    # NOTE: env -u LD_PRELOAD — same jemalloc + JVM SIGSEGV guard as 02_run_gatk.sh
    log "GATK_GDB" "[dataset_${size}] Step A: GenomicsDBImport..."
    env -u LD_PRELOAD "${GATK_BIN}" --java-options "${GATK_JAVA_OPTS}" \
        GenomicsDBImport \
        "${v_args[@]}" \
        --genomicsdb-workspace-path "${workspace}" \
        "${l_args[@]}" \
        --reader-threads "${GENOMICSDB_READER_THREADS:-4}" \
        --batch-size     "${GENOMICSDB_BATCH_SIZE:-50}" \
        --tmp-dir        "${tmp_dir}" \
        >> "${stderr_log}" 2>&1 \
    || exit_code=$?

    if [[ "${exit_code}" -ne 0 ]]; then
        warn "GATK_GDB" "[dataset_${size}] GenomicsDBImport failed (exit_code=${exit_code}) — see ${stderr_log}"
    fi

    # ── Step B: GenotypeGVCFs (via gendb://) ─────────────────────────────────
    if [[ "${exit_code}" -eq 0 ]]; then
        log "GATK_GDB" "[dataset_${size}] Step B: GenotypeGVCFs (gendb)..."
        env -u LD_PRELOAD "${GATK_BIN}" --java-options "${GATK_JAVA_OPTS}" \
            GenotypeGVCFs \
            -R  "${REF_GENOME}" \
            -V  "gendb://${workspace}" \
            -O  "${output_vcf}" \
            --tmp-dir "${tmp_dir}" \
            >> "${stderr_log}" 2>&1 \
        || exit_code=$?
    fi

    local end_epoch; end_epoch=$(date +%s)
    local end_iso;   end_iso=$(date -u +%Y-%m-%dT%H:%M:%S)
    local wall_time=$(( end_epoch - start_epoch ))

    stop_monitor

    # Remove workspace and tmp (keep final VCF)
    rm -rf "${workspace}" "${tmp_dir}"

    local status="success"
    local output_valid="false"
    local variant_count=0
    VARIANT_COUNT=0

    if [[ "${exit_code}" -ne 0 ]]; then
        status="failed"
        warn "GATK_GDB" "[dataset_${size}] exit_code=${exit_code} — see ${stderr_log}"
    else
        if validate_vcf "${output_vcf}"; then
            output_valid="true"
            variant_count=${VARIANT_COUNT}
        else
            status="invalid_output"
            warn "GATK_GDB" "[dataset_${size}] Output VCF validation failed"
        fi
    fi

    write_metrics_json \
        "${metrics_json}" "${SOFTWARE}" "${size}" "${status}" "${exit_code}" \
        "${wall_time}"    "${start_iso}" "${end_iso}" \
        "${output_vcf}"   "${output_valid}" "${variant_count}" \
        "${monitor_json}"

    log "GATK_GDB" "[dataset_${size}] Done — status=${status} wall_time=${wall_time}s variants=${variant_count}"

    if [[ "${status}" == "success" ]]; then
        mark_done "${output_dir}"
    fi
}

# ─── Main ─────────────────────────────────────────────────────────────────────
main() {
    log "GATK_GDB" "============================================================"
    log "GATK_GDB" "GATK GenomicsDB Benchmark — sizes: ${DATASET_SIZES[*]}"
    log "GATK_GDB" "============================================================"

    check_prerequisites

    local any_failed=0
    for size in "${DATASET_SIZES[@]}"; do
        run_size "${size}" || { warn "GATK_GDB" "dataset_${size} failed (continuing)"; any_failed=1; }
    done

    log "GATK_GDB" "GATK GenomicsDB benchmark complete."
    return "${any_failed}"
}

main "$@"
