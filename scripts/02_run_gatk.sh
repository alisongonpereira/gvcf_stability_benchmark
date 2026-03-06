#!/usr/bin/env bash
# =============================================================================
# GATK Benchmark Runner — split-by-chromosome strategy
#
# Pipeline per dataset:
#   [Shared preprocessing]  CombineGVCFs per chr, in parallel  (02_preprocess_by_chr.sh)
#   [GATK genotyping]       GenotypeGVCFs per chr, in parallel
#   [Concat]                bcftools concat → output.vcf.gz
#
# Metrics captured:
#   wall_time_preprocess_s  – CombineGVCFs phase (shared with Parabricks)
#   wall_time_genotype_s    – GenotypeGVCFs + concat phase
#   wall_time_s             – total (preprocess + genotype)
#   monitor.json            – CPU/RAM during genotyping
#   preprocessing/monitor.json – CPU/RAM during preprocessing
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../config.sh"
source "${SCRIPT_DIR}/common.sh"
source "${SCRIPTS_DIR}/02_preprocess_by_chr.sh"

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

    if [[ ! -f "${REF_GENOME}.fai" ]]; then
        warn "GATK" "Reference index (.fai) not found: ${REF_GENOME}.fai"
        warn "GATK" "Run: samtools faidx ${REF_GENOME}"
    fi

    local dict_file="${REF_GENOME%.fa*}.dict"
    if [[ ! -f "${dict_file}" ]]; then
        warn "GATK" "Sequence dictionary (.dict) not found — GATK may fail"
        warn "GATK" "Run: gatk CreateSequenceDictionary -R ${REF_GENOME}"
    fi

    if ! command -v bcftools &>/dev/null; then
        warn "GATK" "bcftools not found — concat step will fail"
    fi

    log "GATK" "Using: $(command -v ${GATK_BIN})"
    log "GATK" "Reference: ${REF_GENOME}"
    log "GATK" "Java opts (whole-genome): ${GATK_JAVA_OPTS}"
    log "GATK" "Java opts (per-chr):      ${GATK_CHR_JAVA_OPTS}"
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
    local output_vcf="${output_dir}/output.vcf.gz"
    local stderr_log="${output_dir}/stderr.log"
    local monitor_json="${output_dir}/monitor.json"
    local metrics_json="${BENCHMARK_DIR}/03_metrics/metrics_${SOFTWARE}_${size}_rep${rep}.json"
    local tmp_dir="${output_dir}/tmp"

    mkdir -p "${output_dir}" "${tmp_dir}"

    if is_done "${output_dir}"; then
        log "GATK" "[${label}] Already completed — skipping"
        return 0
    fi

    if [[ ! -f "${manifest}" ]]; then
        warn "GATK" "[${label}] manifest.txt not found — run preparation first"
        return 1
    fi

    local gvcf_count
    gvcf_count=$(wc -l < "${manifest}")
    log "GATK" "[${label}] Starting — ${gvcf_count} GVCFs"

    local overall_start; overall_start=$(date +%s)
    local overall_start_iso; overall_start_iso=$(date -u +%Y-%m-%dT%H:%M:%S)
    local exit_code=0

    # ── Step A: Shared preprocessing (CombineGVCFs per chr, parallel) ────────
    log "GATK" "[${label}] Step A: preprocessing (CombineGVCFs per chr)..."
    run_preprocessing "${size}" "${rep}" "${manifest}" "${preproc_dir}" \
        || { exit_code=$?; warn "GATK" "[${label}] Preprocessing failed"; }

    local wall_time_preprocess="${PREPROC_WALL_TIME}"
    local n_chrs="${#PREPROC_CHROMOSOMES[@]}"

    if [[ "${exit_code}" -ne 0 || "${n_chrs}" -eq 0 ]]; then
        warn "GATK" "[${label}] Aborting — preprocessing did not complete"
        local overall_end; overall_end=$(date +%s)
        local overall_end_iso; overall_end_iso=$(date -u +%Y-%m-%dT%H:%M:%S)
        write_metrics_json \
            "${metrics_json}" "${SOFTWARE}" "${size}" "failed" "${exit_code}" \
            "$(( overall_end - overall_start ))" "${overall_start_iso}" "${overall_end_iso}" \
            "${output_vcf}" "false" "0" ""
        return 1
    fi

    log "GATK" "[${label}] Step A done — ${n_chrs} chrs combined in ${wall_time_preprocess}s"

    # ── Step B: GenotypeGVCFs per chr, in parallel ────────────────────────────
    log "GATK" "[${label}] Step B: GenotypeGVCFs per chr (${n_chrs} jobs in parallel)..."
    start_monitor "${monitor_json}"

    local geno_start; geno_start=$(date +%s)
    local geno_start_iso; geno_start_iso=$(date -u +%Y-%m-%dT%H:%M:%S)

    local pids=()
    for chr in "${PREPROC_CHROMOSOMES[@]}"; do
        local chr_combined="${preproc_dir}/chr_${chr}.combined.g.vcf.gz"
        local chr_vcf="${output_dir}/chr_${chr}.vcf.gz"
        env -u LD_PRELOAD "${GATK_BIN}" --java-options "${GATK_CHR_JAVA_OPTS}" \
            GenotypeGVCFs \
            -R "${REF_GENOME}" \
            -V "${chr_combined}" \
            -O "${chr_vcf}" \
            --tmp-dir "${tmp_dir}" \
            >> "${output_dir}/stderr_geno_${chr}.log" 2>&1 &
        pids+=($!)
    done

    for i in "${!pids[@]}"; do
        if wait "${pids[$i]}"; then
            log "GATK" "[${label}]   chr ${PREPROC_CHROMOSOMES[$i]}: genotyped"
        else
            warn "GATK" "[${label}]   chr ${PREPROC_CHROMOSOMES[$i]}: GenotypeGVCFs FAILED"
            exit_code=1
        fi
    done

    # ── Step C: bcftools concat → output.vcf.gz ───────────────────────────────
    if [[ "${exit_code}" -eq 0 ]]; then
        log "GATK" "[${label}] Step C: bcftools concat..."
        local chr_vcfs=()
        for chr in "${PREPROC_CHROMOSOMES[@]}"; do
            chr_vcfs+=("${output_dir}/chr_${chr}.vcf.gz")
        done
        bcftools concat -a -D -O z \
            --threads "${BCFTOOLS_THREADS}" \
            -o "${output_vcf}" \
            "${chr_vcfs[@]}" \
            >> "${stderr_log}" 2>&1 \
        && bcftools index -t --threads "${BCFTOOLS_THREADS}" "${output_vcf}" \
            >> "${stderr_log}" 2>&1 \
        || exit_code=$?

        # Remove per-chr VCFs (keep only final output)
        for chr in "${PREPROC_CHROMOSOMES[@]}"; do
            rm -f "${output_dir}/chr_${chr}.vcf.gz" \
                  "${output_dir}/chr_${chr}.vcf.gz.tbi"
        done
    fi

    local geno_end; geno_end=$(date +%s)
    local geno_end_iso; geno_end_iso=$(date -u +%Y-%m-%dT%H:%M:%S)
    local wall_time_genotype=$(( geno_end - geno_start ))
    local wall_time_total=$(( wall_time_preprocess + wall_time_genotype ))

    stop_monitor
    rm -rf "${tmp_dir}"

    # ── Validate + write metrics ──────────────────────────────────────────────
    local status="success"
    local output_valid="false"
    local variant_count=0
    VARIANT_COUNT=0

    if [[ "${exit_code}" -ne 0 ]]; then
        status="failed"
        warn "GATK" "[${label}] exit_code=${exit_code} — see ${output_dir}/stderr*.log"
    else
        if validate_vcf "${output_vcf}"; then
            output_valid="true"
            variant_count=${VARIANT_COUNT}
        else
            status="invalid_output"
            warn "GATK" "[${label}] Output VCF validation failed"
        fi
    fi

    write_metrics_json \
        "${metrics_json}" "${SOFTWARE}" "${size}" "${status}" "${exit_code}" \
        "${wall_time_total}" "${overall_start_iso}" "${geno_end_iso}" \
        "${output_vcf}" "${output_valid}" "${variant_count}" \
        "${monitor_json}"

    # Inject per-step timings and preprocessing resource reference
    python3 - "${metrics_json}" \
        "${wall_time_preprocess}" "${wall_time_genotype}" \
        "${PREPROC_MONITOR_JSON}" \
        <<'PYEOF'
import json, os, sys
path, t_pre, t_geno, pre_monitor = sys.argv[1:]
with open(path) as f:
    d = json.load(f)
d["wall_time_preprocess_s"] = int(t_pre)
d["wall_time_genotype_s"]   = int(t_geno)
if pre_monitor and os.path.isfile(pre_monitor):
    with open(pre_monitor) as f:
        mon = json.load(f)
    d["preprocess_resources"] = mon.get("aggregates", {})
with open(path, "w") as f:
    json.dump(d, f, indent=2)
PYEOF

    log "GATK" "[${label}] Done — status=${status} total=${wall_time_total}s" \
        "(preprocess=${wall_time_preprocess}s + genotype=${wall_time_genotype}s)" \
        "variants=${variant_count}"

    if [[ "${status}" == "success" ]]; then
        mark_done "${output_dir}"
    fi
}

# ─── Main ─────────────────────────────────────────────────────────────────────
main() {
    log "GATK" "============================================================"
    log "GATK" "GATK Benchmark — sizes: ${DATASET_SIZES[*]}  replicates: ${BENCHMARK_REPLICATES}"
    log "GATK" "============================================================"

    check_prerequisites

    local any_failed=0
    for size in "${DATASET_SIZES[@]}"; do
        for rep in $(seq 1 "${BENCHMARK_REPLICATES}"); do
            run_size "${size}" "${rep}" \
                || { warn "GATK" "dataset_${size}_rep${rep} failed (continuing)"; any_failed=1; }
        done
    done

    log "GATK" "GATK benchmark complete."
    return "${any_failed}"
}

main "$@"
