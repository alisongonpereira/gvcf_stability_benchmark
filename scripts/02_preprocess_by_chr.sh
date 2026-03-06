#!/usr/bin/env bash
# =============================================================================
# Shared Preprocessing: CombineGVCFs per chromosome in parallel
#
# Source this file from 02_run_gatk.sh and 02_run_parabricks.sh.
# Provides: run_preprocessing SIZE REP MANIFEST PREPROC_DIR
#
# After a successful call the following variables are set in the caller:
#   PREPROC_WALL_TIME    – wall-clock seconds for all CombineGVCFs jobs
#   PREPROC_START_ISO    – ISO-8601 start timestamp
#   PREPROC_END_ISO      – ISO-8601 end timestamp
#   PREPROC_MONITOR_JSON – path to the resource-monitor JSON
#   PREPROC_CHROMOSOMES  – array of chromosome names (in .fai order)
#
# Output directory layout:
#   PREPROC_DIR/
#     chr_<CHR>.combined.g.vcf.gz   – per-chr combined GVCF
#     chr_<CHR>.combined.g.vcf.gz.tbi
#     stderr_<CHR>.log              – per-chr GATK stderr
#     monitor.json                  – resource samples during CombineGVCFs
#     timing.json                   – wall-time + status summary
#     .done                         – created on success (idempotency guard)
# =============================================================================

PREPROC_WALL_TIME=0
PREPROC_START_ISO=""
PREPROC_END_ISO=""
PREPROC_MONITOR_JSON=""
PREPROC_CHROMOSOMES=()

# ─── Detect main chromosomes from .fai ────────────────────────────────────────
# Fills caller-supplied nameref array with chromosome names in .fai order.
# Falls back to all contigs if no standard chr names are found.
get_chromosomes() {
    local -n _chrs=$1   # nameref
    _chrs=()

    if [[ ! -f "${REF_GENOME}.fai" ]]; then
        warn "PREPROC" "Cannot detect chromosomes: ${REF_GENOME}.fai not found"
        return 1
    fi

    # NOTE: chrY is intentionally excluded from preprocessing.
    # CombineGVCFs on chrY is pathologically slow with large cohorts due to:
    #   1. Mixed ploidy (haploid in males, diploid in PAR1/PAR2)
    #   2. Massive annotation conflicts when merging male variants with female
    #      <NON_REF> blocks, generating thousands of WARN events per run
    #   3. Runtime scales super-linearly with sample count (observed: 115+ min
    #      for 50 samples vs ~4 min for autosomes of comparable size)
    # chrY variants are not included in any downstream metric or report.
    # Future work: handle chrY separately with non-PAR intervals only
    # (-L chrY:2781480-56887902 for hg38) or with ploidy-aware genotyping.
    while IFS= read -r chr; do
        _chrs+=("${chr}")
    done < <(awk '{print $1}' "${REF_GENOME}.fai" \
             | grep -E '^(chr[0-9]+|chrX|chrM|[0-9]+|X|MT)$')

    if [[ "${#_chrs[@]}" -eq 0 ]]; then
        warn "PREPROC" "No standard chr names found in .fai — using all contigs"
        while IFS= read -r chr; do
            _chrs+=("${chr}")
        done < <(awk '{print $1}' "${REF_GENOME}.fai")
    fi
}

# ─── Main preprocessing function ──────────────────────────────────────────────
# Usage: run_preprocessing SIZE REP MANIFEST PREPROC_DIR
# Returns 0 on success, 1 on any failure.
# Idempotent: if PREPROC_DIR/.done exists, skips computation and reads timing.
run_preprocessing() {
    local size="$1"
    local rep="$2"
    local manifest="$3"
    local preproc_dir="$4"
    local label="dataset_${size}_rep${rep}"

    PREPROC_MONITOR_JSON="${preproc_dir}/monitor.json"

    # ── Already complete → reuse ───────────────────────────────────────────
    if is_done "${preproc_dir}"; then
        log "PREPROC" "[${label}] Already complete — reusing per-chr combined GVCFs"
        if [[ -f "${preproc_dir}/timing.json" ]]; then
            PREPROC_WALL_TIME=$(python3 -c \
                "import json; d=json.load(open('${preproc_dir}/timing.json')); print(d.get('wall_time_s',0))")
            PREPROC_START_ISO=$(python3 -c \
                "import json; d=json.load(open('${preproc_dir}/timing.json')); print(d.get('start_time',''))")
            PREPROC_END_ISO=$(python3 -c \
                "import json; d=json.load(open('${preproc_dir}/timing.json')); print(d.get('end_time',''))")
        fi
        get_chromosomes PREPROC_CHROMOSOMES
        return 0
    fi

    mkdir -p "${preproc_dir}"

    # ── Detect chromosomes ─────────────────────────────────────────────────
    get_chromosomes PREPROC_CHROMOSOMES
    if [[ "${#PREPROC_CHROMOSOMES[@]}" -eq 0 ]]; then
        warn "PREPROC" "[${label}] No chromosomes detected — aborting"
        return 1
    fi
    log "PREPROC" "[${label}] ${#PREPROC_CHROMOSOMES[@]} chromosomes to process"

    # ── Build per-sample -V args (ensure bgzip + tabix) ───────────────────
    local v_args=()
    while IFS= read -r gvcf; do
        [[ -z "${gvcf}" ]] && continue
        local ready_gvcf
        if command -v bgzip &>/dev/null && command -v tabix &>/dev/null; then
            ready_gvcf=$(ensure_bgzipped "${gvcf}")
        else
            ready_gvcf="${gvcf}"
        fi
        v_args+=(-V "${ready_gvcf}")
    done < "${manifest}"

    local tmp_dir="${preproc_dir}/tmp"
    mkdir -p "${tmp_dir}"

    log "PREPROC" "[${label}] Starting CombineGVCFs per chr in parallel..."
    start_monitor "${PREPROC_MONITOR_JSON}"

    PREPROC_START_ISO=$(date -u +%Y-%m-%dT%H:%M:%S)
    local start_epoch; start_epoch=$(date +%s)

    # ── Launch one CombineGVCFs per chromosome ─────────────────────────────
    local pids=()
    for chr in "${PREPROC_CHROMOSOMES[@]}"; do
        local chr_gvcf="${preproc_dir}/chr_${chr}.combined.g.vcf.gz"
        # NOTE: env -u LD_PRELOAD prevents libjemalloc SIGSEGV in JVM (exit 245)
        env -u LD_PRELOAD "${GATK_BIN}" --java-options "${GATK_CHR_JAVA_OPTS}" \
            CombineGVCFs \
            -R "${REF_GENOME}" \
            "${v_args[@]}" \
            -L "${chr}" \
            -O "${chr_gvcf}" \
            --tmp-dir "${tmp_dir}" \
            >> "${preproc_dir}/stderr_${chr}.log" 2>&1 &
        pids+=($!)
    done

    # ── Wait and collect exit codes ────────────────────────────────────────
    local exit_code=0
    for i in "${!pids[@]}"; do
        if wait "${pids[$i]}"; then
            log "PREPROC" "[${label}]   chr ${PREPROC_CHROMOSOMES[$i]}: done"
        else
            warn "PREPROC" "[${label}]   chr ${PREPROC_CHROMOSOMES[$i]}: FAILED"
            exit_code=1
        fi
    done

    local end_epoch; end_epoch=$(date +%s)
    PREPROC_END_ISO=$(date -u +%Y-%m-%dT%H:%M:%S)
    PREPROC_WALL_TIME=$(( end_epoch - start_epoch ))

    stop_monitor
    rm -rf "${tmp_dir}"

    # ── Write timing JSON ──────────────────────────────────────────────────
    python3 - \
        "${preproc_dir}/timing.json" \
        "${size}" "${exit_code}" \
        "${PREPROC_WALL_TIME}" "${PREPROC_START_ISO}" "${PREPROC_END_ISO}" \
        "${#PREPROC_CHROMOSOMES[@]}" "${PREPROC_MONITOR_JSON}" \
        <<'PYEOF'
import json, os, sys
(path, size, exit_code, wall_time,
 start_iso, end_iso, n_chrs, monitor_json) = sys.argv[1:]
d = {
    "stage":        "preprocessing",
    "dataset_size": int(size),
    "exit_code":    int(exit_code),
    "wall_time_s":  int(wall_time),
    "start_time":   start_iso,
    "end_time":     end_iso,
    "n_chromosomes": int(n_chrs),
    "resources":    {},
    "resource_samples": [],
}
if monitor_json and os.path.isfile(monitor_json):
    with open(monitor_json) as f:
        mon = json.load(f)
    d["resources"]        = mon.get("aggregates", {})
    d["resource_samples"] = mon.get("samples", [])
with open(path, "w") as f:
    json.dump(d, f, indent=2)
print(f"[PREPROC] Timing written: {path}")
PYEOF

    if [[ "${exit_code}" -eq 0 ]]; then
        mark_done "${preproc_dir}"
        log "PREPROC" "[${label}] Done — wall_time=${PREPROC_WALL_TIME}s"
    fi

    return "${exit_code}"
}
