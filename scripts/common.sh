#!/usr/bin/env bash
# =============================================================================
# Common helpers — source this from every benchmark script
# =============================================================================

: "${LOG_FILE:=/dev/stderr}"
: "${MONITOR_INTERVAL:=5}"

# ─── Logging ──────────────────────────────────────────────────────────────────
log() {
    local step="$1"; shift
    local ts; ts=$(date '+%Y-%m-%d %H:%M:%S')
    local msg="[${ts}] [${step}] $*"
    echo "${msg}"
    echo "${msg}" >> "${LOG_FILE}" 2>/dev/null || true
}

warn() {
    local step="$1"; shift
    local ts; ts=$(date '+%Y-%m-%d %H:%M:%S')
    local msg="[${ts}] [WARN:${step}] $*"
    echo "${msg}" >&2
    echo "${msg}" >> "${LOG_FILE}" 2>/dev/null || true
}

err() {
    local step="$1"; shift
    local ts; ts=$(date '+%Y-%m-%d %H:%M:%S')
    local msg="[${ts}] [ERROR:${step}] $*"
    echo "${msg}" >&2
    echo "${msg}" >> "${LOG_FILE}" 2>/dev/null || true
}

# ─── Rerun support ────────────────────────────────────────────────────────────
is_done() { [[ -f "$1/.done" ]]; }
mark_done() { date '+%Y-%m-%d %H:%M:%S' > "$1/.done"; }

# ─── Resource monitor ─────────────────────────────────────────────────────────
MONITOR_PID=""

start_monitor() {
    local output_file="$1"
    local monitor_script
    monitor_script="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/monitor_resources.py"
    python3 "${monitor_script}" "${output_file}" --interval "${MONITOR_INTERVAL}" &
    MONITOR_PID=$!
    log "MONITOR" "Resource monitor started (PID=${MONITOR_PID}) → ${output_file}"
}

stop_monitor() {
    if [[ -n "${MONITOR_PID}" ]] && kill -0 "${MONITOR_PID}" 2>/dev/null; then
        kill -TERM "${MONITOR_PID}"
        wait "${MONITOR_PID}" 2>/dev/null || true
        log "MONITOR" "Resource monitor stopped (PID=${MONITOR_PID})"
    fi
    MONITOR_PID=""
}

# ─── VCF / BCF validation ─────────────────────────────────────────────────────
# Sets VARIANT_COUNT; returns 0 on valid, 1 on invalid
VARIANT_COUNT=0
validate_vcf() {
    local vcf_file="$1"
    VARIANT_COUNT=0

    if [[ ! -f "${vcf_file}" ]]; then
        warn "VALIDATE" "Output file not found: ${vcf_file}"
        return 1
    fi

    local size
    size=$(stat -c%s "${vcf_file}" 2>/dev/null || echo 0)
    if [[ "${size}" -eq 0 ]]; then
        warn "VALIDATE" "Output file is empty: ${vcf_file}"
        return 1
    fi

    if command -v bcftools &>/dev/null; then
        local count
        count=$(bcftools view -H "${vcf_file}" 2>/dev/null | wc -l || echo -1)
        if [[ "${count}" -ge 0 ]]; then
            VARIANT_COUNT="${count}"
            log "VALIDATE" "Valid VCF with ${count} variants: $(basename "${vcf_file}")"
            return 0
        else
            warn "VALIDATE" "bcftools could not parse: ${vcf_file}"
            return 1
        fi
    else
        # Fallback without bcftools
        if grep -q "^#CHROM" "${vcf_file}" 2>/dev/null; then
            VARIANT_COUNT=$(grep -cv "^#" "${vcf_file}" 2>/dev/null || echo 0)
            log "VALIDATE" "Valid VCF with ${VARIANT_COUNT} variants (basic check)"
            return 0
        fi
        warn "VALIDATE" "Cannot validate without bcftools"
        return 1
    fi
}

# ─── GVCF bgzip + tabix ───────────────────────────────────────────────────────
# Ensures a GVCF is bgzipped and tabix-indexed.
# Prints the path to the ready-to-use file.
ensure_bgzipped() {
    local gvcf="$1"
    local ready

    if [[ "${gvcf}" == *.gz ]]; then
        ready="${gvcf}"
    else
        ready="${gvcf}.gz"
        if [[ ! -f "${ready}" ]]; then
            log "BGZIP" "bgzipping ${gvcf}"
            bgzip -c "${gvcf}" > "${ready}"
        fi
    fi

    if [[ ! -f "${ready}.tbi" ]]; then
        log "TABIX" "Indexing ${ready}"
        tabix -p vcf "${ready}"
    fi

    echo "${ready}"
}

# ─── Metrics JSON writer ──────────────────────────────────────────────────────
# write_metrics_json OUTPUT_JSON software size status exit_code
#   wall_time_s start_iso end_iso output_vcf output_valid variant_count monitor_json
write_metrics_json() {
    local json_path="$1"
    local software="$2"
    local size="$3"
    local status="$4"
    local exit_code="$5"
    local wall_time_s="$6"
    local start_time="$7"
    local end_time="$8"
    local output_vcf="$9"
    local output_valid="${10}"
    local variant_count="${11}"
    local monitor_json="${12:-}"

    python3 - \
        "${json_path}" "${software}" "${size}" "${status}" "${exit_code}" \
        "${wall_time_s}" "${start_time}" "${end_time}" \
        "${output_vcf}" "${output_valid}" "${variant_count}" \
        "${monitor_json}" \
        <<'PYEOF'
import json, os, sys

(json_path, software, size, status, exit_code,
 wall_time_s, start_time, end_time,
 output_vcf, output_valid, variant_count, monitor_json) = sys.argv[1:]

data = {
    "software": software,
    "dataset_size": int(size),
    "status": status,
    "exit_code": int(exit_code),
    "wall_time_s": float(wall_time_s),
    "start_time": start_time,
    "end_time": end_time,
    "output_vcf": output_vcf,
    "output_valid": output_valid.lower() == "true",
    "variant_count": int(variant_count),
    "resources": {},
    "resource_samples": [],
}

if monitor_json and os.path.isfile(monitor_json):
    with open(monitor_json) as f:
        mon = json.load(f)
    data["resources"] = mon.get("aggregates", {})
    data["resource_samples"] = mon.get("samples", [])

with open(json_path, "w") as f:
    json.dump(data, f, indent=2)
print(f"[METRICS] Written: {json_path}")
PYEOF
}
