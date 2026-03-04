#!/usr/bin/env python3
"""
Step 4c — Variant Type & Quality Analysis

Extracts per-run SNP/Indel breakdown stratified by GQ window, Ti/Tv ratio
curves, indel size distribution, het/hom ratio, and GQ yield curves.

Outputs to benchmarks/03_metrics/:
    variant_types_{software}_{size}_rep{r}.json

Outputs to --output-dir:
    variant_analysis_report.html

Usage:
    python3 04c_variant_analysis.py \\
        --benchmark-dir benchmarks \\
        --output-dir    benchmarks/04_reports \\
        [--force] [--skip-analysis]

Requirements:
    bcftools in PATH
"""

import argparse
import json
import math
import os
import re
import subprocess
import sys
from datetime import datetime
from pathlib import Path

# ─── Configuration ────────────────────────────────────────────────────────────

SOFTWARES  = ["glnexus", "parabricks", "gatk", "gatk_genomicsdb"]
SIZES      = [10, 25, 50, 75, 100]
REPLICATES = [1, 2, 3]

SOFT_LABEL = {
    "glnexus":         "GLnexus",
    "parabricks":      "Parabricks",
    "gatk":            "GATK (CombineGVCFs)",
    "gatk_genomicsdb": "GATK (GenomicsDB)",
}
SW_COLOR = {
    "glnexus":         "#1f77b4",
    "parabricks":      "#ff7f0e",
    "gatk":            "#2ca02c",
    "gatk_genomicsdb": "#9467bd",
}

PLOTLY_CDN = "https://cdn.plot.ly/plotly-2.35.0.min.js"

# GQ windows: [0,10), [10,20), ..., [90,100]
GQ_WINDOWS = [f"{i*10}-{i*10+10}" for i in range(10)]

# Yield curve thresholds
GQ_THRESHOLDS = [0, 10, 20, 30, 40, 50]

# Transition pairs (as frozensets of upper-case nucleotides)
_TRANSITIONS = {frozenset({"A", "G"}), frozenset({"C", "T"})}

INDEL_SIZE_CAP = 30  # cap indel size tracking at ±N bp

# Quality thresholds for auto-interpretation
_TITV_NOISE  = 1.8   # Ti/Tv below this = likely noise
_TITV_WGS    = 2.0
_TITV_WES    = 3.0
_LOW_GQ_WARN = 0.50  # fraction of genotypes in [0-10] → warning


# ─── Helpers ──────────────────────────────────────────────────────────────────

def _mean(vals):
    return sum(vals) / len(vals) if vals else 0.0


def _std(vals):
    if len(vals) < 2:
        return 0.0
    m = _mean(vals)
    return math.sqrt(sum((v - m) ** 2 for v in vals) / (len(vals) - 1))


def _is_transition(ref: str, alt: str) -> bool:
    return frozenset({ref.upper(), alt.upper()}) in _TRANSITIONS


def _classify(ref: str, alt: str) -> str:
    """Return 'SNP', 'INDEL', or 'MULTI'."""
    if "," in alt:
        return "MULTI"
    if len(ref) == 1 and len(alt) == 1:
        return "SNP"
    if len(ref) != len(alt):
        return "INDEL"
    return "MULTI"  # MNP


def _gq_window_idx(gq: int) -> int:
    """Map GQ value to window index 0–9."""
    return min(gq // 10, 9)


# ─── bcftools wrappers ────────────────────────────────────────────────────────

def _check_bcftools() -> bool:
    try:
        r = subprocess.run(
            ["bcftools", "--version"], capture_output=True, text=True, check=True
        )
        print(f"[04c] {r.stdout.splitlines()[0]}")
        return True
    except (subprocess.CalledProcessError, FileNotFoundError):
        return False


def _run_bcftools_stats(vcf_path: str) -> dict:
    """
    Run `bcftools stats` and parse key summary lines.
    Returns dict with keys: snps, indels, multiallelic, mnps,
    ts, tv, titv, het_total, hom_alt_total.
    Fast operation (minutes, not hours).
    """
    result = {"snps": 0, "indels": 0, "multiallelic": 0, "mnps": 0,
              "ts": 0, "tv": 0, "titv": None,
              "het_total": 0, "hom_alt_total": 0}
    try:
        proc = subprocess.run(
            ["bcftools", "stats", str(vcf_path)],
            capture_output=True, text=True, timeout=1800
        )
        if proc.returncode != 0:
            print(f"[WARN] bcftools stats failed: {proc.stderr[:200].strip()}", file=sys.stderr)
            return result
    except subprocess.TimeoutExpired:
        print("[WARN] bcftools stats timed out", file=sys.stderr)
        return result

    for line in proc.stdout.splitlines():
        if line.startswith("SN"):
            parts = line.split("\t")
            if len(parts) < 4:
                continue
            key, val_str = parts[2].strip().rstrip(":"), parts[3].strip()
            try:
                val = int(val_str)
            except ValueError:
                continue
            if "SNPs" in key or "snps" in key.lower():
                result["snps"] = val
            elif "indels" in key.lower():
                result["indels"] = val
            elif "MNPs" in key or "mnps" in key.lower():
                result["mnps"] = val
            elif "multiallelic" in key.lower():
                result["multiallelic"] = val
        elif line.startswith("TSTV"):
            parts = line.split("\t")
            # TSTV  0  ts  tv  ts/tv  ts(1st ALT)  tv(1st ALT)  ts/tv(1st ALT)
            if len(parts) >= 5:
                try:
                    result["ts"]   = int(parts[2])
                    result["tv"]   = int(parts[3])
                    result["titv"] = round(float(parts[4]), 4)
                except (ValueError, IndexError):
                    pass
        elif line.startswith("PSC"):
            # PSC  id  sample  hom_RR  het  hom_AA  ...
            parts = line.split("\t")
            if len(parts) >= 6:
                try:
                    result["hom_alt_total"] += int(parts[4])  # hom_RR (ref-ref) at idx 3
                    result["het_total"]     += int(parts[4])
                    # Correct parsing: PSC col4=hom_RR, col5=het, col6=hom_AA
                    # Re-parse carefully
                except (ValueError, IndexError):
                    pass

    # Re-parse PSC correctly: col indices 3=hom_RR, 4=het, 5=hom_AA
    result["het_total"] = 0
    result["hom_alt_total"] = 0
    for line in proc.stdout.splitlines():
        if line.startswith("PSC"):
            parts = line.split("\t")
            if len(parts) >= 6:
                try:
                    result["het_total"]     += int(parts[4])  # het
                    result["hom_alt_total"] += int(parts[5])  # hom_AA
                except (ValueError, IndexError):
                    pass

    return result


# ─── Core streaming analysis ──────────────────────────────────────────────────

def _stream_gq_analysis(vcf_path: str) -> dict | None:
    """
    Stream `bcftools query -f "%REF\\t%ALT\\t[%GQ,]\\n"` line by line.

    Returns dict with:
      gq_windows  – {label: {snps, indels, transitions, transversions, titv}}
      yield_curve – {str(thresh): {snps, indels}}
      indel_size_dist – {str(delta): count}  (per variant SITE, not genotype)
      total_snp_sites, total_indel_sites, total_multi_sites  (per-site counts)
      total_snp_gts, total_indel_gts  (per-genotype GQ calls counted)
    """
    # List-based counters indexed by window index 0–9 (faster than dict access in hot loop)
    N = 10
    snp_gq  = [0] * N
    ind_gq  = [0] * N
    ti_gq   = [0] * N
    tv_gq   = [0] * N

    total_snp_sites  = 0
    total_indel_sites = 0
    total_multi_sites = 0
    indel_size_dist: dict[int, int] = {}

    # Yield accumulators — count per-GQ-call surviving each threshold
    yield_snp  = [0] * len(GQ_THRESHOLDS)
    yield_ind  = [0] * len(GQ_THRESHOLDS)
    # Precompute minimum threshold index for each window (window i covers GQ [i*10, i*10+10))
    # A call with window index w satisfies threshold t if t <= w*10
    thresh_min_window = [t // 10 for t in GQ_THRESHOLDS]  # window idx must be >= this

    cmd = ["bcftools", "query", "-f", r"%REF\t%ALT\t[%GQ,]\n", str(vcf_path)]
    try:
        proc = subprocess.Popen(
            cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            text=True, bufsize=1 << 20
        )
    except FileNotFoundError:
        print("[ERROR] bcftools not found", file=sys.stderr)
        return None

    sites = 0
    try:
        for raw_line in proc.stdout:
            line = raw_line.rstrip("\n")
            tab1 = line.find("\t")
            if tab1 < 0:
                continue
            tab2 = line.find("\t", tab1 + 1)
            if tab2 < 0:
                continue

            ref      = line[:tab1]
            alt      = line[tab1 + 1:tab2]
            gq_field = line[tab2 + 1:]

            vtype = _classify(ref, alt)

            if vtype == "MULTI":
                total_multi_sites += 1
                continue

            is_snp = (vtype == "SNP")
            if is_snp:
                total_snp_sites += 1
                is_ti = _is_transition(ref, alt)
            else:
                total_indel_sites += 1
                delta = len(alt) - len(ref)
                delta = max(-INDEL_SIZE_CAP, min(INDEL_SIZE_CAP, delta))
                indel_size_dist[delta] = indel_size_dist.get(delta, 0) + 1

            # Parse comma-separated GQ values
            for gq_str in gq_field.split(","):
                if not gq_str or gq_str == ".":
                    continue
                try:
                    gq = int(gq_str)
                except ValueError:
                    continue
                if gq < 0:
                    continue
                w = _gq_window_idx(gq)

                if is_snp:
                    snp_gq[w] += 1
                    if is_ti:
                        ti_gq[w] += 1
                    else:
                        tv_gq[w] += 1
                else:
                    ind_gq[w] += 1

                # Yield: for each threshold, accumulate if GQ >= threshold
                # (= window index >= thresh_min_window[t_idx])
                for t_idx, min_w in enumerate(thresh_min_window):
                    if w >= min_w:
                        if is_snp:
                            yield_snp[t_idx] += 1
                        else:
                            yield_ind[t_idx] += 1

            sites += 1
            if sites % 1_000_000 == 0:
                print(f"[04c]   … {sites:,} sites processed", flush=True)

    finally:
        proc.stdout.close()

    proc.wait()
    if proc.returncode not in (0, None):
        stderr_msg = proc.stderr.read(300).strip()
        print(f"[WARN] bcftools query rc={proc.returncode}: {stderr_msg}", file=sys.stderr)

    # Build structured output
    gq_windows = {}
    for i, label in enumerate(GQ_WINDOWS):
        tv = tv_gq[i]
        ti = ti_gq[i]
        gq_windows[label] = {
            "snps":          snp_gq[i],
            "indels":        ind_gq[i],
            "transitions":   ti,
            "transversions": tv,
            "titv":          round(ti / tv, 4) if tv > 0 else None,
        }

    yield_curve = {
        str(GQ_THRESHOLDS[i]): {"snps": yield_snp[i], "indels": yield_ind[i]}
        for i in range(len(GQ_THRESHOLDS))
    }

    return {
        "gq_windows":        gq_windows,
        "yield_curve":       yield_curve,
        "indel_size_dist":   {str(k): v for k, v in sorted(indel_size_dist.items())},
        "total_snp_sites":   total_snp_sites,
        "total_indel_sites": total_indel_sites,
        "total_multi_sites": total_multi_sites,
        "total_snp_gts":     sum(snp_gq),
        "total_indel_gts":   sum(ind_gq),
    }


# ─── Per-run orchestration ────────────────────────────────────────────────────

def _load_run_metrics(metrics_dir: Path, software: str, size: int, rep: int | None) -> dict | None:
    """Load metrics JSON; handles new (rep suffix) and legacy (no suffix) formats."""
    if rep is not None:
        p = metrics_dir / f"metrics_{software}_{size}_rep{rep}.json"
    else:
        p = metrics_dir / f"metrics_{software}_{size}.json"
    if not p.exists():
        return None
    try:
        with open(p) as f:
            return json.load(f)
    except Exception as e:
        print(f"[WARN] Could not read {p}: {e}", file=sys.stderr)
        return None


def _out_path(metrics_dir: Path, software: str, size: int, rep: int | None) -> Path:
    if rep is not None:
        return metrics_dir / f"variant_types_{software}_{size}_rep{rep}.json"
    return metrics_dir / f"variant_types_{software}_{size}.json"


def analyze_one_run(
    metrics_dir: Path,
    software: str,
    size: int,
    rep: int | None,
    force: bool,
) -> dict | None:
    out  = _out_path(metrics_dir, software, size, rep)
    tag  = f"rep{rep}" if rep is not None else "legacy"
    full = f"{software} N={size} {tag}"

    if out.exists() and not force:
        print(f"[04c] {full}: already done — skip")
        try:
            with open(out) as f:
                return json.load(f)
        except Exception:
            pass

    rec = _load_run_metrics(metrics_dir, software, size, rep)
    if rec is None:
        return None
    if rec.get("status") != "success":
        return None
    if not rec.get("output_valid", False):
        return None

    vcf_path = rec.get("output_vcf", "")
    if not vcf_path or not os.path.isfile(vcf_path):
        print(f"[04c] {full}: VCF not found ({vcf_path}) — skip", file=sys.stderr)
        return None

    print(f"[04c] {full}: running bcftools stats …")
    stats = _run_bcftools_stats(vcf_path)

    print(f"[04c] {full}: streaming GQ analysis …")
    gq_data = _stream_gq_analysis(vcf_path)
    if gq_data is None:
        return None

    het  = stats.get("het_total", 0)
    hom  = stats.get("hom_alt_total", 0)

    data = {
        "software":       software,
        "dataset_size":   size,
        "replicate":      rep if rep is not None else 1,
        "vcf_path":       vcf_path,
        # Per-site totals (from bcftools stats — authoritative)
        "total_snps":         stats.get("snps", gq_data["total_snp_sites"]),
        "total_indels":       stats.get("indels", gq_data["total_indel_sites"]),
        "total_multiallelic": stats.get("multiallelic", gq_data["total_multi_sites"]),
        "titv_overall":       stats.get("titv"),
        # Het/hom (from bcftools stats PSC lines, summed across samples)
        "het_total":    het,
        "hom_alt_total": hom,
        "het_hom_ratio": round(het / hom, 4) if hom > 0 else None,
        # GQ-stratified data (from streaming query)
        "gq_windows":      gq_data["gq_windows"],
        "yield_curve":     gq_data["yield_curve"],
        "indel_size_dist": gq_data["indel_size_dist"],
    }

    with open(out, "w") as f:
        json.dump(data, f, indent=2)
    print(f"[04c] Written: {out}")
    return data


def discover_runs(metrics_dir: Path) -> list[tuple[str, int, int | None]]:
    """Scan metrics_dir for metrics_*.json; return (software, size, rep|None) tuples."""
    runs = []
    for path in sorted(metrics_dir.glob("metrics_*.json")):
        m = re.match(r"metrics_(\w+)_(\d+)_rep(\d+)\.json", path.name)
        if m:
            runs.append((m.group(1), int(m.group(2)), int(m.group(3))))
            continue
        m = re.match(r"metrics_(\w+)_(\d+)\.json", path.name)
        if m:
            runs.append((m.group(1), int(m.group(2)), None))
    return runs


# ─── Data loading ─────────────────────────────────────────────────────────────

def load_type_data(metrics_dir: Path) -> dict:
    """Load all variant_types_*.json → {software: {size: [rec, ...]}}."""
    data: dict[str, dict[int, list]] = {s: {} for s in SOFTWARES}
    for path in sorted(metrics_dir.glob("variant_types_*.json")):
        m = re.match(r"variant_types_(\w+)_(\d+)_rep(\d+)\.json", path.name)
        if m:
            sw, size, rep = m.group(1), int(m.group(2)), int(m.group(3))
        else:
            m = re.match(r"variant_types_(\w+)_(\d+)\.json", path.name)
            if not m:
                continue
            sw, size, rep = m.group(1), int(m.group(2)), 1
        if sw not in data:
            data[sw] = {}
        try:
            with open(path) as f:
                rec = json.load(f)
        except Exception as e:
            print(f"[WARN] Could not read {path}: {e}", file=sys.stderr)
            continue
        data[sw].setdefault(size, []).append(rec)
    for sw in data:
        for size in data[sw]:
            data[sw][size].sort(key=lambda r: r.get("replicate", 1))
    return data


# ─── CSS / JS ─────────────────────────────────────────────────────────────────

_CSS = """
* { box-sizing: border-box; margin: 0; padding: 0; }
body { font-family: 'Segoe UI', Arial, sans-serif; background: #f5f6fa; color: #333; }
header { background: linear-gradient(135deg,#1A1A2E,#16213E); color: #fff; padding: 24px 32px; }
header h1 { font-size: 1.6rem; margin-bottom: 6px; }
header p  { opacity: .7; font-size: .9rem; }
.tabs { display: flex; gap: 4px; padding: 12px 32px 0; background: #fff;
        border-bottom: 2px solid #e0e0e0; flex-wrap: wrap; }
.tab-btn { padding: 8px 18px; border: none; border-radius: 6px 6px 0 0;
           cursor: pointer; background: #f0f0f0; font-size: .9rem; }
.tab-btn.active { background: #1A1A2E; color: #fff; }
.tab-content { display: none; padding: 24px 32px; }
.tab-content.active { display: block; }
.card { background: #fff; border-radius: 8px; padding: 20px;
        box-shadow: 0 1px 4px rgba(0,0,0,.08); margin-bottom: 16px; }
.section-title { font-size: 1.1rem; font-weight: 600; color: #1A1A2E;
                 margin: 24px 0 12px; border-left: 4px solid #4a90e2; padding-left: 10px; }
table { border-collapse: collapse; width: 100%; font-size: .875rem; }
th, td { padding: 8px 12px; border: 1px solid #e0e0e0; text-align: center; }
th { background: #1A1A2E; color: #fff; font-weight: 600; }
tr:nth-child(even) { background: #f9f9f9; }
.alert-box { border-radius: 6px; padding: 12px 18px; margin-bottom: 8px; border-left: 5px solid; }
.aw  { background: #fff8e1; border-color: #ff9800; }
.ac  { background: #fff0f0; border-color: #d62728; }
.ai  { background: #e8f5e9; border-color: #4caf50; }
.grid-2 { display: grid; grid-template-columns: 1fr 1fr; gap: 16px; }
.grid-3 { display: grid; grid-template-columns: 1fr 1fr 1fr; gap: 16px; }
.chart-wrap { background: #fff; border-radius: 8px; padding: 16px;
              box-shadow: 0 1px 4px rgba(0,0,0,.08); margin-bottom: 16px; }
select.sz-sel { padding: 6px 12px; border: 1px solid #ccc; border-radius: 4px;
                font-size: .875rem; cursor: pointer; }
.sel-bar { display: flex; align-items: center; gap: 12px; margin-bottom: 16px; }
@media (max-width:900px) { .grid-2,.grid-3 { grid-template-columns:1fr; } }
"""

_JS_TAB = """
function showTab(id, btn) {
  document.querySelectorAll('.tab-content').forEach(function(t){ t.classList.remove('active'); });
  document.querySelectorAll('.tab-btn').forEach(function(b){ b.classList.remove('active'); });
  document.getElementById(id).classList.add('active');
  btn.classList.add('active');
}
"""

_LAYOUT_BASE = json.dumps({
    "paper_bgcolor": "rgba(0,0,0,0)",
    "plot_bgcolor":  "#fafafa",
    "font": {"family": "Segoe UI, Arial", "size": 12},
    "legend": {"orientation": "h", "y": -0.28},
    "hovermode": "x unified",
    "margin": {"l": 70, "r": 20, "t": 48, "b": 80},
})


def _j(obj) -> str:
    return json.dumps(obj)


# ─── Interpretation ───────────────────────────────────────────────────────────

def _build_interpretation(data: dict) -> str:
    cards = []
    for sw in SOFTWARES:
        for size, recs in data.get(sw, {}).items():
            lbl = f"{SOFT_LABEL.get(sw, sw)} N={size}"
            for rec in recs:
                rep = rec.get("replicate", "?")
                tag = f"{lbl} rep{rep}"

                # Fraction of genotype calls in GQ [0-10]
                total_gts = rec.get("total_snp_gts", 0) + rec.get("total_indel_gts", 0)
                low_gts   = (rec.get("gq_windows", {}).get("0-10", {}).get("snps", 0) +
                             rec.get("gq_windows", {}).get("0-10", {}).get("indels", 0))
                if total_gts > 0 and low_gts / total_gts > _LOW_GQ_WARN:
                    pct = 100 * low_gts / total_gts
                    cards.append(
                        f'<div class="alert-box aw"><b>ATENÇÃO — {tag}</b>: '
                        f'{pct:.0f}% dos genótipos têm GQ 0–10 '
                        f'(possível excesso de chamadas de baixa qualidade).</div>'
                    )

                # Ti/Tv in GQ [0-10]
                w0  = rec.get("gq_windows", {}).get("0-10", {})
                ti0 = w0.get("transitions", 0)
                tv0 = w0.get("transversions", 0)
                if tv0 > 0 and ti0 / tv0 < _TITV_NOISE:
                    cards.append(
                        f'<div class="alert-box ac"><b>CRÍTICO — {tag}</b>: '
                        f'Ti/Tv = {ti0/tv0:.2f} na janela GQ 0–10 '
                        f'(abaixo de {_TITV_NOISE} → ruído; esperado ≥ {_TITV_WGS} WGS / ≥ {_TITV_WES} WES).</div>'
                    )

                # HQ fraction info
                total_snps = rec.get("total_snps", 0)
                hq_snps    = rec.get("yield_curve", {}).get("20", {}).get("snps", 0)
                if total_snps > 0:
                    pct_hq = 100 * hq_snps / total_snps
                    cards.append(
                        f'<div class="alert-box ai"><b>INFO — {tag}</b>: '
                        f'{pct_hq:.1f}% dos SNPs têm GQ ≥ 20 ({hq_snps:,} de {total_snps:,}).</div>'
                    )

    if not cards:
        return '<div class="alert-box ai">Sem dados suficientes para interpretação automática.</div>'
    return "\n".join(cards)


# ─── Tab 1: Visão Geral ───────────────────────────────────────────────────────

def _build_overview(data: dict) -> tuple[str, str]:
    all_sizes = sorted({s for sd in data.values() for s in sd})

    traces_snp, traces_ind = [], []
    for sw in SOFTWARES:
        xs_s, ys_s, es_s = [], [], []
        xs_i, ys_i, es_i = [], [], []
        for size in all_sizes:
            recs = data.get(sw, {}).get(size, [])
            sv = [r["total_snps"]   for r in recs if "total_snps"   in r]
            iv = [r["total_indels"] for r in recs if "total_indels" in r]
            if sv:
                xs_s.append(size); ys_s.append(round(_mean(sv))); es_s.append(round(_std(sv)))
            if iv:
                xs_i.append(size); ys_i.append(round(_mean(iv))); es_i.append(round(_std(iv)))

        c = SW_COLOR.get(sw, "#888")
        l = SOFT_LABEL.get(sw, sw)
        if xs_s:
            t = {"type": "bar", "name": l, "x": xs_s, "y": ys_s,
                 "marker": {"color": c},
                 "hovertemplate": l + " N=%{x}: %{y:,} SNPs<extra></extra>"}
            if any(e > 0 for e in es_s):
                t["error_y"] = {"type": "data", "array": es_s, "visible": True}
            traces_snp.append(t)
        if xs_i:
            t = {"type": "bar", "name": l, "x": xs_i, "y": ys_i,
                 "marker": {"color": c, "opacity": 0.75},
                 "hovertemplate": l + " N=%{x}: %{y:,} INDELs<extra></extra>"}
            if any(e > 0 for e in es_i):
                t["error_y"] = {"type": "data", "array": es_i, "visible": True}
            traces_ind.append(t)

    cat_axis = {"title": "Tamanho do dataset (# GVCFs)", "type": "category", "gridcolor": "#eee"}
    layout_snp = {"barmode": "group", "title": "Total de SNPs por Ferramenta",
                  "xaxis": cat_axis, "yaxis": {"title": "SNPs", "gridcolor": "#eee"}}
    layout_ind = {"barmode": "group", "title": "Total de INDELs por Ferramenta",
                  "xaxis": cat_axis, "yaxis": {"title": "INDELs", "gridcolor": "#eee"}}

    js = (
        f"Plotly.newPlot('ov-snp',{_j(traces_snp)},Object.assign({{}},{_LAYOUT_BASE},{_j(layout_snp)}),{{responsive:true}});"
        f"Plotly.newPlot('ov-ind',{_j(traces_ind)},Object.assign({{}},{_LAYOUT_BASE},{_j(layout_ind)}),{{responsive:true}});"
    )

    # SNP:Indel ratio table
    hdr = "".join(f"<th>N={s}</th>" for s in all_sizes)
    rows = []
    for sw in SOFTWARES:
        cells = [f'<td style="font-weight:600;text-align:left">{SOFT_LABEL.get(sw,sw)}</td>']
        for size in all_sizes:
            recs = data.get(sw, {}).get(size, [])
            sm = _mean([r["total_snps"]   for r in recs if "total_snps"   in r])
            im = _mean([r["total_indels"] for r in recs if "total_indels" in r])
            cells.append(f"<td>{sm/im:.1f}</td>" if im > 0 else "<td>—</td>")
        rows.append("<tr>" + "".join(cells) + "</tr>")

    ratio_tbl = (
        f'<table><thead><tr><th style="text-align:left">Ferramenta</th>{hdr}</tr></thead>'
        f'<tbody>{"".join(rows)}</tbody></table>'
    )

    html = (
        '<div class="grid-2">'
        '<div class="chart-wrap"><div id="ov-snp" style="height:380px"></div></div>'
        '<div class="chart-wrap"><div id="ov-ind" style="height:380px"></div></div>'
        '</div>'
        '<p class="section-title">Ratio SNP:Indel — esperado 4–8× para WGS humano (média entre réplicas)</p>'
        f'<div class="card">{ratio_tbl}</div>'
    )
    return html, js


# ─── Tab 2: Qualidade por Tipo ────────────────────────────────────────────────

def _build_quality(data: dict) -> tuple[str, str]:
    all_sizes = sorted({s for sd in data.values() for s in sd})
    if not all_sizes:
        return "<p>Sem dados.</p>", ""

    # Embed all data as JS dict
    qdata: dict = {}
    for size in all_sizes:
        qdata[size] = {}
        for sw in SOFTWARES:
            recs = data.get(sw, {}).get(size, [])
            if not recs:
                continue
            avg: dict = {}
            for w in GQ_WINDOWS:
                sv = [r.get("gq_windows", {}).get(w, {}).get("snps",   0) for r in recs]
                iv = [r.get("gq_windows", {}).get(w, {}).get("indels", 0) for r in recs]
                avg[w] = {"snps": round(_mean(sv)), "indels": round(_mean(iv))}
            qdata[size][sw] = avg

    opts = "".join(f'<option value="{s}">N = {s}</option>' for s in all_sizes)
    html = (
        '<div class="sel-bar"><label><b>Dataset size:</b></label>'
        f'<select class="sz-sel" onchange="updateQuality(this.value)">{opts}</select></div>'
        '<div class="grid-2">'
        '<div class="chart-wrap"><div id="q-snp" style="height:380px"></div></div>'
        '<div class="chart-wrap"><div id="q-ind" style="height:380px"></div></div>'
        '</div>'
    )
    js = f"""
var qData={_j(qdata)};
function updateQuality(sz){{
  sz=parseInt(sz);
  var d=qData[sz]||{{}};
  var wins={_j(GQ_WINDOWS)};
  var cols={_j(SW_COLOR)};
  var lbls={_j(SOFT_LABEL)};
  var sws={_j(SOFTWARES)};
  var ts=[],ti=[];
  sws.forEach(function(sw){{
    if(!d[sw])return;
    var ys=wins.map(function(w){{return d[sw][w]?d[sw][w].snps:0;}});
    var yi=wins.map(function(w){{return d[sw][w]?d[sw][w].indels:0;}});
    ts.push({{x:wins,y:ys,name:lbls[sw]||sw,type:'bar',marker:{{color:cols[sw]||'#888'}}}});
    ti.push({{x:wins,y:yi,name:lbls[sw]||sw,type:'bar',marker:{{color:cols[sw]||'#888'}}}});
  }});
  var base={{barmode:'group',xaxis:{{title:'Janela de GQ'}},yaxis:{{title:'Genótipos',gridcolor:'#eee'}},legend:{{orientation:'h',y:-0.3}},paper_bgcolor:'rgba(0,0,0,0)',plot_bgcolor:'#fafafa',margin:{{l:70,r:20,t:48,b:90}}}};
  Plotly.react('q-snp',ts,Object.assign({{}},base,{{title:'SNPs por Janela de GQ — N='+sz}}),{{responsive:true}});
  Plotly.react('q-ind',ti,Object.assign({{}},base,{{title:'INDELs por Janela de GQ — N='+sz}}),{{responsive:true}});
}}
updateQuality({all_sizes[-1]});
"""
    return html, js


# ─── Tab 3: Ti/Tv por Qualidade ───────────────────────────────────────────────

def _build_titv(data: dict) -> tuple[str, str]:
    all_sizes = sorted({s for sd in data.values() for s in sd})
    if not all_sizes:
        return "<p>Sem dados.</p>", ""

    tvdata: dict = {}
    for size in all_sizes:
        tvdata[size] = {}
        for sw in SOFTWARES:
            recs = data.get(sw, {}).get(size, [])
            if not recs:
                continue
            row = {}
            for w in GQ_WINDOWS:
                vals = [r.get("gq_windows", {}).get(w, {}).get("titv")
                        for r in recs
                        if r.get("gq_windows", {}).get(w, {}).get("titv") is not None]
                row[w] = round(_mean(vals), 4) if vals else None
            tvdata[size][sw] = row

    opts = "".join(f'<option value="{s}">N = {s}</option>' for s in all_sizes)
    html = (
        '<div class="card"><p>Ti/Tv próximo de 1,0 indica variantes aleatórias (ruído). '
        'Esperado ≥ 2,0 para WGS e ≥ 3,0 para WES. '
        'Linha de Ruído (vermelho) = 1,8. '
        'Uma curva crescente da esquerda para a direita confirma que '
        'variantes de GQ alto são mais confiáveis.</p></div>'
        '<div class="sel-bar"><label><b>Dataset size:</b></label>'
        f'<select class="sz-sel" onchange="updateTiTv(this.value)">{opts}</select></div>'
        '<div class="chart-wrap"><div id="tv-chart" style="height:450px"></div></div>'
    )
    js = f"""
var tvData={_j(tvdata)};
function updateTiTv(sz){{
  sz=parseInt(sz);
  var d=tvData[sz]||{{}};
  var wins={_j(GQ_WINDOWS)};
  var cols={_j(SW_COLOR)};
  var lbls={_j(SOFT_LABEL)};
  var sws={_j(SOFTWARES)};
  var traces=[];
  sws.forEach(function(sw){{
    if(!d[sw])return;
    var ys=wins.map(function(w){{return d[sw][w];}});
    if(ys.every(function(v){{return v===null;}}))return;
    traces.push({{x:wins,y:ys,name:lbls[sw]||sw,mode:'lines+markers',connectgaps:true,
      line:{{color:cols[sw]||'#888',width:2.5}},marker:{{size:8}},
      hovertemplate:(lbls[sw]||sw)+' %{{x}}: Ti/Tv=%{{y:.3f}}<extra></extra>'}});
  }});
  var layout={{
    title:'Ti/Tv por Janela de GQ — N='+sz,
    xaxis:{{title:'Janela de GQ'}},
    yaxis:{{title:'Ti/Tv ratio',range:[0,4.5],gridcolor:'#eee'}},
    legend:{{orientation:'h',y:-0.28}},
    hovermode:'x unified',
    paper_bgcolor:'rgba(0,0,0,0)',plot_bgcolor:'#fafafa',
    margin:{{l:70,r:20,t:48,b:80}},
    shapes:[
      {{type:'line',x0:wins[0],x1:wins[wins.length-1],y0:2.0,y1:2.0,line:{{color:'#607D8B',dash:'dash',width:1.5}}}},
      {{type:'line',x0:wins[0],x1:wins[wins.length-1],y0:3.0,y1:3.0,line:{{color:'#795548',dash:'dot',width:1.5}}}},
      {{type:'line',x0:wins[0],x1:wins[wins.length-1],y0:1.8,y1:1.8,line:{{color:'#d62728',dash:'dash',width:1.0}}}}
    ],
    annotations:[
      {{x:wins[wins.length-1],y:2.08,xanchor:'right',text:'WGS (2.0)',showarrow:false,font:{{color:'#607D8B',size:10}}}},
      {{x:wins[wins.length-1],y:3.08,xanchor:'right',text:'WES (3.0)',showarrow:false,font:{{color:'#795548',size:10}}}},
      {{x:wins[wins.length-1],y:1.72,xanchor:'right',text:'Ruído (<1.8)',showarrow:false,font:{{color:'#d62728',size:10}}}}
    ]
  }};
  Plotly.react('tv-chart',traces,layout,{{responsive:true}});
}}
updateTiTv({all_sizes[-1]});
"""
    return html, js


# ─── Tab 4: Curva de Rendimento ───────────────────────────────────────────────

def _build_yield(data: dict) -> tuple[str, str]:
    all_sizes = sorted({s for sd in data.values() for s in sd})
    if not all_sizes:
        return "<p>Sem dados.</p>", ""

    ydata: dict = {}
    for size in all_sizes:
        ydata[size] = {}
        for sw in SOFTWARES:
            recs = data.get(sw, {}).get(size, [])
            if not recs:
                continue
            row: dict = {}
            for t in GQ_THRESHOLDS:
                sv = [r.get("yield_curve", {}).get(str(t), {}).get("snps",   0) for r in recs]
                iv = [r.get("yield_curve", {}).get(str(t), {}).get("indels", 0) for r in recs]
                row[t] = {"snps": round(_mean(sv)), "indels": round(_mean(iv))}
            ydata[size][sw] = row

    opts = "".join(f'<option value="{s}">N = {s}</option>' for s in all_sizes)
    html = (
        '<div class="card"><p>Variantes retidas ao aplicar um filtro mínimo de GQ. '
        'Queda abrupta = muitas variantes de baixa qualidade. '
        'Curva mais plana indica variantes concentradas acima do threshold.</p></div>'
        '<div class="sel-bar"><label><b>Dataset size:</b></label>'
        f'<select class="sz-sel" onchange="updateYield(this.value)">{opts}</select></div>'
        '<div class="grid-2">'
        '<div class="chart-wrap"><div id="yd-snp" style="height:380px"></div></div>'
        '<div class="chart-wrap"><div id="yd-ind" style="height:380px"></div></div>'
        '</div>'
    )
    js = f"""
var ydData={_j(ydata)};
function updateYield(sz){{
  sz=parseInt(sz);
  var d=ydData[sz]||{{}};
  var thrs={_j(GQ_THRESHOLDS)};
  var cols={_j(SW_COLOR)};
  var lbls={_j(SOFT_LABEL)};
  var sws={_j(SOFTWARES)};
  var ts=[],ti=[];
  sws.forEach(function(sw){{
    if(!d[sw])return;
    var ys=thrs.map(function(t){{return d[sw][t]?d[sw][t].snps:0;}});
    var yi=thrs.map(function(t){{return d[sw][t]?d[sw][t].indels:0;}});
    ts.push({{x:thrs,y:ys,name:lbls[sw]||sw,mode:'lines+markers',
      line:{{color:cols[sw]||'#888',width:2}},marker:{{size:7}},
      hovertemplate:(lbls[sw]||sw)+' GQ≥%{{x}}: %{{y:,}}<extra></extra>'}});
    ti.push({{x:thrs,y:yi,name:lbls[sw]||sw,mode:'lines+markers',
      line:{{color:cols[sw]||'#888',width:2}},marker:{{size:7}},
      hovertemplate:(lbls[sw]||sw)+' GQ≥%{{x}}: %{{y:,}}<extra></extra>'}});
  }});
  var base={{xaxis:{{title:'Threshold mínimo de GQ'}},yaxis:{{title:'Genótipos retidos',gridcolor:'#eee'}},
    legend:{{orientation:'h',y:-0.28}},hovermode:'x unified',
    paper_bgcolor:'rgba(0,0,0,0)',plot_bgcolor:'#fafafa',margin:{{l:70,r:20,t:48,b:80}}}};
  Plotly.react('yd-snp',ts,Object.assign({{}},base,{{title:'SNPs retidos — N='+sz}}),{{responsive:true}});
  Plotly.react('yd-ind',ti,Object.assign({{}},base,{{title:'INDELs retidos — N='+sz}}),{{responsive:true}});
}}
updateYield({all_sizes[-1]});
"""
    return html, js


# ─── Tab 5: Distribuição de Indels ────────────────────────────────────────────

def _build_indel_dist(data: dict) -> tuple[str, str]:
    all_sizes = sorted({s for sd in data.values() for s in sd})
    if not all_sizes:
        return "<p>Sem dados.</p>", ""

    idata: dict = {}
    for size in all_sizes:
        idata[size] = {}
        for sw in SOFTWARES:
            recs = data.get(sw, {}).get(size, [])
            if not recs:
                continue
            # Collect all deltas seen across reps
            all_deltas: set[int] = set()
            for r in recs:
                all_deltas |= {int(k) for k in r.get("indel_size_dist", {})}
            if not all_deltas:
                continue
            row: dict = {}
            for d in sorted(all_deltas):
                vals = [r.get("indel_size_dist", {}).get(str(d), 0) for r in recs]
                row[d] = round(_mean(vals))
            idata[size][sw] = row

    opts = "".join(f'<option value="{s}">N = {s}</option>' for s in all_sizes)
    html = (
        '<div class="card"><p>Distribuição de tamanho de indels '
        '(negativo = deleção, positivo = inserção). '
        'Picos anômalos em tamanhos específicos indicam artefatos sistemáticos. '
        'Deleções de 1 bp são as mais comuns em Illumina.</p></div>'
        '<div class="sel-bar"><label><b>Dataset size:</b></label>'
        f'<select class="sz-sel" onchange="updateIndel(this.value)">{opts}</select></div>'
        '<div class="chart-wrap"><div id="id-chart" style="height:480px"></div></div>'
    )
    js = f"""
var idData={_j(idata)};
function updateIndel(sz){{
  sz=parseInt(sz);
  var d=idData[sz]||{{}};
  var cols={_j(SW_COLOR)};
  var lbls={_j(SOFT_LABEL)};
  var sws={_j(SOFTWARES)};
  var traces=[];
  sws.forEach(function(sw){{
    if(!d[sw])return;
    var xs=Object.keys(d[sw]).map(Number).sort(function(a,b){{return a-b;}});
    var ys=xs.map(function(x){{return d[sw][x]||0;}});
    traces.push({{x:xs,y:ys,name:lbls[sw]||sw,type:'bar',
      marker:{{color:cols[sw]||'#888',opacity:0.75}},
      hovertemplate:(lbls[sw]||sw)+' Δ=%{{x}}: %{{y:,}}<extra></extra>'}});
  }});
  var layout={{
    barmode:'overlay',
    title:'Distribuição de Tamanho de Indels — N='+sz,
    xaxis:{{title:'Tamanho do Indel (bp)',gridcolor:'#eee'}},
    yaxis:{{title:'Contagem (média entre réplicas)',gridcolor:'#eee'}},
    legend:{{orientation:'h',y:-0.28}},
    paper_bgcolor:'rgba(0,0,0,0)',plot_bgcolor:'#fafafa',
    margin:{{l:70,r:20,t:48,b:90}},
    shapes:[{{type:'line',x0:0,x1:0,y0:0,y1:1,yref:'paper',
              line:{{color:'#aaa',dash:'dot',width:1}}}}]
  }};
  Plotly.react('id-chart',traces,layout,{{responsive:true}});
}}
updateIndel({all_sizes[-1]});
"""
    return html, js


# ─── Tab 6: Het/Hom ──────────────────────────────────────────────────────────

def _build_hetphom(data: dict) -> tuple[str, str]:
    all_sizes = sorted({s for sd in data.values() for s in sd})
    traces = []
    for sw in SOFTWARES:
        xs, ys, es = [], [], []
        for size in all_sizes:
            recs = data.get(sw, {}).get(size, [])
            vals = [r["het_hom_ratio"] for r in recs
                    if r.get("het_hom_ratio") is not None]
            if vals:
                xs.append(size); ys.append(round(_mean(vals), 3)); es.append(round(_std(vals), 3))
        if xs:
            t = {"type": "scatter", "mode": "lines+markers", "name": SOFT_LABEL.get(sw, sw),
                 "x": xs, "y": ys,
                 "line": {"color": SW_COLOR.get(sw, "#888"), "width": 2},
                 "marker": {"size": 8},
                 "hovertemplate": SOFT_LABEL.get(sw, sw) + " N=%{x}: Het/Hom=%{y:.3f}<extra></extra>"}
            if any(e > 0 for e in es):
                t["error_y"] = {"type": "data", "array": es, "visible": True}
            traces.append(t)

    layout = {
        "title": "Ratio Het/Hom por Ferramenta e Tamanho de Dataset",
        "xaxis": {"title": "Tamanho do dataset (# GVCFs)", "type": "category", "gridcolor": "#eee"},
        "yaxis": {"title": "Ratio Het/Hom", "gridcolor": "#eee"},
        "shapes": [{"type": "line", "x0": 0, "x1": 1, "xref": "paper",
                    "y0": 2.0, "y1": 2.0,
                    "line": {"color": "#607D8B", "dash": "dash", "width": 1.5}}],
        "annotations": [{"x": 0.98, "y": 2.06, "xref": "paper", "text": "Referência ~2.0",
                         "showarrow": False, "font": {"color": "#607D8B", "size": 10}}],
    }
    js = (f"Plotly.newPlot('hh-chart',{_j(traces)},"
          f"Object.assign({{}},{_LAYOUT_BASE},{_j(layout)}),{{responsive:true}});")

    html = (
        '<div class="card"><p>Ratio het/hom esperado ~2 para populações humanas '
        '(Hardy-Weinberg). Valores muito baixos (&lt;1,0) indicam super-homozigotização '
        '(possível contaminação ou chamadas incorretas). Valores muito altos indicam '
        'excesso de heterozigosidade (artefatos ou amostras relacionadas).</p></div>'
        '<div class="chart-wrap"><div id="hh-chart" style="height:420px"></div></div>'
    )
    return html, js


# ─── Tab 7: Dados Brutos ──────────────────────────────────────────────────────

def _build_rawdata(data: dict) -> str:
    rows = []
    for sw in SOFTWARES:
        for size in sorted(data.get(sw, {})):
            for rec in data[sw][size]:
                rep    = rec.get("replicate", "?")
                snps   = rec.get("total_snps", 0)
                indels = rec.get("total_indels", 0)
                multi  = rec.get("total_multiallelic", 0)
                titv   = rec.get("titv_overall")
                ratio  = rec.get("het_hom_ratio")
                hq_s   = rec.get("yield_curve", {}).get("20", {}).get("snps")
                hq_i   = rec.get("yield_curve", {}).get("20", {}).get("indels")
                td = lambda v, fmt=",": (f"{v:{fmt}}" if isinstance(v, (int, float)) else "—")
                rows.append(
                    f"<tr>"
                    f'<td style="text-align:left">{SOFT_LABEL.get(sw,sw)}</td>'
                    f"<td>{size}</td><td>{rep}</td>"
                    f"<td>{td(snps)}</td><td>{td(indels)}</td><td>{td(multi)}</td>"
                    f"<td>{f'{titv:.3f}' if titv else '—'}</td>"
                    f"<td>{f'{ratio:.2f}' if ratio else '—'}</td>"
                    f"<td>{td(hq_s)}</td><td>{td(hq_i)}</td>"
                    f"</tr>"
                )
    hdr = (
        '<tr><th style="text-align:left">Ferramenta</th>'
        "<th>N</th><th>Rep</th>"
        "<th>SNPs totais</th><th>INDELs totais</th><th>Multialélicos</th>"
        "<th>Ti/Tv global</th><th>Het/Hom</th>"
        "<th>SNPs GQ≥20</th><th>INDELs GQ≥20</th></tr>"
    )
    tbl = f"<table><thead>{hdr}</thead><tbody>{''.join(rows) or '<tr><td colspan=10>—</td></tr>'}</tbody></table>"
    return f'<div class="card">{tbl}</div>'


# ─── HTML report ──────────────────────────────────────────────────────────────

def generate_html(data: dict, output_path: Path) -> None:
    interp          = _build_interpretation(data)
    ov_h,   ov_js   = _build_overview(data)
    qlt_h,  qlt_js  = _build_quality(data)
    tv_h,   tv_js   = _build_titv(data)
    yld_h,  yld_js  = _build_yield(data)
    id_h,   id_js   = _build_indel_dist(data)
    hh_h,   hh_js   = _build_hetphom(data)
    raw_h           = _build_rawdata(data)

    all_js = "\n".join([ov_js, qlt_js, tv_js, yld_js, id_js, hh_js])

    html = f"""<!DOCTYPE html>
<html lang="pt-BR">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width,initial-scale=1.0">
  <title>GVCF Benchmark — Análise de Variantes por Tipo</title>
  <script src="{PLOTLY_CDN}"></script>
  <style>{_CSS}</style>
</head>
<body>

<header>
  <h1>Análise de Variantes: SNPs, INDELs e Qualidade GQ</h1>
  <p>Gerado: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')} &nbsp;|&nbsp;
     Stratificação por janelas de GQ [0–10, 10–20, …, 90–100]</p>
</header>

<div style="padding:20px 32px 0">
  <p class="section-title" style="margin-top:0">Interpretação Automática</p>
  {interp}
</div>

<div class="tabs">
  <button class="tab-btn active" onclick="showTab('t-ov',this)">Visão Geral</button>
  <button class="tab-btn" onclick="showTab('t-q',this)">Qualidade por Tipo</button>
  <button class="tab-btn" onclick="showTab('t-tv',this)">Ti/Tv por Qualidade ★</button>
  <button class="tab-btn" onclick="showTab('t-yd',this)">Curva de Rendimento</button>
  <button class="tab-btn" onclick="showTab('t-id',this)">Distribuição Indels</button>
  <button class="tab-btn" onclick="showTab('t-hh',this)">Het/Hom</button>
  <button class="tab-btn" onclick="showTab('t-raw',this)">Dados Brutos</button>
</div>

<div id="t-ov"  class="tab-content active">{ov_h}</div>
<div id="t-q"   class="tab-content">{qlt_h}</div>
<div id="t-tv"  class="tab-content">{tv_h}</div>
<div id="t-yd"  class="tab-content">{yld_h}</div>
<div id="t-id"  class="tab-content">{id_h}</div>
<div id="t-hh"  class="tab-content">{hh_h}</div>
<div id="t-raw" class="tab-content">{raw_h}</div>

<footer style="padding:16px 32px;text-align:center;color:#999;font-size:.8rem;margin-top:32px">
  GVCF Scalability Benchmark — Análise de Tipos de Variantes &mdash; {datetime.now().year}
</footer>

<script>
{_JS_TAB}
(function(){{
{all_js}
}})();
</script>

</body>
</html>"""

    output_path.write_text(html, encoding="utf-8")
    print(f"[04c] HTML report: {output_path}")


# ─── Main ─────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(
        description=(
            "Variant type & quality analysis: SNP/Indel stratified by GQ, "
            "Ti/Tv ratio, indel size distribution, het/hom ratio, GQ yield curves."
        )
    )
    ap.add_argument("--benchmark-dir", required=True,
                    help="Root benchmarks directory (contains 03_metrics/)")
    ap.add_argument("--output-dir",    required=True,
                    help="Directory for variant_analysis_report.html")
    ap.add_argument("--force", action="store_true",
                    help="Re-run analysis even if output JSON already exists")
    ap.add_argument("--skip-analysis", action="store_true",
                    help="Skip VCF analysis; only regenerate HTML from existing JSONs")
    args = ap.parse_args()

    metrics_dir = Path(args.benchmark_dir) / "03_metrics"
    output_dir  = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    if not metrics_dir.exists():
        print(f"[ERROR] Metrics directory not found: {metrics_dir}", file=sys.stderr)
        sys.exit(1)

    # ── Phase 1: VCF analysis
    if not args.skip_analysis:
        if not _check_bcftools():
            print("[ERROR] bcftools not found in PATH", file=sys.stderr)
            sys.exit(1)

        runs = discover_runs(metrics_dir)
        if not runs:
            print("[WARN] No metrics_*.json files found", file=sys.stderr)

        for sw, size, rep in runs:
            analyze_one_run(metrics_dir, sw, size, rep, args.force)

    # ── Phase 2: HTML report
    print("[04c] Loading results for report …")
    data = load_type_data(metrics_dir)
    total = sum(sum(len(v) for v in sd.values()) for sd in data.values())
    print(f"[04c] Loaded {total} variant-type records")

    if total == 0:
        print("[WARN] No variant_types_*.json found — report will be empty")

    html_path = output_dir / "variant_analysis_report.html"
    generate_html(data, html_path)
    print("[04c] Done.")


if __name__ == "__main__":
    main()
