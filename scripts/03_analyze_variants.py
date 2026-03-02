#!/usr/bin/env python3
"""
Step 3b — Variant Quality Analysis

Extracts per-run GQ histograms and inter-tool variant overlap statistics.
Must be run after the benchmark tools have produced VCF output.

Outputs to benchmarks/03_metrics/:
    variant_analysis_{software}_{size}.json  — GQ histogram (density)
    variant_overlap_{size}.json              — inter-tool overlap counts

Usage:
    python3 03_analyze_variants.py --benchmark-dir /path/to/benchmarks [--force]

Requirements:
    bcftools in PATH
"""

import argparse
import json
import os
import subprocess
import sys
from pathlib import Path

SOFTWARES = ["glnexus", "parabricks", "gatk"]
SIZES     = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
BIN_WIDTH = 5   # GQ histogram bins: [0,5), [5,10), … [95,100]


# ─── bcftools wrappers ────────────────────────────────────────────────────────

def _bcftools(args: list[str], vcf_path: str) -> list[str]:
    """Run bcftools <args> <vcf_path>; return stdout lines or [] on error."""
    cmd = ["bcftools"] + args + [str(vcf_path)]
    try:
        result = subprocess.run(
            cmd, capture_output=True, text=True, check=True
        )
        return result.stdout.splitlines()
    except subprocess.CalledProcessError as exc:
        print(
            f"[WARN] bcftools failed for {vcf_path}: "
            f"{exc.stderr[:300].strip()}",
            file=sys.stderr,
        )
        return []


def extract_gq_histogram(vcf_path: str) -> tuple[dict, int]:
    """
    Extract per-sample GQ values and return a density histogram.

    Returns:
        hist   – {str(bin_start): density}  e.g. {"0": 0.012, "5": 0.034, …}
        total  – number of non-missing genotype calls processed
    """
    lines = _bcftools(["query", "-f", "[%GQ\\n]"], vcf_path)

    counts: dict[int, int] = {b: 0 for b in range(0, 100, BIN_WIDTH)}
    total = 0
    for line in lines:
        v = line.strip()
        if not v or v == ".":
            continue
        try:
            gq = int(float(v))
        except ValueError:
            continue
        if gq < 0:
            continue
        b = min((gq // BIN_WIDTH) * BIN_WIDTH, 95)
        counts[b] = counts.get(b, 0) + 1
        total += 1

    if total == 0:
        return {str(k): 0.0 for k in sorted(counts)}, 0

    hist = {str(k): round(v / total, 6) for k, v in sorted(counts.items())}
    return hist, total


def extract_variant_ids(vcf_path: str) -> set[str]:
    """Return set of CHROM:POS:REF:ALT identifiers (in-memory only)."""
    lines = _bcftools(["query", "-f", "%CHROM:%POS:%REF:%ALT\\n"], vcf_path)
    return {line.strip() for line in lines if line.strip()}


# ─── Overlap computation ──────────────────────────────────────────────────────

def compute_overlap(id_sets: dict[str, set]) -> dict:
    """
    Compute exclusive, pairwise and all-common variant counts.

    id_sets: {software: set_of_variant_ids}  (only tools with data)
    """
    tools = sorted(id_sets.keys())
    all_ids = set().union(*id_sets.values()) if id_sets else set()

    result: dict = {
        "tools":         tools,
        "total_unique":  len(all_ids),
        "counts":        {t: len(id_sets[t]) for t in tools},
        "exclusive":     {},
        "pairwise_only": {},
        "all_common":    0,
    }

    # Per-tool exclusive (in this tool only)
    for t in tools:
        others = (
            set().union(*(id_sets[o] for o in tools if o != t))
            if len(tools) > 1 else set()
        )
        result["exclusive"][t] = len(id_sets[t] - others)

    # All-common
    if tools:
        common = id_sets[tools[0]].copy()
        for t in tools[1:]:
            common &= id_sets[t]
        result["all_common"] = len(common)

    # Pairwise intersections
    for i in range(len(tools)):
        for j in range(i + 1, len(tools)):
            t1, t2 = tools[i], tools[j]
            in_both = id_sets[t1] & id_sets[t2]
            others  = [t for t in tools if t not in (t1, t2)]
            others_union = (
                set().union(*(id_sets[t] for t in others))
                if others else set()
            )
            result["pairwise_only"][f"{t1}_{t2}"] = len(in_both - others_union)
            result[f"count_{t1}_{t2}"]            = len(in_both)

    return result


# ─── Per-size runners ─────────────────────────────────────────────────────────

def _load_run_metrics(metrics_dir: Path, software: str, size: int) -> dict | None:
    p = metrics_dir / f"metrics_{software}_{size}.json"
    if not p.exists():
        return None
    try:
        with open(p) as f:
            return json.load(f)
    except Exception:
        return None


def analyze_gq(
    metrics_dir: Path, software: str, size: int, force: bool
) -> dict | None:
    """Extract GQ histogram for one run. Returns the data dict or None."""
    out = metrics_dir / f"variant_analysis_{software}_{size}.json"

    if out.exists() and not force:
        print(f"[ANALYZE] {software} N={size}: already done — skip")
        try:
            with open(out) as f:
                return json.load(f)
        except Exception:
            pass

    rec = _load_run_metrics(metrics_dir, software, size)
    if rec is None:
        return None
    if rec.get("status") != "success":
        return None

    vcf_path = rec.get("output_vcf", "")
    if not vcf_path or not os.path.isfile(vcf_path):
        print(f"[ANALYZE] {software} N={size}: VCF not found ({vcf_path}) — skip")
        return None

    print(f"[ANALYZE] {software} N={size}: extracting GQ …")
    hist, total = extract_gq_histogram(vcf_path)

    if total == 0:
        print(f"[WARN] {software} N={size}: no GQ values extracted", file=sys.stderr)
        return None

    data = {
        "software":       software,
        "dataset_size":   size,
        "gq_histogram":   hist,
        "total_genotypes": total,
        "vcf_path":       vcf_path,
    }
    with open(out, "w") as f:
        json.dump(data, f, indent=2)
    print(f"[ANALYZE] Written: {out}")
    return data


def analyze_overlap(metrics_dir: Path, size: int, force: bool) -> None:
    """Compute variant overlap across all tools for one dataset size."""
    out = metrics_dir / f"variant_overlap_{size}.json"

    if out.exists() and not force:
        print(f"[OVERLAP] N={size}: already done — skip")
        return

    id_sets: dict[str, set] = {}
    for software in SOFTWARES:
        rec = _load_run_metrics(metrics_dir, software, size)
        if rec is None or rec.get("status") != "success":
            continue
        vcf_path = rec.get("output_vcf", "")
        if not vcf_path or not os.path.isfile(vcf_path):
            continue
        print(f"[OVERLAP] {software} N={size}: extracting variant IDs …")
        ids = extract_variant_ids(vcf_path)
        if ids:
            id_sets[software] = ids
            print(f"[OVERLAP] {software} N={size}: {len(ids):,} variants")

    if len(id_sets) < 2:
        print(f"[OVERLAP] N={size}: fewer than 2 tools with results — skip")
        return

    overlap = compute_overlap(id_sets)
    overlap["dataset_size"] = size

    with open(out, "w") as f:
        json.dump(overlap, f, indent=2)
    print(f"[OVERLAP] Written: {out}")


# ─── Main ─────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(
        description="Extract GQ histograms and variant overlap for benchmark report"
    )
    ap.add_argument("--benchmark-dir", required=True,
                    help="Root benchmarks directory (contains 03_metrics/)")
    ap.add_argument("--force", action="store_true",
                    help="Re-run even if output JSON files already exist")
    args = ap.parse_args()

    metrics_dir = Path(args.benchmark_dir) / "03_metrics"
    if not metrics_dir.exists():
        print(f"[ERROR] Metrics directory not found: {metrics_dir}", file=sys.stderr)
        sys.exit(1)

    # Verify bcftools is available
    try:
        ver = subprocess.run(
            ["bcftools", "--version"], capture_output=True, text=True, check=True
        )
        print(f"[ANALYZE] {ver.stdout.splitlines()[0]}")
    except (subprocess.CalledProcessError, FileNotFoundError):
        print("[ERROR] bcftools not found in PATH — required for variant analysis",
              file=sys.stderr)
        sys.exit(1)

    # GQ histograms (independent per software × size)
    for sw in SOFTWARES:
        for size in SIZES:
            analyze_gq(metrics_dir, sw, size, args.force)

    # Inter-tool variant overlap
    for size in SIZES:
        analyze_overlap(metrics_dir, size, args.force)

    print("[ANALYZE] Done.")


if __name__ == "__main__":
    main()
