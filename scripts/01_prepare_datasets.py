#!/usr/bin/env python3
"""
Step 1 — Dataset Preparation

Randomly selects GVCFs from the input pool and creates per-size dataset
directories with symlinks and a manifest (manifest.txt).  Reproducible via
--seed.  Skips sizes that are already prepared unless --force is given.

Usage:
    python3 01_prepare_datasets.py \
        --input-dir   /path/to/input_gvcfs \
        --output-dir  /path/to/benchmarks/01_prep \
        --sizes 10 20 30 40 50 60 70 80 90 100 \
        --seed 42 \
        --log  /path/to/preparation.log
"""

import argparse
import json
import logging
import os
import random
import sys
from datetime import datetime
from pathlib import Path


# ─── Logging setup ────────────────────────────────────────────────────────────

def _setup_logging(log_file: str) -> logging.Logger:
    logger = logging.getLogger("prep")
    logger.setLevel(logging.DEBUG)
    fmt = logging.Formatter("[%(asctime)s] [PREP] %(message)s",
                            datefmt="%Y-%m-%d %H:%M:%S")

    sh = logging.StreamHandler(sys.stdout)
    sh.setFormatter(fmt)
    logger.addHandler(sh)

    Path(log_file).parent.mkdir(parents=True, exist_ok=True)
    fh = logging.FileHandler(log_file)
    fh.setFormatter(fmt)
    logger.addHandler(fh)

    return logger


# ─── GVCF discovery ───────────────────────────────────────────────────────────

GVCF_EXTENSIONS = {".g.vcf", ".gvcf", ".g.vcf.gz", ".gvcf.gz"}


def discover_gvcfs(input_dir: Path) -> list[Path]:
    """Return sorted list of GVCF paths in input_dir (non-recursive)."""
    found = []
    for f in sorted(input_dir.iterdir()):
        if f.is_file():
            suffixes = "".join(f.suffixes)
            if suffixes in GVCF_EXTENSIONS or any(
                str(f).endswith(ext) for ext in GVCF_EXTENSIONS
            ):
                found.append(f)
    return found


# ─── Per-dataset creation ─────────────────────────────────────────────────────

def prepare_dataset(
    size: int,
    pool: list[Path],
    output_dir: Path,
    rng: random.Random,
    logger: logging.Logger,
    force: bool = False,
) -> dict:
    """
    Create dataset_<size>/ with symlinks and manifest.txt.
    Returns a dict describing the dataset.
    """
    dataset_dir = output_dir / f"dataset_{size}"
    done_flag   = dataset_dir / ".done"

    if done_flag.exists() and not force:
        logger.info(f"dataset_{size}: already prepared — skipping (use --force to redo)")
        # Re-read manifest
        manifest = dataset_dir / "manifest.txt"
        selected = manifest.read_text().splitlines() if manifest.exists() else []
        return {"size": size, "status": "skipped", "files": selected}

    if size > len(pool):
        logger.warning(
            f"dataset_{size}: requested {size} GVCFs but pool only has {len(pool)} — skipping"
        )
        return {"size": size, "status": "insufficient_pool", "files": []}

    dataset_dir.mkdir(parents=True, exist_ok=True)

    selected = rng.sample(pool, size)
    selected_sorted = sorted(selected, key=lambda p: p.name)

    symlink_dir = dataset_dir / "gvcfs"
    symlink_dir.mkdir(exist_ok=True)

    symlink_paths = []
    for gvcf in selected_sorted:
        link = symlink_dir / gvcf.name
        if link.exists() or link.is_symlink():
            link.unlink()
        link.symlink_to(gvcf.resolve())
        # Also symlink companion index files (.tbi, .csi) if present
        for idx_ext in (".tbi", ".csi"):
            idx = Path(str(gvcf) + idx_ext)
            if idx.exists():
                idx_link = symlink_dir / (gvcf.name + idx_ext)
                if idx_link.exists() or idx_link.is_symlink():
                    idx_link.unlink()
                idx_link.symlink_to(idx.resolve())
        symlink_paths.append(str(link.resolve()))

    # Write manifest (one absolute path per line)
    manifest_path = dataset_dir / "manifest.txt"
    manifest_path.write_text("\n".join(symlink_paths) + "\n")

    # Write JSON metadata
    meta = {
        "dataset_size": size,
        "created_at": datetime.now().isoformat(),
        "manifest": symlink_paths,
        "source_files": [str(p) for p in selected_sorted],
    }
    (dataset_dir / "dataset_info.json").write_text(
        json.dumps(meta, indent=2) + "\n"
    )

    # Mark done
    done_flag.write_text(datetime.now().isoformat() + "\n")

    logger.info(f"dataset_{size}: created {size} symlinks in {symlink_dir}")
    for p in selected_sorted:
        logger.info(f"  → {p.name}")

    return {"size": size, "status": "created", "files": symlink_paths}


# ─── Main ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="Prepare GVCF benchmark datasets")
    parser.add_argument("--input-dir",  required=True, help="Directory with GVCF pool")
    parser.add_argument("--output-dir", required=True, help="benchmarks/01_prep directory")
    parser.add_argument("--sizes", type=int, nargs="+", required=True,
                        help="Dataset sizes, e.g. 10 20 30 ... 100")
    parser.add_argument("--seed",  type=int, default=42, help="Random seed")
    parser.add_argument("--log",   default="preparation.log", help="Log file path")
    parser.add_argument("--force", action="store_true",
                        help="Recreate datasets even if already done")
    args = parser.parse_args()

    logger = _setup_logging(args.log)
    logger.info("=" * 60)
    logger.info("GVCF Scalability Benchmark — Dataset Preparation")
    logger.info(f"Input directory : {args.input_dir}")
    logger.info(f"Output directory: {args.output_dir}")
    logger.info(f"Dataset sizes   : {args.sizes}")
    logger.info(f"Random seed     : {args.seed}")
    logger.info("=" * 60)

    input_dir  = Path(args.input_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    if not input_dir.exists():
        logger.error(f"Input directory does not exist: {input_dir}")
        sys.exit(1)

    pool = discover_gvcfs(input_dir)
    logger.info(f"Discovered {len(pool)} GVCF files in pool")

    if len(pool) == 0:
        logger.error("No GVCF files found. Expected extensions: " +
                     ", ".join(GVCF_EXTENSIONS))
        sys.exit(1)

    # Write the full pool list
    pool_list = output_dir / "all_gvcfs.txt"
    pool_list.write_text("\n".join(str(p) for p in pool) + "\n")
    logger.info(f"Pool list written to {pool_list}")

    rng = random.Random(args.seed)

    results = []
    for size in sorted(args.sizes):
        result = prepare_dataset(size, pool, output_dir, rng, logger, args.force)
        results.append(result)

    # Summary
    logger.info("-" * 60)
    for r in results:
        logger.info(f"  dataset_{r['size']:>3}: {r['status']} ({len(r['files'])} files)")

    created = sum(1 for r in results if r["status"] == "created")
    skipped = sum(1 for r in results if r["status"] == "skipped")
    failed  = sum(1 for r in results if r["status"] not in ("created", "skipped"))
    logger.info(f"Summary: {created} created, {skipped} skipped, {failed} failed")
    logger.info("Dataset preparation complete.")

    if failed > 0:
        sys.exit(1)


if __name__ == "__main__":
    main()
