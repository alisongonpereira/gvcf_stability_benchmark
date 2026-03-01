#!/usr/bin/env bash
# =============================================================================
# GVCF Scalability Benchmark — Configuration
# Edit this file to match your environment before running run_benchmark.sh
# =============================================================================

# ─── Reference Genome ─────────────────────────────────────────────────────────
# Required for GATK and Parabricks. GLnexus does not need a reference.
# File must have a .fai index and a .dict sequence dictionary.
REF_GENOME="${REF_GENOME:-}"

# ─── GLnexus ──────────────────────────────────────────────────────────────────
# Preset config: DeepVariantWGS, DeepVariantWES, gatk, gatk_unphased, etc.
GLNEXUS_CONFIG="${GLNEXUS_CONFIG:-DeepVariantWGS}"
# Path to glnexus_cli (leave as-is if it is in PATH)
GLNEXUS_BIN="${GLNEXUS_BIN:-glnexus_cli}"

# ─── NVIDIA Parabricks ────────────────────────────────────────────────────────
# GPU device(s) to use, e.g. "0" for the first A100
PARABRICKS_GPU="${PARABRICKS_GPU:-0}"
# pbrun binary path
PARABRICKS_BIN="${PARABRICKS_BIN:-pbrun}"

# ─── GATK ─────────────────────────────────────────────────────────────────────
# Java options — tune heap based on available RAM
GATK_JAVA_OPTS="${GATK_JAVA_OPTS:--Xmx32g -XX:+UseParallelGC}"
# gatk wrapper path
GATK_BIN="${GATK_BIN:-gatk}"

# ─── Benchmark Control ────────────────────────────────────────────────────────
# Sizes (number of GVCFs per dataset run)
DATASET_SIZES=(10 20 30 40 50 60 70 80 90 100)
# Reproducible random selection seed
RANDOM_SEED="${RANDOM_SEED:-42}"
# Resource-monitor sampling interval (seconds)
MONITOR_INTERVAL="${MONITOR_INTERVAL:-5}"
# Which softwares to benchmark — comment out any to skip
BENCHMARK_SOFTWARES=("glnexus" "parabricks" "gatk")

# ─── Paths (derived; normally no need to change) ──────────────────────────────
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
INPUT_DIR="${SCRIPT_DIR}/input_gvcfs"
BENCHMARK_DIR="${SCRIPT_DIR}/benchmarks"
SCRIPTS_DIR="${SCRIPT_DIR}/scripts"
LOG_FILE="${BENCHMARK_DIR}/04_reports/execution.log"
