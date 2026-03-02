#!/usr/bin/env bash
# =============================================================================
# GVCF Scalability Benchmark — Configuration
# =============================================================================

# ─── Reference Genome ─────────────────────────────────────────────────────────
# Required for GATK and Parabricks. GLnexus does not need a reference.
# File must have a .fai index and a .dict sequence dictionary.
REF_GENOME="/home/alisongonpereira/raid/largefiles/hg38_fasta/Homo_sapiens_assembly38.fasta"

# ─── Resource limits (80% of 128-core / 504G server) ─────────────────────────
# Adjust if running on a different machine.
# 80% of 128 cores = 102 | 80% of 504G RAM = ~400G
BENCHMARK_THREADS="${BENCHMARK_THREADS:-102}"
BENCHMARK_MEM_GB="${BENCHMARK_MEM_GB:-400}"


# jemalloc dramatically speeds up GLnexus's allocator-heavy workload.
# IMPORTANT: scripts/02_run_gatk.sh uses `env -u LD_PRELOAD` before every
# `gatk` invocation to prevent this from being inherited by the JVM —
# libjemalloc + JVM = SIGSEGV (exit 245).  Do NOT remove that env -u guard.
export LD_PRELOAD="${CONDA_PREFIX}/lib/libjemalloc.so:${LD_PRELOAD:-}"

# ─── GLnexus ──────────────────────────────────────────────────────────────────
# Preset config: DeepVariantWGS, DeepVariantWES, gatk, gatk_unphased, etc.
GLNEXUS_CONFIG="${GLNEXUS_CONFIG:-DeepVariantWES}"
# Path to glnexus_cli (leave as-is if it is in PATH)
GLNEXUS_BIN="${GLNEXUS_BIN:-glnexus_cli}"
# Worker threads for glnexus_cli (default: all available → pin to 80%)
GLNEXUS_THREADS="${GLNEXUS_THREADS:-${BENCHMARK_THREADS}}"
# Memory budget in GiB passed to glnexus_cli (0 = unlimited)
GLNEXUS_MEM_GB="${GLNEXUS_MEM_GB:-${BENCHMARK_MEM_GB}}"
# Threads for the bcftools view BCF→VCF pipe
BCFTOOLS_THREADS="${BCFTOOLS_THREADS:-${BENCHMARK_THREADS}}"

# ─── NVIDIA Parabricks ────────────────────────────────────────────────────────
# GPU device(s) to use, e.g. "0" for the first A100
PARABRICKS_GPU="${PARABRICKS_GPU:-0}"

# Parabricks via Docker (pbrun inside container)
# Note: mount /nfs/theseus and $HOME so container sees your data + reference.
PARABRICKS_DOCKER_IMAGE="${PARABRICKS_DOCKER_IMAGE:-nvcr.io/nvidia/clara/clara-parabricks:4.5.1-1}"
PARABRICKS_BIN="${PARABRICKS_BIN:-docker run --rm --gpus \"device=${PARABRICKS_GPU}\" -v /nfs/theseus:/nfs/theseus -v ${HOME}:${HOME} -w ${PWD} ${PARABRICKS_DOCKER_IMAGE} pbrun}"

# ─── GATK ─────────────────────────────────────────────────────────────────────
# Java heap = 80% of RAM; ParallelGC uses 80% of CPU cores for GC.
GATK_JAVA_OPTS="${GATK_JAVA_OPTS:--Xmx${BENCHMARK_MEM_GB}g -XX:+UseParallelGC -XX:ParallelGCThreads=${BENCHMARK_THREADS}}"
# gatk wrapper path
GATK_BIN="${GATK_BIN:-gatk}"

# ─── GATK GenomicsDBImport ────────────────────────────────────────────────────
# Optional: path to an intervals file (BED or .list) for GenomicsDBImport.
# If empty, main chromosomes (chr[0-9]+, chrX, chrY, chrM) are auto-detected
# from ${REF_GENOME}.fai.  For WES set this to your capture-kit BED.
GENOMICSDB_INTERVALS="${GENOMICSDB_INTERVALS:-}"
# Parallel reader threads for GenomicsDBImport — feed all 80% cores
GENOMICSDB_READER_THREADS="${GENOMICSDB_READER_THREADS:-${BENCHMARK_THREADS}}"
# Batch size for GenomicsDBImport (reduce below 50 if you hit OOM)
GENOMICSDB_BATCH_SIZE="${GENOMICSDB_BATCH_SIZE:-50}"

# ─── Benchmark Control ────────────────────────────────────────────────────────
# Sizes (number of GVCFs per dataset run)
DATASET_SIZES=(10 20 30 40 50 60 70 80 90 100)
# Reproducible random selection seed
RANDOM_SEED="${RANDOM_SEED:-42}"
# Resource-monitor sampling interval (seconds)
MONITOR_INTERVAL="${MONITOR_INTERVAL:-5}"
# Which softwares to benchmark — comment out any to skip
BENCHMARK_SOFTWARES=("glnexus" "parabricks" "gatk" "gatk_genomicsdb")

# ─── Paths (derived; normally no need to change) ──────────────────────────────
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
INPUT_DIR="${SCRIPT_DIR}/input_gvcfs"
BENCHMARK_DIR="${SCRIPT_DIR}/benchmarks"
SCRIPTS_DIR="${SCRIPT_DIR}/scripts"
LOG_FILE="${BENCHMARK_DIR}/04_reports/execution.log"