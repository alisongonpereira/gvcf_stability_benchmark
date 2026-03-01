#!/usr/bin/env python3
"""
Real-time resource monitor for GVCF benchmark runs.

Usage:
    python3 monitor_resources.py output.json [--interval 5] [--pid 1234]

Runs until SIGTERM/SIGINT; writes JSON with per-sample readings and aggregates.
"""
import argparse
import json
import os
import signal
import subprocess
import sys
import time
from datetime import datetime

# Optional: psutil provides richer metrics
try:
    import psutil
    HAS_PSUTIL = True
except ImportError:
    HAS_PSUTIL = False


# ─── GPU helpers ──────────────────────────────────────────────────────────────

def _gpu_stats():
    """Return list of per-GPU dicts via nvidia-smi; empty list on failure."""
    try:
        result = subprocess.run(
            [
                "nvidia-smi",
                "--query-gpu=name,utilization.gpu,utilization.memory,"
                "memory.used,memory.total,temperature.gpu",
                "--format=csv,noheader,nounits",
            ],
            capture_output=True,
            text=True,
            timeout=15,
        )
        if result.returncode != 0:
            return []
        gpus = []
        for line in result.stdout.strip().splitlines():
            parts = [p.strip() for p in line.split(",")]
            if len(parts) >= 6:
                try:
                    gpus.append(
                        {
                            "name": parts[0],
                            "util_pct": float(parts[1]),
                            "mem_util_pct": float(parts[2]),
                            "mem_used_mb": float(parts[3]),
                            "mem_total_mb": float(parts[4]),
                            "temperature_c": float(parts[5]),
                        }
                    )
                except ValueError:
                    pass
        return gpus
    except Exception:
        return []


# ─── CPU / RAM helpers (fallback when psutil missing) ─────────────────────────

def _cpu_pct_proc():
    """Rough system CPU % from /proc/stat (non-idle fraction)."""
    try:
        with open("/proc/stat") as f:
            line = f.readline()
        vals = list(map(int, line.split()[1:]))
        total = sum(vals)
        idle = vals[3]
        if not hasattr(_cpu_pct_proc, "_prev"):
            _cpu_pct_proc._prev = (total, idle)
            return 0.0
        prev_total, prev_idle = _cpu_pct_proc._prev
        d_total = total - prev_total
        d_idle = idle - prev_idle
        _cpu_pct_proc._prev = (total, idle)
        return round(100.0 * (1 - d_idle / d_total) if d_total else 0.0, 1)
    except Exception:
        return None


def _ram_gb_proc():
    """Return (used_gb, available_gb, total_gb) from /proc/meminfo."""
    try:
        info = {}
        with open("/proc/meminfo") as f:
            for line in f:
                key, val = line.split(":", 1)
                info[key.strip()] = int(val.split()[0])  # kB
        total = info.get("MemTotal", 0)
        avail = info.get("MemAvailable", info.get("MemFree", 0))
        used = total - avail
        gb = 1024 * 1024
        return round(used / gb, 3), round(avail / gb, 3), round(total / gb, 3)
    except Exception:
        return None, None, None


def _disk_gb(path="/"):
    """Return (used_gb, free_gb, total_gb) for filesystem containing path."""
    try:
        st = os.statvfs(path)
        total = st.f_frsize * st.f_blocks / (1024**3)
        free  = st.f_frsize * st.f_bfree  / (1024**3)
        used  = total - free
        return round(used, 3), round(free, 3), round(total, 3)
    except Exception:
        return None, None, None


# ─── Sample collection ────────────────────────────────────────────────────────

def collect_sample(start_time: float) -> dict:
    sample: dict = {
        "timestamp": datetime.now().isoformat(),
        "elapsed_s": round(time.time() - start_time, 2),
    }

    if HAS_PSUTIL:
        sample["cpu_pct"] = psutil.cpu_percent(interval=None)
        mem = psutil.virtual_memory()
        gb = 1024**3
        sample["ram_used_gb"]      = round(mem.used      / gb, 3)
        sample["ram_available_gb"] = round(mem.available / gb, 3)
        sample["ram_total_gb"]     = round(mem.total     / gb, 3)
        try:
            du = psutil.disk_usage("/")
            sample["disk_used_gb"]  = round(du.used  / gb, 3)
            sample["disk_free_gb"]  = round(du.free  / gb, 3)
            sample["disk_total_gb"] = round(du.total / gb, 3)
        except Exception:
            pass
    else:
        sample["cpu_pct"] = _cpu_pct_proc()
        used, avail, total = _ram_gb_proc()
        if used is not None:
            sample["ram_used_gb"]      = used
            sample["ram_available_gb"] = avail
            sample["ram_total_gb"]     = total
        du, df, dt = _disk_gb("/")
        if du is not None:
            sample["disk_used_gb"]  = du
            sample["disk_free_gb"]  = df
            sample["disk_total_gb"] = dt

    sample["gpus"] = _gpu_stats()
    return sample


# ─── Aggregation ─────────────────────────────────────────────────────────────

def _agg(vals):
    if not vals:
        return {}
    return {
        "avg": round(sum(vals) / len(vals), 3),
        "max": round(max(vals), 3),
        "min": round(min(vals), 3),
    }


def compute_aggregates(samples: list) -> dict:
    agg = {}
    def _extract(key):
        return [s[key] for s in samples if isinstance(s.get(key), (int, float))]

    for metric in ("cpu_pct", "ram_used_gb", "ram_available_gb",
                   "disk_used_gb", "disk_free_gb"):
        vals = _extract(metric)
        if vals:
            agg[metric] = _agg(vals)

    # GPU (aggregate across first GPU for simplicity; all GPUs stored in samples)
    for gkey in ("util_pct", "mem_util_pct", "mem_used_mb"):
        vals = [s["gpus"][0][gkey] for s in samples
                if s.get("gpus") and isinstance(s["gpus"][0].get(gkey), (int, float))]
        if vals:
            agg[f"gpu_{gkey}"] = _agg(vals)

    # Convenience flat aliases expected by the report generator
    if "cpu_pct" in agg:
        agg["cpu_pct_avg"] = agg["cpu_pct"]["avg"]
        agg["cpu_pct_max"] = agg["cpu_pct"]["max"]
    if "ram_used_gb" in agg:
        agg["ram_used_gb_avg"] = agg["ram_used_gb"]["avg"]
        agg["ram_used_gb_max"] = agg["ram_used_gb"]["max"]
    if "ram_available_gb" in agg:
        agg["ram_available_gb_min"] = agg["ram_available_gb"]["min"]
    if "gpu_util_pct" in agg:
        agg["gpu_util_pct_avg"] = agg["gpu_util_pct"]["avg"]
        agg["gpu_util_pct_max"] = agg["gpu_util_pct"]["max"]
    if "gpu_mem_used_mb" in agg:
        agg["gpu_mem_used_mb_avg"] = agg["gpu_mem_used_mb"]["avg"]
        agg["gpu_mem_used_mb_max"] = agg["gpu_mem_used_mb"]["max"]
    if "disk_free_gb" in agg:
        agg["disk_free_gb_final"] = agg["disk_free_gb"]["min"]
    if "disk_used_gb" in agg:
        agg["disk_used_gb_max"] = agg["disk_used_gb"]["max"]

    return agg


# ─── Main monitor loop ────────────────────────────────────────────────────────

def monitor(output_file: str, interval: float, target_pid: int | None):
    samples = []
    running = True
    start_time = time.time()

    def _stop(signum, frame):
        nonlocal running
        running = False

    signal.signal(signal.SIGTERM, _stop)
    signal.signal(signal.SIGINT, _stop)

    # Prime psutil cpu_percent (first call always returns 0)
    if HAS_PSUTIL:
        psutil.cpu_percent(interval=None)

    tmp_file = output_file + ".tmp"

    while running:
        sample = collect_sample(start_time)

        # If monitoring a specific PID, check it is still alive
        if target_pid:
            try:
                if HAS_PSUTIL:
                    proc = psutil.Process(target_pid)
                    children = proc.children(recursive=True)
                    all_procs = [proc] + children
                    proc_cpu = sum(p.cpu_percent() for p in all_procs if p.is_running())
                    proc_mem = sum(p.memory_info().rss for p in all_procs if p.is_running())
                    sample["proc_cpu_pct"] = round(proc_cpu, 1)
                    sample["proc_ram_gb"]  = round(proc_mem / 1024**3, 3)
                else:
                    # Check existence via /proc
                    if not os.path.exists(f"/proc/{target_pid}"):
                        running = False
            except Exception:
                pass

        samples.append(sample)

        # Periodic intermediate save
        if len(samples) % 6 == 0:
            with open(tmp_file, "w") as f:
                json.dump({"samples": samples}, f)

        time.sleep(interval)

    aggregates = compute_aggregates(samples)
    result = {
        "start_time":   datetime.fromtimestamp(start_time).isoformat(),
        "end_time":     datetime.now().isoformat(),
        "duration_s":   round(time.time() - start_time, 2),
        "sample_count": len(samples),
        "interval_s":   interval,
        "has_psutil":   HAS_PSUTIL,
        "aggregates":   aggregates,
        "samples":      samples,
    }
    with open(output_file, "w") as f:
        json.dump(result, f, indent=2)

    if os.path.exists(tmp_file):
        os.remove(tmp_file)


# ─── Entry point ─────────────────────────────────────────────────────────────

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Benchmark resource monitor")
    parser.add_argument("output_file", help="Path for the output JSON file")
    parser.add_argument("--interval", type=float, default=5.0,
                        help="Sampling interval in seconds (default: 5)")
    parser.add_argument("--pid", type=int, default=None,
                        help="PID of process to monitor (optional)")
    args = parser.parse_args()

    monitor(args.output_file, args.interval, args.pid)
