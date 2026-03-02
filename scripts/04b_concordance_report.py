#!/usr/bin/env python3
"""
Step 4b — Inter-tool Concordance Report

Reads variant overlap JSON files produced by 03_analyze_variants.py and generates:
  • benchmarks/04_reports/concordance_report.html

Focus: for each dataset size × replicate, how much do the tools agree on which
       variants to call when given the **same set of input samples**?

Key metrics:
  • Jaccard similarity  = |A ∩ B| / |A ∪ B|  (per tool pair per size per replicate)
  • Exclusive variants  = variants called by exactly one tool
  • Concordance stability = mean ± SD Jaccard across replicates
  • Alerts when Jaccard < WARNING_THRESHOLD (0.80) or < CRITICAL_THRESHOLD (0.70)

Usage:
    python3 04b_concordance_report.py \
        --benchmark-dir /path/to/benchmarks \
        --output-dir    /path/to/benchmarks/04_reports
"""

import argparse
import json
import math
import re
import sys
from datetime import datetime
from itertools import combinations
from pathlib import Path

# ─── Configuration ────────────────────────────────────────────────────────────

SOFTWARES = ["glnexus", "parabricks", "gatk", "gatk_genomicsdb"]
SOFT_LABEL = {
    "glnexus":         "GLnexus",
    "parabricks":      "Parabricks",
    "gatk":            "GATK (CombineGVCFs)",
    "gatk_genomicsdb": "GATK (GenomicsDB)",
}
SOFT_COLOR = {
    "glnexus":         "#1f77b4",
    "parabricks":      "#ff7f0e",
    "gatk":            "#2ca02c",
    "gatk_genomicsdb": "#9467bd",
}

WARNING_THRESHOLD  = 0.80
CRITICAL_THRESHOLD = 0.70

PLOTLY_CDN = "https://cdn.plot.ly/plotly-2.35.0.min.js"

SIZES      = [10, 25, 50, 75, 100]
REPLICATES = [1, 2, 3]


# ─── Statistics helpers ───────────────────────────────────────────────────────

def _mean(vals):
    return sum(vals) / len(vals) if vals else 0.0


def _std(vals):
    if len(vals) < 2:
        return 0.0
    m = _mean(vals)
    return math.sqrt(sum((v - m) ** 2 for v in vals) / (len(vals) - 1))


# ─── Data loading ─────────────────────────────────────────────────────────────

def load_overlap_data(metrics_dir: Path) -> dict:
    """
    Returns {size: {rep: overlap_dict}}.
    overlap_dict keys (from 03_analyze_variants.py):
        tools, counts, exclusive, pairwise_only, all_common,
        total_unique, count_{t1}_{t2}, dataset_size, replicate
    """
    result: dict[int, dict[int, dict]] = {}

    for path in sorted(metrics_dir.glob("variant_overlap_*.json")):
        # New format: variant_overlap_{size}_rep{r}.json
        m = re.match(r"variant_overlap_(\d+)_rep(\d+)\.json", path.name)
        if m:
            size, rep = int(m.group(1)), int(m.group(2))
        else:
            # Legacy format: variant_overlap_{size}.json → rep 1
            m = re.match(r"variant_overlap_(\d+)\.json", path.name)
            if not m:
                continue
            size, rep = int(m.group(1)), 1

        try:
            with open(path) as f:
                rec = json.load(f)
        except Exception as e:
            print(f"[WARN] Could not read {path}: {e}", file=sys.stderr)
            continue

        result.setdefault(size, {})[rep] = rec

    return result


# ─── Jaccard computation ──────────────────────────────────────────────────────

def jaccard(overlap: dict, t1: str, t2: str) -> float | None:
    """
    Compute Jaccard(t1, t2) = |t1 ∩ t2| / |t1 ∪ t2|
    from the overlap dict produced by compute_overlap().
    """
    counts = overlap.get("counts", {})
    n1 = counts.get(t1)
    n2 = counts.get(t2)
    if n1 is None or n2 is None:
        return None

    # Intersection count is stored as count_{t1}_{t2} or count_{t2}_{t1}
    inter = overlap.get(f"count_{t1}_{t2}") or overlap.get(f"count_{t2}_{t1}")
    if inter is None:
        return None

    union = n1 + n2 - inter
    if union == 0:
        return 1.0
    return round(inter / union, 6)


def all_pairs() -> list[tuple[str, str]]:
    """All ordered pairs of softwares present in SOFTWARES."""
    return list(combinations(SOFTWARES, 2))


def compute_jaccard_table(overlap_data: dict) -> dict:
    """
    Returns {(t1,t2): {size: {rep: jaccard_value}}}.
    """
    table: dict[tuple, dict] = {}
    for pair in all_pairs():
        t1, t2 = pair
        table[pair] = {}
        for size, rep_dict in overlap_data.items():
            table[pair][size] = {}
            for rep, ov in rep_dict.items():
                tools_present = ov.get("tools", [])
                if t1 not in tools_present or t2 not in tools_present:
                    continue
                j = jaccard(ov, t1, t2)
                if j is not None:
                    table[pair][size][rep] = j
    return table


# ─── Alert detection ──────────────────────────────────────────────────────────

def find_alerts(jaccard_table: dict) -> list[dict]:
    """Return list of alert dicts for low-concordance pairs."""
    alerts = []
    for (t1, t2), size_dict in jaccard_table.items():
        for size, rep_dict in size_dict.items():
            for rep, j in rep_dict.items():
                if j < CRITICAL_THRESHOLD:
                    level = "CRÍTICO"
                    color = "#d62728"
                elif j < WARNING_THRESHOLD:
                    level = "AVISO"
                    color = "#ff7f0e"
                else:
                    continue
                alerts.append({
                    "level":  level,
                    "color":  color,
                    "t1":     t1,
                    "t2":     t2,
                    "size":   size,
                    "rep":    rep,
                    "jaccard": j,
                })
    alerts.sort(key=lambda a: a["jaccard"])
    return alerts


# ─── HTML/JS helpers ──────────────────────────────────────────────────────────

_CSS = """
* { box-sizing: border-box; margin: 0; padding: 0; }
body { font-family: 'Segoe UI', Arial, sans-serif; background: #f5f6fa; color: #333; }
header { background: linear-gradient(135deg,#1A1A2E,#16213E); color: #fff;
         padding: 24px 32px; }
header h1 { font-size: 1.6rem; margin-bottom: 6px; }
header p  { opacity: .7; font-size: .9rem; }
.tabs { display: flex; gap: 4px; padding: 12px 32px 0;
        background: #fff; border-bottom: 2px solid #e0e0e0; flex-wrap: wrap; }
.tab-btn { padding: 8px 18px; border: none; border-radius: 6px 6px 0 0;
           cursor: pointer; background: #f0f0f0; font-size: .9rem; }
.tab-btn.active { background: #1A1A2E; color: #fff; }
.tab-content { display: none; padding: 24px 32px; }
.tab-content.active { display: block; }
.card { background: #fff; border-radius: 8px; padding: 20px;
        box-shadow: 0 1px 4px rgba(0,0,0,.08); margin-bottom: 16px; }
.section-title { font-size: 1.1rem; font-weight: 600; color: #1A1A2E;
                 margin: 24px 0 12px; border-left: 4px solid #4a90e2;
                 padding-left: 10px; }
table { border-collapse: collapse; width: 100%; font-size: .875rem; }
th, td { padding: 8px 12px; border: 1px solid #e0e0e0; text-align: center; }
th { background: #1A1A2E; color: #fff; font-weight: 600; }
tr:nth-child(even) { background: #f9f9f9; }
.alert-box { border-radius: 6px; padding: 12px 18px; margin-bottom: 8px;
             border-left: 5px solid; }
.grid-2 { display: grid; grid-template-columns: 1fr 1fr; gap: 16px; }
.grid-3 { display: grid; grid-template-columns: 1fr 1fr 1fr; gap: 16px; }
@media (max-width: 900px) { .grid-2, .grid-3 { grid-template-columns: 1fr; } }
"""

_JS_TAB = """
function showTab(id, btn) {
  document.querySelectorAll('.tab-content').forEach(t => t.classList.remove('active'));
  document.querySelectorAll('.tab-btn').forEach(b => b.classList.remove('active'));
  document.getElementById(id).classList.add('active');
  btn.classList.add('active');
}
"""

_CHART_LAYOUT_BASE = {
    "margin": {"l": 60, "r": 20, "t": 40, "b": 60},
    "paper_bgcolor": "rgba(0,0,0,0)",
    "plot_bgcolor": "#fafafa",
    "font": {"family": "Segoe UI, Arial", "size": 12},
    "xaxis": {"title": "Tamanho do dataset (# GVCFs)", "gridcolor": "#eee"},
    "yaxis": {"gridcolor": "#eee"},
    "legend": {"orientation": "h", "y": -0.25},
    "hovermode": "x unified",
}

PAIR_COLORS = [
    "#e41a1c", "#377eb8", "#4daf4a", "#984ea3",
    "#ff7f00", "#a65628",
]


def _pair_label(t1: str, t2: str) -> str:
    return f"{SOFT_LABEL.get(t1,t1)} vs {SOFT_LABEL.get(t2,t2)}"


# ─── Tab builders ─────────────────────────────────────────────────────────────

def build_overview_tab(overlap_data: dict, jaccard_table: dict) -> tuple[str, str]:
    """Summary heatmap-style table + mean Jaccard line chart."""
    html_parts = [
        '<div class="section-title">Concordância inter-ferramenta — Visão Geral</div>',
        '<p style="margin-bottom:16px;color:#555;font-size:.9rem">'
        'Jaccard = |A∩B|/|A∪B|. Valor 1.0 = concordância perfeita; '
        f'valores abaixo de {WARNING_THRESHOLD} são marcados como Aviso e abaixo de '
        f'{CRITICAL_THRESHOLD} como Crítico.</p>',
    ]
    js_parts = []

    # Mean Jaccard table per size (averaged across reps)
    sizes_with_data = sorted(overlap_data.keys())

    html_parts.append('<div class="card"><table>')
    header = ["Par de Ferramentas"] + [f"N={s}" for s in sizes_with_data] + ["Média", "DP"]
    html_parts.append("<tr>" + "".join(f"<th>{h}</th>" for h in header) + "</tr>")

    for ci, (t1, t2) in enumerate(all_pairs()):
        pair_key = (t1, t2)
        row_cells = [f'<td style="text-align:left;font-weight:600">{_pair_label(t1,t2)}</td>']
        all_j = []
        for size in sizes_with_data:
            rep_dict = jaccard_table.get(pair_key, {}).get(size, {})
            vals = list(rep_dict.values())
            if not vals:
                row_cells.append("<td>—</td>")
                continue
            mean_j = _mean(vals)
            all_j.extend(vals)
            color = (
                "#2ca02c" if mean_j >= WARNING_THRESHOLD else
                "#ff7f0e" if mean_j >= CRITICAL_THRESHOLD else
                "#d62728"
            )
            row_cells.append(f'<td style="color:{color};font-weight:600">{mean_j:.3f}</td>')

        if all_j:
            mg = _mean(all_j)
            sg = _std(all_j)
            color_g = (
                "#2ca02c" if mg >= WARNING_THRESHOLD else
                "#ff7f0e" if mg >= CRITICAL_THRESHOLD else
                "#d62728"
            )
            row_cells.append(f'<td style="color:{color_g};font-weight:700">{mg:.3f}</td>')
            row_cells.append(f"<td>{sg:.3f}</td>")
        else:
            row_cells += ["<td>—</td>", "<td>—</td>"]

        html_parts.append("<tr>" + "".join(row_cells) + "</tr>")

    html_parts.append("</table></div>")

    # Jaccard line chart: mean across reps, one line per pair
    div_id = "chart_jaccard_overview"
    html_parts.append(
        f'<div class="card"><div class="section-title" style="font-size:1rem">'
        f'Jaccard médio por tamanho de dataset</div>'
        f'<div id="{div_id}" style="height:350px"></div></div>'
    )

    traces = []
    for ci, (t1, t2) in enumerate(all_pairs()):
        pair_key = (t1, t2)
        xs, ys, es = [], [], []
        for size in sizes_with_data:
            rep_dict = jaccard_table.get(pair_key, {}).get(size, {})
            vals = list(rep_dict.values())
            if not vals:
                continue
            xs.append(size)
            ys.append(round(_mean(vals), 4))
            es.append(round(_std(vals), 4))

        if not xs:
            continue

        color = PAIR_COLORS[ci % len(PAIR_COLORS)]
        d = {
            "x": xs, "y": ys,
            "name": _pair_label(t1, t2),
            "mode": "lines+markers",
            "line": {"color": color, "width": 2},
            "marker": {"size": 7, "color": color},
        }
        if any(e > 0 for e in es):
            d["error_y"] = {
                "type": "data", "array": es, "visible": True,
                "color": color, "thickness": 1.5, "width": 4,
            }
        traces.append(json.dumps(d))

    # Add threshold lines as shapes
    layout = dict(_CHART_LAYOUT_BASE)
    layout["title"] = {"text": "Jaccard médio por tamanho de dataset", "font": {"size": 14}}
    layout["yaxis"] = {"title": "Jaccard similarity", "range": [0, 1.05], "gridcolor": "#eee"}
    layout["shapes"] = [
        {"type": "line", "x0": sizes_with_data[0] if sizes_with_data else 0,
         "x1": sizes_with_data[-1] if sizes_with_data else 100,
         "y0": WARNING_THRESHOLD, "y1": WARNING_THRESHOLD,
         "line": {"color": "#ff7f0e", "dash": "dash", "width": 1.5}},
        {"type": "line", "x0": sizes_with_data[0] if sizes_with_data else 0,
         "x1": sizes_with_data[-1] if sizes_with_data else 100,
         "y0": CRITICAL_THRESHOLD, "y1": CRITICAL_THRESHOLD,
         "line": {"color": "#d62728", "dash": "dot", "width": 1.5}},
    ]
    layout["annotations"] = [
        {"x": sizes_with_data[-1] if sizes_with_data else 100,
         "y": WARNING_THRESHOLD + 0.02, "xanchor": "right",
         "text": f"Aviso ({WARNING_THRESHOLD})", "showarrow": False,
         "font": {"color": "#ff7f0e", "size": 10}},
        {"x": sizes_with_data[-1] if sizes_with_data else 100,
         "y": CRITICAL_THRESHOLD - 0.03, "xanchor": "right",
         "text": f"Crítico ({CRITICAL_THRESHOLD})", "showarrow": False,
         "font": {"color": "#d62728", "size": 10}},
    ]

    js_parts.append(
        f"Plotly.newPlot('{div_id}', [{', '.join(traces)}], "
        f"{json.dumps(layout)}, {{responsive: true}});"
    )

    return "\n".join(html_parts), "\n".join(js_parts)


def build_alerts_tab(alerts: list[dict]) -> str:
    """Alert boxes for all sub-threshold concordance pairs."""
    html_parts = [
        '<div class="section-title">Alertas de Concordância</div>',
    ]

    if not alerts:
        html_parts.append(
            '<div style="background:#d4edda;border-left:4px solid #28a745;'
            'padding:14px 20px;border-radius:4px">'
            '<strong>✓ Todas as concordâncias estão acima do limiar de aviso '
            f'({WARNING_THRESHOLD}).</strong>'
            '</div>'
        )
        return "\n".join(html_parts)

    html_parts.append(
        f'<p style="color:#555;margin-bottom:16px">'
        f'{len(alerts)} par(es) de ferramentas abaixo dos limiares. '
        f'Aviso: Jaccard &lt; {WARNING_THRESHOLD}. '
        f'Crítico: Jaccard &lt; {CRITICAL_THRESHOLD}.</p>'
    )

    for a in alerts:
        icon = "🚨" if a["level"] == "CRÍTICO" else "⚠️"
        html_parts.append(
            f'<div class="alert-box" style="background:#fff8f8;'
            f'border-color:{a["color"]};margin-bottom:10px">'
            f'<strong style="color:{a["color"]}">{icon} {a["level"]}</strong> — '
            f'{_pair_label(a["t1"], a["t2"])} &nbsp;|&nbsp; '
            f'N={a["size"]} rep{a["rep"]} &nbsp;|&nbsp; '
            f'Jaccard = <strong>{a["jaccard"]:.4f}</strong>'
            f'</div>'
        )

    return "\n".join(html_parts)


def build_heatmap_tab(overlap_data: dict, jaccard_table: dict) -> tuple[str, str]:
    """Jaccard heatmaps: one per dataset size (averaged over reps)."""
    html_parts = [
        '<div class="section-title">Heatmap de Concordância por Tamanho</div>',
        '<p style="color:#555;margin-bottom:16px;font-size:.9rem">'
        'Média do Jaccard sobre todas as réplicas disponíveis. '
        'Diagonal = 1.0 (ferramenta comparada com si mesma).</p>',
    ]
    js_parts = []

    sizes_with_data = sorted(overlap_data.keys())
    tools = SOFTWARES

    for size in sizes_with_data:
        div_id = f"chart_heatmap_{size}"
        html_parts.append(
            f'<div class="card" style="margin-bottom:16px">'
            f'<div style="font-weight:600;margin-bottom:8px">N = {size} amostras</div>'
            f'<div id="{div_id}" style="height:320px"></div></div>'
        )

        z_vals = []
        for t1 in tools:
            row = []
            for t2 in tools:
                if t1 == t2:
                    row.append(1.0)
                else:
                    pair_key = (t1, t2) if (t1, t2) in jaccard_table else (t2, t1)
                    rep_dict = jaccard_table.get(pair_key, {}).get(size, {})
                    vals = list(rep_dict.values())
                    row.append(round(_mean(vals), 4) if vals else None)
            z_vals.append(row)

        labels = [SOFT_LABEL.get(t, t) for t in tools]
        trace = json.dumps({
            "z": z_vals,
            "x": labels,
            "y": labels,
            "type": "heatmap",
            "colorscale": "RdYlGn",
            "zmin": 0,
            "zmax": 1,
            "text": [[f"{v:.3f}" if v is not None else "N/A" for v in row] for row in z_vals],
            "texttemplate": "%{text}",
            "hovertemplate": "%{y} vs %{x}: Jaccard=%{z:.4f}<extra></extra>",
            "colorbar": {"title": "Jaccard", "thickness": 15},
        })
        layout = {
            "margin": {"l": 150, "r": 80, "t": 30, "b": 100},
            "paper_bgcolor": "rgba(0,0,0,0)",
            "plot_bgcolor": "#fafafa",
            "font": {"family": "Segoe UI, Arial", "size": 11},
        }
        js_parts.append(
            f"Plotly.newPlot('{div_id}', [{trace}], "
            f"{json.dumps(layout)}, {{responsive: true}});"
        )

    return "\n".join(html_parts), "\n".join(js_parts)


def build_composition_tab(overlap_data: dict) -> tuple[str, str]:
    """
    Stacked bar: exclusive / pairwise-shared / all-common variants per tool per size.
    """
    html_parts = [
        '<div class="section-title">Composição das Variantes por Ferramenta</div>',
        '<p style="color:#555;margin-bottom:16px;font-size:.9rem">'
        'Percentual de variantes exclusivas (apenas nessa ferramenta), '
        'compartilhadas com pelo menos uma outra ferramenta, e presentes em '
        'todas as ferramentas. Médias sobre todas as réplicas.</p>',
    ]
    js_parts = []

    sizes_with_data = sorted(overlap_data.keys())
    div_id = "chart_composition"
    html_parts.append(
        f'<div class="card"><div id="{div_id}" style="height:420px"></div></div>'
    )

    # For each tool: across all sizes, compute mean % exclusive, shared, all-common
    traces = []
    bar_colors_excl   = {sw: SOFT_COLOR.get(sw, "#888") for sw in SOFTWARES}

    for sw in SOFTWARES:
        label = SOFT_LABEL.get(sw, sw)
        color = SOFT_COLOR.get(sw, "#888")

        xs, excl_pcts, allc_pcts = [], [], []
        for size in sizes_with_data:
            rep_dict = overlap_data.get(size, {})
            excl_vals, allc_vals = [], []
            for rep, ov in rep_dict.items():
                tools_present = ov.get("tools", [])
                if sw not in tools_present:
                    continue
                total_sw = ov.get("counts", {}).get(sw, 0)
                if total_sw == 0:
                    continue
                excl = ov.get("exclusive", {}).get(sw, 0)
                allc = ov.get("all_common", 0)
                excl_vals.append(100.0 * excl / total_sw)
                allc_vals.append(100.0 * allc / total_sw)

            if excl_vals:
                xs.append(size)
                excl_pcts.append(round(_mean(excl_vals), 2))
                allc_pcts.append(round(_mean(allc_vals), 2))

        if not xs:
            continue

        traces.append(json.dumps({
            "x": xs, "y": excl_pcts, "name": f"{label} exclusivo",
            "type": "bar", "marker": {"color": color, "opacity": 1.0},
            "hovertemplate": f"{label} exclusivo: %{{y:.1f}}%<extra></extra>",
        }))

    layout = dict(_CHART_LAYOUT_BASE)
    layout["title"]  = {"text": "% Variantes Exclusivas por Ferramenta", "font": {"size": 14}}
    layout["yaxis"]  = {"title": "% de variantes exclusivas", "gridcolor": "#eee"}
    layout["barmode"] = "group"
    layout["legend"] = {"orientation": "h", "y": -0.3}

    js_parts.append(
        f"Plotly.newPlot('{div_id}', [{', '.join(traces)}], "
        f"{json.dumps(layout)}, {{responsive: true}});"
    )

    # All-common chart
    div_id2 = "chart_all_common"
    html_parts.append(
        f'<div class="card"><div class="section-title" style="font-size:1rem">'
        f'Variantes presentes em TODAS as ferramentas</div>'
        f'<div id="{div_id2}" style="height:300px"></div></div>'
    )

    traces2 = []
    for sw in SOFTWARES:
        label = SOFT_LABEL.get(sw, sw)
        color = SOFT_COLOR.get(sw, "#888")
        xs, ys = [], []
        for size in sizes_with_data:
            rep_dict = overlap_data.get(size, {})
            vals = []
            for rep, ov in rep_dict.items():
                tools_present = ov.get("tools", [])
                if sw not in tools_present:
                    continue
                total_sw = ov.get("counts", {}).get(sw, 0)
                if total_sw == 0:
                    continue
                allc = ov.get("all_common", 0)
                vals.append(100.0 * allc / total_sw)
            if vals:
                xs.append(size)
                ys.append(round(_mean(vals), 2))
        if xs:
            traces2.append(json.dumps({
                "x": xs, "y": ys, "name": label,
                "mode": "lines+markers",
                "line": {"color": color, "width": 2},
                "marker": {"size": 7},
                "hovertemplate": f"{label}: %{{y:.1f}}%<extra></extra>",
            }))

    layout2 = dict(_CHART_LAYOUT_BASE)
    layout2["title"] = {"text": "% Variantes Presentes em Todas as Ferramentas", "font": {"size": 14}}
    layout2["yaxis"] = {"title": "% do total da ferramenta", "gridcolor": "#eee"}

    js_parts.append(
        f"Plotly.newPlot('{div_id2}', [{', '.join(traces2)}], "
        f"{json.dumps(layout2)}, {{responsive: true}});"
    )

    return "\n".join(html_parts), "\n".join(js_parts)


def build_stability_tab(jaccard_table: dict) -> tuple[str, str]:
    """
    Cross-replicate concordance stability: for each pair × size,
    show mean ± SD Jaccard across replicates.
    """
    html_parts = [
        '<div class="section-title">Estabilidade da Concordância entre Réplicas</div>',
        '<p style="color:#555;margin-bottom:16px;font-size:.9rem">'
        'Cada ponto mostra o Jaccard médio ± desvio padrão calculado sobre as 3 réplicas. '
        'Barras de erro largas indicam comportamento inconsistente entre execuções '
        'com diferentes conjuntos de amostras.</p>',
    ]
    js_parts = []

    pairs = all_pairs()
    div_id = "chart_stability"
    html_parts.append(f'<div class="card"><div id="{div_id}" style="height:400px"></div></div>')

    traces = []
    for ci, (t1, t2) in enumerate(pairs):
        pair_key = (t1, t2)
        color = PAIR_COLORS[ci % len(PAIR_COLORS)]
        size_dict = jaccard_table.get(pair_key, {})
        sizes = sorted(size_dict.keys())

        xs, ys, es = [], [], []
        for size in sizes:
            vals = list(size_dict[size].values())
            if not vals:
                continue
            xs.append(size)
            ys.append(round(_mean(vals), 4))
            es.append(round(_std(vals), 4))

        if not xs:
            continue

        d = {
            "x": xs, "y": ys,
            "name": _pair_label(t1, t2),
            "mode": "lines+markers",
            "line": {"color": color, "width": 2},
            "marker": {"size": 7},
        }
        if any(e > 0 for e in es):
            d["error_y"] = {
                "type": "data", "array": es, "visible": True,
                "color": color, "thickness": 1.5, "width": 4,
            }
        traces.append(json.dumps(d))

    layout = dict(_CHART_LAYOUT_BASE)
    layout["title"] = {"text": "Estabilidade do Jaccard por tamanho (média ± DP das réplicas)",
                       "font": {"size": 14}}
    layout["yaxis"] = {"title": "Jaccard similarity", "range": [0, 1.05], "gridcolor": "#eee"}

    js_parts.append(
        f"Plotly.newPlot('{div_id}', [{', '.join(traces)}], "
        f"{json.dumps(layout)}, {{responsive: true}});"
    )

    # Stability summary table
    html_parts.append('<div class="section-title">Tabela de Estabilidade</div>')
    html_parts.append('<div class="card"><table>')
    sizes_all = sorted({s for _, sd in jaccard_table.items() for s in sd})
    hdr = ["Par"] + [f"N={s} (média ± DP)" for s in sizes_all]
    html_parts.append("<tr>" + "".join(f"<th>{h}</th>" for h in hdr) + "</tr>")

    for ci, (t1, t2) in enumerate(pairs):
        pair_key = (t1, t2)
        row = [f'<td style="text-align:left;font-weight:600">{_pair_label(t1,t2)}</td>']
        for size in sizes_all:
            vals = list(jaccard_table.get(pair_key, {}).get(size, {}).values())
            if not vals:
                row.append("<td>—</td>")
                continue
            m = _mean(vals)
            s = _std(vals)
            color = (
                "#2ca02c" if m >= WARNING_THRESHOLD else
                "#ff7f0e" if m >= CRITICAL_THRESHOLD else
                "#d62728"
            )
            row.append(
                f'<td style="color:{color}">'
                f'<strong>{m:.3f}</strong> ± {s:.3f}</td>'
            )
        html_parts.append("<tr>" + "".join(row) + "</tr>")

    html_parts.append("</table></div>")

    return "\n".join(html_parts), "\n".join(js_parts)


def build_rawdata_tab(overlap_data: dict, jaccard_table: dict) -> str:
    """Full concordance data table."""
    html_parts = [
        '<div class="section-title">Dados Brutos de Concordância</div>',
        '<div class="card"><table>',
    ]

    pairs = all_pairs()
    header = (
        ["Tamanho", "Réplica", "Ferramentas presentes", "Total únicas"]
        + [f"Jaccard: {_pair_label(t1,t2)}" for t1, t2 in pairs]
        + [f"Exclusivas: {SOFT_LABEL.get(sw,sw)}" for sw in SOFTWARES]
        + ["Em comum (todas)"]
    )
    html_parts.append("<tr>" + "".join(f"<th>{h}</th>" for h in header) + "</tr>")

    for size in sorted(overlap_data.keys()):
        for rep in sorted(overlap_data[size].keys()):
            ov = overlap_data[size][rep]
            tools_present = ov.get("tools", [])
            total_unique  = ov.get("total_unique", "—")
            row = [f"<td>{size}</td>", f"<td>{rep}</td>",
                   f"<td>{', '.join(SOFT_LABEL.get(t,t) for t in tools_present)}</td>",
                   f"<td>{total_unique:,}</td>" if isinstance(total_unique, int) else f"<td>{total_unique}</td>"]

            for t1, t2 in pairs:
                j = None
                if t1 in tools_present and t2 in tools_present:
                    j = jaccard(ov, t1, t2)
                if j is None:
                    row.append("<td>—</td>")
                else:
                    color = (
                        "#2ca02c" if j >= WARNING_THRESHOLD else
                        "#ff7f0e" if j >= CRITICAL_THRESHOLD else
                        "#d62728"
                    )
                    row.append(f'<td style="color:{color}">{j:.4f}</td>')

            for sw in SOFTWARES:
                excl = ov.get("exclusive", {}).get(sw)
                row.append(f"<td>{excl:,}</td>" if isinstance(excl, int) else "<td>—</td>")

            allc = ov.get("all_common", "—")
            row.append(f"<td>{allc:,}</td>" if isinstance(allc, int) else f"<td>{allc}</td>")

            html_parts.append("<tr>" + "".join(row) + "</tr>")

    html_parts.append("</table></div>")
    return "\n".join(html_parts)


# ─── Main HTML generator ──────────────────────────────────────────────────────

def generate_concordance_html(
    overlap_data: dict,
    jaccard_table: dict,
    output_path: Path,
):
    alerts = find_alerts(jaccard_table)

    ov_html,   ov_js   = build_overview_tab(overlap_data, jaccard_table)
    alr_html            = build_alerts_tab(alerts)
    hm_html,   hm_js   = build_heatmap_tab(overlap_data, jaccard_table)
    comp_html, comp_js  = build_composition_tab(overlap_data)
    stab_html, stab_js  = build_stability_tab(jaccard_table)
    raw_html            = build_rawdata_tab(overlap_data, jaccard_table)

    n_alerts_critical = sum(1 for a in alerts if a["level"] == "CRÍTICO")
    n_alerts_warning  = sum(1 for a in alerts if a["level"] == "AVISO")

    alert_badge = ""
    if n_alerts_critical:
        alert_badge = (
            f'<span style="background:#d62728;color:#fff;padding:2px 8px;'
            f'border-radius:10px;font-size:.75rem;margin-left:6px">'
            f'{n_alerts_critical} CRÍTICO</span>'
        )
    elif n_alerts_warning:
        alert_badge = (
            f'<span style="background:#ff7f0e;color:#fff;padding:2px 8px;'
            f'border-radius:10px;font-size:.75rem;margin-left:6px">'
            f'{n_alerts_warning} AVISO</span>'
        )

    html = f"""<!DOCTYPE html>
<html lang="pt-BR">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>GVCF Benchmark — Relatório de Concordância</title>
  <script src="{PLOTLY_CDN}"></script>
  <style>{_CSS}</style>
</head>
<body>

<header>
  <h1>Relatório de Concordância Inter-Ferramenta</h1>
  <p>Gerado: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')} &nbsp;|&nbsp;
     Ferramentas: GLnexus, Parabricks, GATK CombineGVCFs, GATK GenomicsDB &nbsp;|&nbsp;
     Métrica: Jaccard similarity por réplica × tamanho</p>
</header>

<div class="tabs">
  <button class="tab-btn active" onclick="showTab('tab-overview',this)">Visão Geral</button>
  <button class="tab-btn" onclick="showTab('tab-alerts',this)">Alertas{alert_badge}</button>
  <button class="tab-btn" onclick="showTab('tab-heatmap',this)">Heatmap</button>
  <button class="tab-btn" onclick="showTab('tab-composition',this)">Composição</button>
  <button class="tab-btn" onclick="showTab('tab-stability',this)">Estabilidade</button>
  <button class="tab-btn" onclick="showTab('tab-rawdata',this)">Dados Brutos</button>
</div>

<div id="tab-overview" class="tab-content active">
{ov_html}
</div>

<div id="tab-alerts" class="tab-content">
{alr_html}
</div>

<div id="tab-heatmap" class="tab-content">
{hm_html}
</div>

<div id="tab-composition" class="tab-content">
{comp_html}
</div>

<div id="tab-stability" class="tab-content">
{stab_html}
</div>

<div id="tab-rawdata" class="tab-content">
{raw_html}
</div>

<footer style="padding:16px 32px;text-align:center;color:#999;font-size:0.8rem;margin-top:32px">
  GVCF Scalability Benchmark — Concordância Inter-Ferramenta &mdash; {datetime.now().year}
</footer>

<script>
{_JS_TAB}

(function() {{
  {ov_js}
  {hm_js}
  {comp_js}
  {stab_js}
}})();
</script>

</body>
</html>"""

    output_path.write_text(html, encoding="utf-8")
    print(f"[CONCORDANCE] Written: {output_path}")


# ─── Main ─────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(
        description="Generate inter-tool concordance report from variant overlap data"
    )
    ap.add_argument("--benchmark-dir", required=True,
                    help="Root benchmarks directory (contains 03_metrics/)")
    ap.add_argument("--output-dir", required=True,
                    help="Directory to write concordance_report.html")
    args = ap.parse_args()

    metrics_dir = Path(args.benchmark_dir) / "03_metrics"
    output_dir  = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    if not metrics_dir.exists():
        print(f"[ERROR] Metrics directory not found: {metrics_dir}", file=sys.stderr)
        sys.exit(1)

    print("[CONCORDANCE] Loading variant overlap data...")
    overlap_data  = load_overlap_data(metrics_dir)

    if not overlap_data:
        print("[WARN] No variant overlap files found — concordance report will be empty.",
              file=sys.stderr)

    jaccard_table = compute_jaccard_table(overlap_data)
    alerts        = find_alerts(jaccard_table)

    if alerts:
        print(f"[CONCORDANCE] ⚠  {len(alerts)} low-concordance pair(s) detected:")
        for a in alerts:
            print(f"  [{a['level']}] {_pair_label(a['t1'],a['t2'])} "
                  f"N={a['size']} rep{a['rep']}: Jaccard={a['jaccard']:.4f}")
    else:
        print("[CONCORDANCE] All concordance values above threshold.")

    out_path = output_dir / "concordance_report.html"
    generate_concordance_html(overlap_data, jaccard_table, out_path)
    print("[CONCORDANCE] Done.")


if __name__ == "__main__":
    main()
