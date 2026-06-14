#!/usr/bin/env python3
"""
concordance_report.py
=====================
Generate an HTML QC report for the genotype concordance analysis.

Reads the outputs of pairwise_concordance.py and renders:
  - Summary statistics card
  - Concordance distribution histogram (SVG, no external deps)
  - Full concordance heatmap (SVG)
  - Flagged duplicates table
  - Per-sample worst-pair table

Usage:
    python3 concordance_report.py \\
        --matrix   <prefix>_concordance_matrix.csv \\
        --pairs    <prefix>_concordance_pairs.tsv \\
        --flagged  <prefix>_flagged_duplicates.tsv \\
        --summary  <prefix>_concordance_summary.txt \\
        --out      <prefix>_concordance_report.html \\
        --threshold 0.99
"""

import argparse
import csv
import html
import math
import os
import sys
from typing import Dict, List, Optional, Tuple

import numpy as np


# ---------------------------------------------------------------------------
# Readers
# ---------------------------------------------------------------------------

def read_matrix(path: str) -> Tuple[List[str], np.ndarray]:
    with open(path, newline="") as fh:
        reader = csv.reader(fh)
        header = next(reader)
        sample_ids = header[1:]
        rows = []
        for line in reader:
            vals = []
            for v in line[1:]:
                try:
                    vals.append(float(v))
                except (ValueError, TypeError):
                    vals.append(float("nan"))
            rows.append(vals)
    return sample_ids, np.array(rows)


def read_pairs(path: str) -> List[Dict]:
    rows = []
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            try:
                row["concordance"] = float(row["concordance"])
            except (ValueError, TypeError):
                row["concordance"] = float("nan")
            row["n_compared"] = int(row.get("n_compared", 0))
            rows.append(row)
    return rows


def read_flagged(path: str) -> List[Dict]:
    rows = []
    if not os.path.exists(path):
        return rows
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            try:
                row["concordance"] = float(row["concordance"])
            except (ValueError, TypeError):
                row["concordance"] = float("nan")
            rows.append(row)
    return rows


def read_summary(path: str) -> str:
    if not os.path.exists(path):
        return ""
    with open(path) as fh:
        return fh.read()


# ---------------------------------------------------------------------------
# SVG helpers
# ---------------------------------------------------------------------------

def _colour_scale(val: float, vmin: float = 0.4, vmax: float = 1.0) -> str:
    """Map a concordance value to an RGB colour (blue-white-red)."""
    if math.isnan(val):
        return "#cccccc"
    t = max(0.0, min(1.0, (val - vmin) / max(vmax - vmin, 1e-9)))
    # Low → blue, mid → white, high → red
    if t < 0.5:
        s = t * 2
        r = int(255 * s)
        g = int(255 * s)
        b = 255
    else:
        s = (t - 0.5) * 2
        r = 255
        g = int(255 * (1 - s))
        b = int(255 * (1 - s))
    return f"#{r:02x}{g:02x}{b:02x}"


def make_heatmap_svg(
    sample_ids: List[str],
    mat: np.ndarray,
    threshold: float,
    max_dim: int = 600,
) -> str:
    n = len(sample_ids)
    cell = max(4, min(20, max_dim // max(n, 1)))
    margin = 140
    size = n * cell + margin

    lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" '
        f'width="{size}" height="{size}" '
        f'style="font-family:monospace;font-size:10px">',
    ]

    # Row labels
    for i, sid in enumerate(sample_ids):
        y = margin + i * cell + cell // 2 + 4
        label = html.escape(sid[-24:]) if len(sid) > 24 else html.escape(sid)
        lines.append(
            f'<text x="{margin - 4}" y="{y}" '
            f'text-anchor="end" font-size="9">{label}</text>'
        )

    # Column labels (rotated)
    for j, sid in enumerate(sample_ids):
        x = margin + j * cell + cell // 2
        label = html.escape(sid[-24:]) if len(sid) > 24 else html.escape(sid)
        lines.append(
            f'<text x="{x}" y="{margin - 4}" '
            f'text-anchor="start" font-size="9" '
            f'transform="rotate(-60 {x} {margin - 4})">{label}</text>'
        )

    # Cells
    for i in range(n):
        for j in range(n):
            val = mat[i, j]
            colour = _colour_scale(val)
            x = margin + j * cell
            y = margin + i * cell
            stroke = "red" if (
                i != j
                and not math.isnan(val)
                and val >= threshold
            ) else "none"
            sw = "1.5" if stroke != "none" else "0"
            tip = f"{sample_ids[i]} × {sample_ids[j]}: {val:.4f}" if not math.isnan(val) else "N/A"
            lines.append(
                f'<rect x="{x}" y="{y}" width="{cell}" height="{cell}" '
                f'fill="{colour}" stroke="{stroke}" stroke-width="{sw}">'
                f'<title>{html.escape(tip)}</title></rect>'
            )

    lines.append("</svg>")
    return "\n".join(lines)


def make_histogram_svg(
    concordances: List[float],
    threshold: float,
    width: int = 560,
    height: int = 260,
) -> str:
    if not concordances:
        return "<p><em>No reliable pairs to plot.</em></p>"

    bins = 40
    arr = np.array(concordances)
    counts, edges = np.histogram(arr, bins=bins, range=(0.0, 1.0))
    max_count = max(counts) if counts.max() > 0 else 1

    pad_l, pad_r, pad_t, pad_b = 50, 20, 20, 40
    inner_w = width - pad_l - pad_r
    inner_h = height - pad_t - pad_b
    bar_w = inner_w / bins

    lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" '
        f'width="{width}" height="{height}" '
        f'style="font-family:sans-serif;font-size:11px">',
    ]

    # Axes
    lines.append(
        f'<line x1="{pad_l}" y1="{pad_t}" x2="{pad_l}" '
        f'y2="{pad_t + inner_h}" stroke="#333" stroke-width="1"/>'
    )
    lines.append(
        f'<line x1="{pad_l}" y1="{pad_t + inner_h}" '
        f'x2="{pad_l + inner_w}" y2="{pad_t + inner_h}" '
        f'stroke="#333" stroke-width="1"/>'
    )

    # Bars
    for k, cnt in enumerate(counts):
        x = pad_l + k * bar_w
        bar_h = (cnt / max_count) * inner_h
        y = pad_t + inner_h - bar_h
        mid = (edges[k] + edges[k + 1]) / 2
        colour = "#c0392b" if mid >= threshold else "#3498db"
        lines.append(
            f'<rect x="{x:.1f}" y="{y:.1f}" '
            f'width="{bar_w:.1f}" height="{bar_h:.1f}" '
            f'fill="{colour}" opacity="0.8"/>'
        )

    # Threshold line
    tx = pad_l + threshold * inner_w
    lines.append(
        f'<line x1="{tx:.1f}" y1="{pad_t}" '
        f'x2="{tx:.1f}" y2="{pad_t + inner_h}" '
        f'stroke="red" stroke-width="1.5" stroke-dasharray="5,3"/>'
    )
    lines.append(
        f'<text x="{tx + 4:.1f}" y="{pad_t + 14}" '
        f'fill="red" font-size="10">threshold={threshold:.2f}</text>'
    )

    # X-axis ticks
    for v in np.arange(0, 1.01, 0.2):
        x = pad_l + v * inner_w
        y_axis = pad_t + inner_h
        lines.append(
            f'<line x1="{x:.1f}" y1="{y_axis}" '
            f'x2="{x:.1f}" y2="{y_axis + 4}" stroke="#333"/>'
        )
        lines.append(
            f'<text x="{x:.1f}" y="{y_axis + 16}" '
            f'text-anchor="middle">{v:.1f}</text>'
        )

    # Y-axis ticks (just 0 and max)
    for frac, label in [(0.0, "0"), (0.5, str(max_count // 2)), (1.0, str(max_count))]:
        y = pad_t + inner_h - frac * inner_h
        lines.append(
            f'<line x1="{pad_l - 4}" y1="{y:.1f}" '
            f'x2="{pad_l}" y2="{y:.1f}" stroke="#333"/>'
        )
        lines.append(
            f'<text x="{pad_l - 8}" y="{y + 4:.1f}" '
            f'text-anchor="end">{label}</text>'
        )

    # Axis labels
    lines.append(
        f'<text x="{pad_l + inner_w / 2:.1f}" y="{height - 4}" '
        f'text-anchor="middle">Concordance</text>'
    )
    lines.append(
        f'<text x="12" y="{pad_t + inner_h / 2:.1f}" '
        f'text-anchor="middle" '
        f'transform="rotate(-90 12 {pad_t + inner_h / 2:.1f})">Pairs</text>'
    )

    lines.append("</svg>")
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# HTML assembly
# ---------------------------------------------------------------------------

CSS = """
body{font-family:'Segoe UI',Tahoma,Geneva,Verdana,sans-serif;
     margin:0;padding:0;background:#f5f5f5}
.container{max-width:1300px;margin:0 auto;padding:24px}
.header{background:linear-gradient(135deg,#1a6b8a,#2ecc71);
        color:#fff;padding:28px;border-radius:10px;margin-bottom:24px}
h1{margin:0;font-size:2em}
.subtitle{opacity:.9;margin-top:8px}
.card{background:#fff;padding:24px;border-radius:10px;
      box-shadow:0 2px 8px rgba(0,0,0,.08);margin-bottom:20px}
h2{color:#1a6b8a;border-bottom:3px solid #1a6b8a;padding-bottom:8px}
.stat-grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(180px,1fr));gap:16px;margin:16px 0}
.stat-box{background:linear-gradient(135deg,#1a6b8a,#2ecc71);
          color:#fff;padding:18px;border-radius:8px;text-align:center}
.stat-value{font-size:2.2em;font-weight:bold;margin:8px 0}
.stat-label{opacity:.9;font-size:.85em}
table{width:100%;border-collapse:collapse;margin:16px 0;font-size:.92em}
th{background:#1a6b8a;color:#fff;padding:10px 12px;text-align:left}
td{padding:9px 12px;border-bottom:1px solid #e8e8e8}
tr:hover td{background:#f0f7fa}
.flag{background:#fdecea;color:#c0392b;font-weight:bold;padding:2px 8px;border-radius:4px}
.ok{color:#27ae60}
pre{background:#f4f4f4;padding:16px;border-radius:6px;
    font-size:.85em;overflow-x:auto;white-space:pre-wrap}
.viz{overflow-x:auto;margin:16px 0}
.warn{background:#fff3cd;border-left:4px solid #f39c12;
      padding:14px;border-radius:4px;margin:12px 0}
.info{background:#e3f2fd;border-left:4px solid #2196f3;
      padding:14px;border-radius:4px;margin:12px 0}
"""


def _td(val, cls=""):
    c = f' class="{cls}"' if cls else ""
    return f"<td{c}>{html.escape(str(val))}</td>"


def build_html(
    sample_ids: List[str],
    mat: np.ndarray,
    pairs: List[Dict],
    flagged: List[Dict],
    summary_text: str,
    threshold: float,
    n_snps: int,
) -> str:
    reliable = [
        p for p in pairs
        if p.get("reliable", "True") in (True, "True")
        and not math.isnan(p["concordance"])
    ]
    concordances = [p["concordance"] for p in reliable]

    n_samples = len(sample_ids)
    n_pairs   = len(pairs)
    n_flagged = len(flagged)
    mean_conc = f"{np.mean(concordances):.4f}" if concordances else "N/A"

    heatmap_svg   = make_heatmap_svg(sample_ids, mat, threshold)
    histogram_svg = make_histogram_svg(concordances, threshold)

    # Per-sample max concordance with any other sample
    idx = {s: i for i, s in enumerate(sample_ids)}
    per_sample_max: Dict[str, Tuple[float, str]] = {}
    for p in reliable:
        a, b, c = p["sample_a"], p["sample_b"], p["concordance"]
        if a not in per_sample_max or c > per_sample_max[a][0]:
            per_sample_max[a] = (c, b)
        if b not in per_sample_max or c > per_sample_max[b][0]:
            per_sample_max[b] = (c, a)

    # ---- assemble HTML ----
    h = f"""<!DOCTYPE html>
<html lang="en"><head><meta charset="UTF-8">
<title>Genotype Concordance Report</title>
<style>{CSS}</style>
</head><body>
<div class="container">

<div class="header">
  <h1>🧬 Genotype Concordance Report</h1>
  <div class="subtitle">
    EPIC Array rs-probe Identity Check — {n_samples} samples, {n_snps} SNP probes
  </div>
</div>

<div class="card">
  <h2>📊 Summary</h2>
  <div class="stat-grid">
    <div class="stat-box"><div class="stat-value">{n_samples}</div>
      <div class="stat-label">Samples</div></div>
    <div class="stat-box"><div class="stat-value">{n_snps}</div>
      <div class="stat-label">SNP Probes</div></div>
    <div class="stat-box"><div class="stat-value">{n_pairs}</div>
      <div class="stat-label">Sample Pairs</div></div>
    <div class="stat-box"><div class="stat-value">{mean_conc}</div>
      <div class="stat-label">Mean Concordance</div></div>
    <div class="stat-box"><div class="stat-value {'style="color:#e74c3c"' if n_flagged else ''}">
      {n_flagged}</div>
      <div class="stat-label">Flagged Duplicates</div></div>
  </div>
  <div class="info">
    <strong>Method:</strong> Concordance is computed from the {n_snps} rs-prefixed probes
    on the Illumina EPIC array. These probes target biallelic SNP positions and are
    not subject to bisulfite-induced signal perturbation. Genotype calls are
    discretised from beta values: AA&nbsp;(β&lt;0.25)&nbsp;=&nbsp;0,
    AB&nbsp;(0.25≤β≤0.75)&nbsp;=&nbsp;1, BB&nbsp;(β&gt;0.75)&nbsp;=&nbsp;2.
    Pairs with concordance ≥ {threshold} are flagged as likely duplicates or sample swaps.
  </div>
</div>
"""

    # Flagged duplicates
    if flagged:
        h += f"""
<div class="card">
  <h2>⚠️ Flagged Duplicates (concordance ≥ {threshold})</h2>
  <div class="warn">
    {n_flagged} pair(s) exceed the duplicate threshold and require investigation.
    These may be genuine technical replicates, sample swaps, or mislabelled aliquots.
  </div>
  <table>
    <tr><th>Sample A</th><th>Sample B</th>
        <th>SNPs Compared</th><th>SNPs Concordant</th>
        <th>Concordance</th><th>Flag</th></tr>
"""
        for r in sorted(flagged, key=lambda x: -x["concordance"]):
            conc_str = f"{float(r['concordance']):.4f}"
            h += (
                f"<tr>{_td(r['sample_a'])}{_td(r['sample_b'])}"
                f"{_td(r['n_compared'])}{_td(r['n_concordant'])}"
                f"{_td(conc_str)}"
                f"<td><span class='flag'>LIKELY DUPLICATE</span></td></tr>\n"
            )
        h += "</table></div>\n"
    else:
        h += f"""
<div class="card">
  <h2>✅ No Duplicates Detected</h2>
  <p class="ok">No sample pairs exceeded the concordance threshold of {threshold}.
  All samples appear to be genetically distinct.</p>
</div>
"""

    # Distribution histogram
    h += f"""
<div class="card">
  <h2>📈 Concordance Distribution</h2>
  <p>Distribution of pairwise concordance values across {len(reliable)} reliable pairs.
     Bars in red exceed the duplicate threshold ({threshold}).
     Unrelated samples typically cluster around 0.50–0.65.</p>
  <div class="viz">{histogram_svg}</div>
</div>
"""

    # Heatmap
    h += f"""
<div class="card">
  <h2>🔥 Concordance Heatmap</h2>
  <p>Colour scale: blue = low concordance, red = high concordance.
     Red-outlined cells exceed the duplicate threshold.
     Diagonal is always 1.0 (self-comparison).</p>
  <div class="viz">{heatmap_svg}</div>
</div>
"""

    # Per-sample worst/best table
    h += """
<div class="card">
  <h2>🔍 Per-Sample Highest Concordance Partner</h2>
  <p>For each sample, the highest concordance with any other sample is shown.
     Values close to 1.0 indicate a likely duplicate.</p>
  <table>
    <tr><th>Sample</th><th>Best Match</th><th>Concordance</th><th>Status</th></tr>
"""
    for sid in sorted(sample_ids):
        if sid in per_sample_max:
            best_c, best_s = per_sample_max[sid]
            status = (
                "<span class='flag'>CHECK</span>"
                if best_c >= threshold
                else "<span class='ok'>OK</span>"
            )
            h += (
                f"<tr>{_td(sid)}{_td(best_s)}"
                f"{_td(f'{best_c:.4f}')}<td>{status}</td></tr>\n"
            )
        else:
            h += f"<tr>{_td(sid)}<td>N/A</td><td>N/A</td><td>N/A</td></tr>\n"
    h += "</table></div>\n"

    # Raw summary text
    h += f"""
<div class="card">
  <h2>📋 Full Summary Log</h2>
  <pre>{html.escape(summary_text)}</pre>
</div>

</div></body></html>
"""
    return h


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Generate HTML concordance report"
    )
    p.add_argument("--matrix",    required=True)
    p.add_argument("--pairs",     required=True)
    p.add_argument("--flagged",   required=True)
    p.add_argument("--summary",   required=True)
    p.add_argument("--out",       required=True, help="Output HTML file path")
    p.add_argument("--threshold", type=float, default=0.99)
    p.add_argument("--n_snps",    type=int,   default=65)
    return p.parse_args()


def main() -> None:
    args = parse_args()

    sample_ids, mat = read_matrix(args.matrix)
    pairs            = read_pairs(args.pairs)
    flagged          = read_flagged(args.flagged)
    summary_text     = read_summary(args.summary)

    html_doc = build_html(
        sample_ids=sample_ids,
        mat=mat,
        pairs=pairs,
        flagged=flagged,
        summary_text=summary_text,
        threshold=args.threshold,
        n_snps=args.n_snps,
    )

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, "w", encoding="utf-8") as fh:
        fh.write(html_doc)

    print(f"✓ Concordance HTML report written: {args.out}")


if __name__ == "__main__":
    main()
