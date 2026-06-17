#!/usr/bin/env python3
"""
cross_study_concordance.py
==========================
All EPIC-methylation-samples vs all H3Africa-genotype-samples identity check,
using the 65 rs-probe SNPs as the common DNA fingerprint.

METHOD
------
Allele-call concordance on binarised genotype calls (0/1/2).
  - EPIC side : integer calls from extract_snp_betas.r (.traw, col 7+)
  - H3A side  : integer calls from PLINK .traw via match_h3a_snps.r

No allele-flip correction is applied.  The H3A .traw is produced from raw
PLINK .bed/.bim/.fam files (not imputed data), so PLINK allele coding is
consistent and flip issues do not arise.

ASSIGNED vs BEST FRAMEWORK
----------------------------
When an identity map is provided (--id_map), for each EPIC sample the script
records TWO values:
  assigned_score – concordance with the H3A sample that the samplesheet says
                   belongs to this individual
  best_score     – concordance with the best-matching H3A sample across the
                   entire H3A cohort (which is larger than the EPIC cohort)

ID MAP FORMAT
-------------
A CSV file with exactly two columns and a header row:
    h3a_id,epic_id
    G007,G007
    G009,G009
    ...

Both IDs come from the Sample_ID column of the respective samplesheets
(H3Aflow and EPICflow use the same format; Sample_ID values match).
The H3A cohort is larger, so not every H3A ID will have an EPIC counterpart.

FOUR-CATEGORY CLASSIFICATION
------------------------------
  correct_assignment      – best match IS the assigned H3A ID; score ≥ 0.90
  swap_high_confidence    – best match is a DIFFERENT H3A ID; score ≥ 0.90
  swap_moderate_confidence– best match is different; score 0.70–0.90
  outside_dataset         – best match is an H3A ID that has no corresponding
                            EPIC sample (genotype-only participant); suggests
                            the EPIC sample's DNA belongs to someone who was
                            genotyped but never enrolled for methylation
  low_confidence          – best score < 0.70; no reliable match found

Without --id_map, only best-match is reported (no classification).

Usage:
    python3 cross_study_concordance.py \\
        --epic_traw  concordance/methylation_snps_snp_calls.traw \\
        --h3a_traw   cross_study/h3a_snps_h3a_snp_calls.traw \\
        --out        cross_study \\
        --prefix     cross_study \\
        [--id_map    samplesheets/id_bridge.csv] \\
        [--threshold 0.80] \\
        [--min_snps  20]

Outputs:
    <prefix>_concordance_matrix.csv      – EPIC(rows) × H3A(cols) matrix
    <prefix>_concordance_pairs.tsv       – all pairs, long form
    <prefix>_per_sample_verdict.tsv      – one row per EPIC sample
    <prefix>_classification_summary.tsv  – category counts
    <prefix>_concordance_summary.txt     – plain-text run summary
    <prefix>_concordance_report.html     – self-contained HTML report
"""

import argparse
import csv
import logging
import math
import os
import sys
from collections import defaultdict
from typing import Dict, List, Optional, Tuple

import numpy as np

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)

THRESH_HIGH = 0.90
THRESH_MOD  = 0.70


# =============================================================================
# I/O helpers
# =============================================================================

def parse_traw(path: str) -> Tuple[List[str], List[str], np.ndarray]:
    """
    Parse a PLINK .traw file.
    Returns (snp_ids, sample_ids, int8 matrix[n_snps, n_samples]; -1=missing).
    """
    log.info("Parsing .traw: %s", path)
    snp_ids: List[str] = []
    sample_ids: List[str] = []
    rows = []

    with open(path, "r", newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader)
        sample_ids = header[6:]
        for line in reader:
            if not line:
                continue
            snp_ids.append(line[1])
            calls = []
            for val in line[6:]:
                v = val.strip()
                if v in ("NA", "nan", "", "."):
                    calls.append(-1)
                else:
                    try:
                        calls.append(int(float(v)))
                    except ValueError:
                        calls.append(-1)
            rows.append(calls)

    gt = np.array(rows, dtype=np.int8)
    log.info("  %d SNPs × %d samples", gt.shape[0], gt.shape[1])
    return snp_ids, sample_ids, gt


def parse_id_map(path: str) -> Tuple[Dict[str, str], set]:
    """
    Parse the id_map CSV: h3a_id,epic_id (header required).
    Returns:
      epic_to_h3a : dict  epic_id -> h3a_id
      epic_ids_in_map : set of EPIC IDs that appear in the map
    """
    log.info("Parsing ID map: %s", path)
    epic_to_h3a: Dict[str, str] = {}
    with open(path, "r", newline="") as fh:
        reader = csv.DictReader(fh)
        # Normalise column names to lowercase stripped
        for row in reader:
            row_norm = {k.strip().lower(): v.strip() for k, v in row.items()}
            h3a_id  = row_norm.get("h3a_id")  or row_norm.get("h3a")
            epic_id = row_norm.get("epic_id") or row_norm.get("epic")
            if h3a_id and epic_id:
                epic_to_h3a[epic_id] = h3a_id
    log.info("  %d epic→h3a mappings loaded", len(epic_to_h3a))
    return epic_to_h3a, set(epic_to_h3a.keys())


# =============================================================================
# SNP intersection
# =============================================================================

def intersect_snps(
    snps_a: List[str], mat_a: np.ndarray,
    snps_b: List[str], mat_b: np.ndarray,
) -> Tuple[List[str], np.ndarray, np.ndarray]:
    set_b = {s: i for i, s in enumerate(snps_b)}
    common = [s for s in snps_a if s in set_b]
    if not common:
        raise ValueError(
            "No SNP IDs in common between the two .traw files. "
            "Check that match_h3a_snps.r relabelled H3A variants to EPIC probe names."
        )
    idx_a = [snps_a.index(s) for s in common]
    idx_b = [set_b[s]        for s in common]
    return common, mat_a[idx_a, :], mat_b[idx_b, :]


# =============================================================================
# Concordance computation
# =============================================================================

def compute_concordance(
    mat_epic: np.ndarray,
    samples_epic: List[str],
    mat_h3a: np.ndarray,
    samples_h3a: List[str],
    min_snps: int,
) -> List[Dict]:
    n_epic = mat_epic.shape[1]
    n_h3a  = mat_h3a.shape[1]
    log.info("Concordance: %d EPIC × %d H3A = %d pairs", n_epic, n_h3a, n_epic * n_h3a)

    results = []
    for i, esid in enumerate(samples_epic):
        col_e = mat_epic[:, i]
        for j, hsid in enumerate(samples_h3a):
            col_h  = mat_h3a[:, j]
            valid  = (col_e >= 0) & (col_h >= 0)
            n_comp = int(valid.sum())
            if n_comp == 0:
                score, n_conc = float("nan"), 0
            else:
                n_conc = int((col_e[valid] == col_h[valid]).sum())
                score  = n_conc / n_comp
            results.append({
                "epic_sample":  esid,
                "h3a_sample":   hsid,
                "n_compared":   n_comp,
                "n_concordant": n_conc,
                "concordance":  score,
                "reliable":     n_comp >= min_snps,
            })
    return results


# =============================================================================
# Assigned-vs-best verdict and classification
# =============================================================================

def build_verdicts(
    results:      List[Dict],
    samples_epic: List[str],
    samples_h3a:  List[str],
    epic_to_h3a:  Optional[Dict[str, str]],
) -> List[Dict]:
    """
    For each EPIC sample:
      - Find the best-matching H3A sample (highest concordance, reliable pairs only)
      - If epic_to_h3a is given, also retrieve the assigned H3A sample's score
      - Classify into one of the five categories

    H3A samples whose ID does NOT appear as a value in epic_to_h3a are
    'outside dataset' candidates: they exist in the H3A genotyping cohort but
    were never enrolled in the methylation study.
    """
    enrolled_h3a: set = set(epic_to_h3a.values()) if epic_to_h3a else set()

    by_epic: Dict[str, List[Dict]] = defaultdict(list)
    for r in results:
        if r["reliable"] and not (
            isinstance(r["concordance"], float) and math.isnan(r["concordance"])
        ):
            by_epic[r["epic_sample"]].append(r)

    verdicts = []
    for esid in samples_epic:
        pairs = by_epic.get(esid, [])

        if not pairs:
            verdicts.append({
                "epic_sample":     esid,
                "assigned_h3a":    epic_to_h3a.get(esid, "N/A") if epic_to_h3a else "N/A",
                "assigned_score":  "NA",
                "best_h3a":        "NA",
                "best_score":      "NA",
                "n_snps_best":     0,
                "classification":  "no_reliable_pairs",
            })
            continue

        best = max(pairs, key=lambda x: x["concordance"])
        best_score = best["concordance"]
        best_h3a   = best["h3a_sample"]

        assigned_h3a   = epic_to_h3a.get(esid) if epic_to_h3a else None
        assigned_score = float("nan")
        if assigned_h3a:
            assigned_pairs = [p for p in pairs if p["h3a_sample"] == assigned_h3a]
            if assigned_pairs:
                assigned_score = assigned_pairs[0]["concordance"]

        # Classification
        if epic_to_h3a is None:
            cat = "unclassified_no_id_map"
        elif math.isnan(best_score):
            cat = "low_confidence"
        elif assigned_h3a and best_h3a == assigned_h3a and best_score >= THRESH_HIGH:
            cat = "correct_assignment"
        elif best_h3a not in enrolled_h3a:
            # Best match is a genotype-only participant — never enrolled for methylation
            cat = "outside_dataset"
        elif best_score >= THRESH_HIGH:
            cat = "swap_high_confidence"
        elif best_score >= THRESH_MOD:
            cat = "swap_moderate_confidence"
        else:
            cat = "low_confidence"

        verdicts.append({
            "epic_sample":    esid,
            "assigned_h3a":   assigned_h3a if assigned_h3a else "N/A",
            "assigned_score": round(assigned_score, 4) if not math.isnan(assigned_score) else "NA",
            "best_h3a":       best_h3a,
            "best_score":     round(best_score, 4),
            "n_snps_best":    best["n_compared"],
            "classification": cat,
        })

    return verdicts


# =============================================================================
# Output writers
# =============================================================================

def write_matrix(
    results: List[Dict],
    samples_epic: List[str],
    samples_h3a: List[str],
    path: str,
) -> None:
    idx_e = {s: i for i, s in enumerate(samples_epic)}
    idx_h = {s: i for i, s in enumerate(samples_h3a)}
    mat = np.full((len(samples_epic), len(samples_h3a)), np.nan)
    for r in results:
        i = idx_e.get(r["epic_sample"])
        j = idx_h.get(r["h3a_sample"])
        if i is not None and j is not None:
            v = r["concordance"]
            if isinstance(v, float) and not math.isnan(v):
                mat[i, j] = v
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["EPIC_sample \\ H3A_sample"] + samples_h3a)
        for i, sid in enumerate(samples_epic):
            row = ["" if np.isnan(mat[i, j]) else f"{mat[i, j]:.6f}"
                   for j in range(len(samples_h3a))]
            w.writerow([sid] + row)
    log.info("Matrix written: %s", path)


def write_pairs(results: List[Dict], path: str) -> None:
    fields = ["epic_sample", "h3a_sample", "n_compared",
              "n_concordant", "concordance", "reliable"]
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        for r in results:
            row = dict(r)
            if isinstance(row["concordance"], float) and not math.isnan(row["concordance"]):
                row["concordance"] = f"{row['concordance']:.6f}"
            w.writerow(row)
    log.info("Pairs written: %s  (%d pairs)", path, len(results))


def write_verdicts(verdicts: List[Dict], path: str) -> None:
    if not verdicts:
        return
    fields = ["epic_sample", "assigned_h3a", "assigned_score",
              "best_h3a", "best_score", "n_snps_best", "classification"]
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(verdicts)
    log.info("Verdicts written: %s", path)


def write_classification_summary(verdicts: List[Dict], path: str) -> Dict[str, int]:
    counts: Dict[str, int] = defaultdict(int)
    for v in verdicts:
        counts[v["classification"]] += 1

    order = [
        "correct_assignment",
        "outside_dataset",
        "swap_high_confidence",
        "swap_moderate_confidence",
        "low_confidence",
        "no_reliable_pairs",
        "unclassified_no_id_map",
    ]
    sorted_counts = {k: counts[k] for k in order if k in counts}
    sorted_counts.update({k: v for k, v in counts.items() if k not in order})

    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["classification", "count"])
        for k, v in sorted_counts.items():
            w.writerow([k, v])

    lines = ["Classification summary:"]
    icons = {"correct_assignment": "✓", "swap_high_confidence": "⚠",
             "swap_moderate_confidence": "⚠", "outside_dataset": "🔍",
             "low_confidence": "✗"}
    for k, v in sorted_counts.items():
        lines.append(f"  {icons.get(k,' ')} {k:<40s}: {v}")
    print("\n".join(lines))
    log.info("Classification summary written: %s", path)
    return dict(sorted_counts)


def write_summary_txt(
    results: List[Dict],
    verdicts: List[Dict],
    class_counts: Dict[str, int],
    common_snps: List[str],
    n_epic: int,
    n_h3a: int,
    threshold: float,
    min_snps: int,
    path: str,
) -> None:
    reliable = [r for r in results
                if r["reliable"]
                and isinstance(r["concordance"], float)
                and not math.isnan(r["concordance"])]
    concs = np.array([r["concordance"] for r in reliable]) if reliable else np.array([])

    lines = [
        "=" * 65,
        "CROSS-STUDY IDENTITY CHECK SUMMARY",
        "  EPIC methylation array vs H3Africa genotyping array",
        "  Method: allele-call concordance",
        "=" * 65,
        f"  EPIC samples        : {n_epic}",
        f"  H3A samples         : {n_h3a}",
        f"  SNPs in common      : {len(common_snps)}",
        f"  Min SNPs (reliable) : {min_snps}",
        f"  Report threshold    : {threshold}",
        f"  Reliable pairs      : {len(reliable)}",
        f"  Unreliable pairs    : {n_epic * n_h3a - len(reliable)}",
        "",
        "  Concordance distribution (reliable pairs):",
    ]

    if len(concs) > 0:
        lines += [
            f"    Min    : {concs.min():.4f}",
            f"    Median : {np.median(concs):.4f}",
            f"    Mean   : {concs.mean():.4f}",
            f"    Max    : {concs.max():.4f}",
            f"    Std    : {concs.std():.4f}",
            "",
            "  Concordance buckets:",
        ]
        for lo, hi, label in [
            (0.00, 0.60, "unrelated      (0.00–0.60)"),
            (0.60, 0.75, "low agreement  (0.60–0.75)"),
            (0.75, 0.90, "likely match   (0.75–0.90)"),
            (0.90, 1.01, "strong match   (≥ 0.90)   "),
        ]:
            count = int(((concs >= lo) & (concs < hi)).sum())
            lines.append(f"    {label}: {count}")
        lines.append("")

    if class_counts:
        lines.append("  Classification (per EPIC sample):")
        icons = {"correct_assignment": "✓", "swap_high_confidence": "⚠",
                 "swap_moderate_confidence": "⚠", "outside_dataset": "🔍",
                 "low_confidence": "✗"}
        for k, v in class_counts.items():
            lines.append(f"    {icons.get(k,' ')} {k:<40s}: {v}")
        lines.append("")

    lines.append("=" * 65)

    text = "\n".join(lines)
    with open(path, "w") as fh:
        fh.write(text + "\n")
    print("\n" + text)
    log.info("Summary written: %s", path)


# =============================================================================
# HTML report
# =============================================================================

def _hist_svg(scores: List[float], threshold: float, width=560, height=160) -> str:
    n_bins, bin_size = 20, 0.05
    hist = [0] * n_bins
    for c in scores:
        if not math.isnan(c):
            b = min(int(c / bin_size), n_bins - 1)
            hist[b] += 1
    max_h  = max(hist) if hist else 1
    bar_w  = (width - 50) // n_bins
    chart_h = height - 40
    bars = ""
    for i, cnt in enumerate(hist):
        bh = int(cnt / max_h * chart_h) if max_h > 0 else 0
        x  = 40 + i * bar_w
        y  = 20 + (chart_h - bh)
        lo = i * bin_size
        color = "#e74c3c" if lo >= threshold else "#3498db"
        bars += (f'<rect x="{x}" y="{y}" width="{bar_w-1}" height="{bh}" '
                 f'fill="{color}" opacity="0.8">'
                 f'<title>{cnt} pairs [{lo:.2f}–{lo+bin_size:.2f})</title></rect>\n')
    for i in range(0, n_bins + 1, 4):
        x = 40 + i * bar_w
        bars += f'<text x="{x}" y="{20+chart_h+14}" font-size="9" text-anchor="middle">{i*bin_size:.1f}</text>\n'
    bars += (f'<line x1="40" y1="20" x2="40" y2="{20+chart_h}" stroke="#999" stroke-width="1"/>'
             f'<line x1="40" y1="{20+chart_h}" x2="{width-10}" y2="{20+chart_h}" stroke="#999" stroke-width="1"/>')
    tx = int(40 + threshold / bin_size * bar_w)
    bars += (f'<line x1="{tx}" y1="20" x2="{tx}" y2="{20+chart_h}" stroke="#e74c3c" '
             f'stroke-width="1.5" stroke-dasharray="4"/>'
             f'<text x="{tx+3}" y="30" font-size="9" fill="#e74c3c">threshold</text>')
    bars += f'<text x="{width//2}" y="{height-4}" font-size="10" text-anchor="middle">Allele-call concordance</text>'
    return f'<svg width="{width}" height="{height}" style="display:block;margin:8px auto">{bars}</svg>'


def _cat_color(cat: str) -> str:
    return {
        "correct_assignment":         "#27ae60",
        "swap_high_confidence":       "#e67e22",
        "swap_moderate_confidence":   "#f39c12",
        "outside_dataset":            "#8e44ad",
        "low_confidence":             "#e74c3c",
        "no_reliable_pairs":          "#95a5a6",
        "unclassified_no_id_map":     "#7f8c8d",
    }.get(cat, "#2c3e50")


def write_html_report(
    results: List[Dict],
    verdicts: List[Dict],
    class_counts: Dict[str, int],
    common_snps: List[str],
    n_epic: int,
    n_h3a: int,
    threshold: float,
    min_snps: int,
    path: str,
) -> None:
    reliable_scores = [r["concordance"] for r in results
                       if r["reliable"]
                       and isinstance(r["concordance"], float)
                       and not math.isnan(r["concordance"])]
    hist_svg = _hist_svg(reliable_scores, threshold) if reliable_scores else ""

    cat_cards = "".join(
        f"<div class='card' style='border-left-color:{_cat_color(cat)}'>"
        f"<h3>{cat.replace('_',' ')}</h3><p>{cnt}</p></div>\n"
        for cat, cnt in class_counts.items()
    )

    def _verdict_row(v):
        color = _cat_color(v["classification"])
        return (
            f"<tr>"
            f"<td>{v['epic_sample']}</td>"
            f"<td>{v['assigned_h3a']}</td>"
            f"<td>{v.get('assigned_score','NA')}</td>"
            f"<td>{v['best_h3a']}</td>"
            f"<td>{v['best_score']}</td>"
            f"<td>{v['n_snps_best']}</td>"
            f"<td><span style='color:{color};font-weight:bold'>"
            f"{v['classification']}</span></td>"
            f"</tr>\n"
        )
    verdict_rows = "".join(_verdict_row(v) for v in verdicts)

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>Cross-Study Identity Check Report</title>
<style>
body  {{ font-family: Arial, sans-serif; margin: 30px; color: #333; }}
h1    {{ color: #2c3e50; border-bottom: 2px solid #3498db; padding-bottom: 8px; }}
h2    {{ color: #2980b9; margin-top: 30px; }}
.cards {{ display: flex; flex-wrap: wrap; gap: 14px; margin: 18px 0; }}
.card  {{ background: #f8f9fa; border-left: 4px solid #3498db; padding: 10px 16px;
          border-radius: 4px; min-width: 150px; }}
.card h3 {{ margin: 0 0 4px; font-size: 11px; color: #666; text-transform: uppercase; }}
.card p  {{ margin: 0; font-size: 22px; font-weight: bold; color: #2c3e50; }}
table {{ border-collapse: collapse; width: 100%; margin-top: 10px; font-size: 12px; }}
th    {{ background: #2980b9; color: white; padding: 7px 10px; text-align: left; }}
td    {{ padding: 5px 10px; border-bottom: 1px solid #eee; }}
tr:hover td {{ background: #f0f7ff; }}
.note {{ background:#fff3cd; border:1px solid #ffc107; padding:10px 14px;
         border-radius:4px; font-size:13px; margin:12px 0; }}
</style>
</head>
<body>
<h1>Cross-Study Identity Check Report</h1>
<p style="color:#666">EPIC methylation array SNP fingerprint vs H3Africa genotyping array<br>
Method: allele-call concordance (genotype calls 0/1/2)</p>

<div class="cards">
  <div class="card"><h3>EPIC samples</h3><p>{n_epic}</p></div>
  <div class="card"><h3>H3A samples</h3><p>{n_h3a}</p></div>
  <div class="card"><h3>Common SNPs</h3><p>{len(common_snps)}</p></div>
  <div class="card"><h3>Total pairs</h3><p>{n_epic * n_h3a:,}</p></div>
  <div class="card"><h3>Reliable pairs</h3><p>{len(reliable_scores):,}</p></div>
</div>

{"<h2>Classification Summary</h2><div class='cards'>" + cat_cards + "</div>" if class_counts else ""}

<div class="note">
  <strong>Categories:</strong><br>
  <span style="color:#27ae60">&#9679;</span> <b>correct_assignment</b>
    – best match is the assigned H3A ID; concordance ≥ {THRESH_HIGH}<br>
  <span style="color:#e67e22">&#9679;</span> <b>swap_high_confidence</b>
    – best match is a <em>different</em> H3A ID; concordance ≥ {THRESH_HIGH}<br>
  <span style="color:#f39c12">&#9679;</span> <b>swap_moderate_confidence</b>
    – best match is different; concordance {THRESH_MOD}–{THRESH_HIGH}<br>
  <span style="color:#8e44ad">&#9679;</span> <b>outside_dataset</b>
    – best match is an H3A participant who was never enrolled for methylation
    (genotype-only individual)<br>
  <span style="color:#e74c3c">&#9679;</span> <b>low_confidence</b>
    – no reliable match found; concordance &lt; {THRESH_MOD} or DNA absent from H3A dataset<br>
</div>

{"<h2>Concordance Distribution</h2>" + hist_svg if hist_svg else ""}

<h2>Per-Sample Verdict</h2>
<table>
<thead><tr>
  <th>EPIC Sample</th><th>Assigned H3A</th><th>Assigned Score</th>
  <th>Best H3A Match</th><th>Best Score</th><th>SNPs</th><th>Classification</th>
</tr></thead>
<tbody>{verdict_rows}</tbody>
</table>

<p style="margin-top:28px;font-size:11px;color:#999">
  Generated by cross_study_concordance.py —
  min_snps={min_snps} | report_threshold={threshold} |
  high_conf≥{THRESH_HIGH} | mod_conf≥{THRESH_MOD}
</p>
</body>
</html>"""

    with open(path, "w") as fh:
        fh.write(html)
    log.info("HTML report written: %s", path)


# =============================================================================
# CLI
# =============================================================================

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Cross-study identity check: EPIC methylation vs H3Africa genotyping"
    )
    p.add_argument("--epic_traw",  required=True,
                   help="EPIC methylation .traw (from extract_snp_betas.r)")
    p.add_argument("--h3a_traw",   required=True,
                   help="H3A genotype .traw (relabelled by match_h3a_snps.r)")
    p.add_argument("--out",        required=True,  help="Output directory")
    p.add_argument("--prefix",     required=True,  help="Output file prefix")
    p.add_argument("--id_map",     default=None,
                   help="CSV with header 'h3a_id,epic_id'. "
                        "Maps each enrolled participant's H3A Sample_ID to their "
                        "EPIC Sample_ID.  Enables assigned-vs-best classification. "
                        "Without this, only best-match is reported.")
    p.add_argument("--threshold",  type=float, default=0.80,
                   help="Concordance threshold used in summary reporting (default: 0.80). "
                        "Does not affect classification thresholds "
                        f"(high≥{THRESH_HIGH}, mod≥{THRESH_MOD}).")
    p.add_argument("--min_snps",   type=int,   default=20,
                   help="Min jointly-called SNPs for a reliable pair (default: 20)")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    os.makedirs(args.out, exist_ok=True)

    log.info("Cross-study identity check")
    log.info("EPIC .traw  : %s", args.epic_traw)
    log.info("H3A .traw   : %s", args.h3a_traw)
    log.info("ID map      : %s", args.id_map or "(none — best-match only)")
    log.info("Threshold   : %.3f", args.threshold)
    log.info("Min SNPs    : %d",   args.min_snps)

    snps_e, samps_e, mat_e = parse_traw(args.epic_traw)
    snps_h, samps_h, mat_h = parse_traw(args.h3a_traw)

    if not samps_e or not samps_h:
        log.error("Need at least 1 sample in each dataset.")
        sys.exit(1)

    epic_to_h3a = None
    if args.id_map:
        epic_to_h3a, _ = parse_id_map(args.id_map)

    common_snps, mat_e_f, mat_h_f = intersect_snps(snps_e, mat_e, snps_h, mat_h)
    log.info("SNP intersection: %d EPIC / %d H3A → %d in common",
             len(snps_e), len(snps_h), len(common_snps))

    if len(common_snps) < args.min_snps:
        log.warning("Only %d common SNPs < min_snps=%d. "
                    "Check match_h3a_snps.r output.",
                    len(common_snps), args.min_snps)

    results = compute_concordance(mat_e_f, samps_e, mat_h_f, samps_h, args.min_snps)
    verdicts     = build_verdicts(results, samps_e, samps_h, epic_to_h3a)
    class_counts = write_classification_summary(verdicts,
                       os.path.join(args.out, f"{args.prefix}_classification_summary.tsv"))

    def out(suffix: str) -> str:
        return os.path.join(args.out, f"{args.prefix}_{suffix}")

    write_matrix(results, samps_e, samps_h,  out("concordance_matrix.csv"))
    write_pairs(results,                      out("concordance_pairs.tsv"))
    write_verdicts(verdicts,                  out("per_sample_verdict.tsv"))
    write_summary_txt(results, verdicts, class_counts,
                      common_snps=common_snps, n_epic=len(samps_e), n_h3a=len(samps_h),
                      threshold=args.threshold, min_snps=args.min_snps,
                      path=out("concordance_summary.txt"))
    write_html_report(results, verdicts, class_counts,
                      common_snps=common_snps, n_epic=len(samps_e), n_h3a=len(samps_h),
                      threshold=args.threshold, min_snps=args.min_snps,
                      path=out("concordance_report.html"))

    log.info("Done. %d EPIC samples processed.", len(samps_e))


if __name__ == "__main__":
    main()
