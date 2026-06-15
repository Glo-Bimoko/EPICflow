#!/usr/bin/env python3
"""
methylation_duplicate_check.py
===============================
For each pair flagged as a likely duplicate by pairwise_concordance.py,
extract their genome-wide beta values from the normalised GDS file and
compute methylation profile similarity.

Three complementary metrics are calculated per pair:
  - Pearson r      : sensitive to linear relationship; ~1.0 for true duplicates
  - Median |Δβ|    : absolute beta difference; should be near 0 for duplicates
  - CpGs |Δβ|>0.2  : fraction of probes with a large discordant swing; should
                     be near 0 for duplicates (biological noise stays <0.05)

A pair is called a CONFIRMED_DUPLICATE when ALL three criteria are met:
  Pearson r >= --r_threshold   AND
  Median |Δβ| <= --delta_threshold   AND
  Fraction |Δβ|>0.2 <= --large_delta_frac

Usage:
    python3 methylation_duplicate_check.py \\
        --flagged  <prefix>_flagged_duplicates.tsv \\
        --gds      beta_normalised.gds \\
        --out      <out_dir> \\
        --prefix   <run_prefix> \\
        [--r_threshold         0.999] \\
        [--delta_threshold     0.01 ] \\
        [--large_delta_frac    0.01 ] \\
        [--chunk_size          5000 ]

Outputs:
    <prefix>_methylation_duplicate_evidence.tsv   – per-pair metrics + verdict
    <prefix>_methylation_duplicate_summary.txt    – plain-text summary
    <prefix>_methylation_duplicate_report.html    – interactive HTML report

Dependencies:
    pip install numpy pandas scipy h5py

Notes:
    h5py is used to read the GDS file directly (GDS files are HDF5 under the
    hood), which avoids the rpy2/R dependency entirely.  The beta matrix is
    stored in GDS as a 2-D float array under the node "genotype" (meffil
    convention: sites × samples).  The file is read in chunks (--chunk_size
    CpGs) to stay memory-safe for large datasets (EPICv2 ~900k probes).
"""

import argparse
import csv
import logging
import os
import sys
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from scipy.stats import pearsonr

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# GDS reader — reads meffil beta_normalised.gds directly via h5py
#
# meffil writes a GDS file that is valid HDF5.  The relevant nodes are:
#   /sample.id   — 1-D array of sample name strings  (length = n_samples)
#   /snp.id      — 1-D array of CpG site name strings (length = n_sites)
#   /genotype    — 2-D float array, shape (n_sites, n_samples), beta values
#
# h5py reads these natively with no R or rpy2 dependency.
# ---------------------------------------------------------------------------

def _open_gds(gds_path: str):
    """Open the GDS file as HDF5 and return the h5py File handle."""
    try:
        import h5py
    except ImportError:
        log.error("h5py is required for GDS access.  Install with:  conda install h5py")
        sys.exit(1)
    return h5py.File(gds_path, "r")


def _decode(arr) -> List[str]:
    """Convert an h5py dataset to a plain Python list of strings."""
    return [v.decode() if isinstance(v, bytes) else str(v) for v in arr]


def list_gds_samples(gds_path: str) -> List[str]:
    """Return sample names stored in the GDS file."""
    with _open_gds(gds_path) as f:
        return _decode(f["sample.id"][:])


def list_gds_sites(gds_path: str) -> List[str]:
    """Return all CpG site names stored in the GDS file."""
    with _open_gds(gds_path) as f:
        return _decode(f["snp.id"][:])


def fetch_betas_for_samples(
    gds_path: str,
    sites: List[str],
    samples: List[str],
) -> pd.DataFrame:
    """
    Load a beta sub-matrix for the given sites and samples.

    Returns a DataFrame of shape (n_sites, n_samples), index = site names,
    columns = sample names.  Values are beta [0, 1]; NaN = missing.

    The GDS /genotype node is shape (n_all_sites, n_all_samples).  We select
    only the rows (sites) and columns (samples) we need using integer indices
    so that h5py performs the slice in a single read — no full matrix in RAM.
    """
    with _open_gds(gds_path) as f:
        all_samples = _decode(f["sample.id"][:])
        all_sites   = _decode(f["snp.id"][:])

        # Build index maps
        sample_idx = {s: i for i, s in enumerate(all_samples)}
        site_idx   = {s: i for i, s in enumerate(all_sites)}

        row_indices = np.array([site_idx[s]   for s in sites   if s in site_idx],   dtype=np.intp)
        col_indices = np.array([sample_idx[s] for s in samples if s in sample_idx], dtype=np.intp)

        valid_sites   = [s for s in sites   if s in site_idx]
        valid_samples = [s for s in samples if s in sample_idx]

        if len(row_indices) == 0 or len(col_indices) == 0:
            return pd.DataFrame(index=valid_sites, columns=valid_samples, dtype=float)

        # h5py fancy-index: rows first, then columns
        geno = f["genotype"]
        # Fancy indexing on both axes simultaneously is not supported in h5py;
        # slice rows first (cheaper — sites are the outer dimension), then cols.
        mat = geno[np.sort(row_indices), :][:, np.sort(col_indices)]

        # Restore requested order (row_indices / col_indices may not be sorted)
        row_order = np.argsort(np.argsort(row_indices))
        col_order = np.argsort(np.argsort(col_indices))
        mat = mat[row_order, :][:, col_order]

        # Replace sentinel missing values (meffil uses NaN directly for floats)
        mat = mat.astype(float)

    return pd.DataFrame(mat, index=valid_sites, columns=valid_samples)


# ---------------------------------------------------------------------------
# Per-pair methylation metrics
# ---------------------------------------------------------------------------

def methylation_metrics(
    beta_a: np.ndarray,
    beta_b: np.ndarray,
    large_delta_cutoff: float = 0.2,
) -> Dict:
    """
    Given two 1-D arrays of beta values (NaN = missing), compute:
      - n_cpgs        : number of CpGs where BOTH samples have a value
      - pearson_r     : Pearson correlation
      - pearson_p     : two-sided p-value
      - median_abs_delta : median |beta_a − beta_b|
      - mean_abs_delta
      - frac_large_delta : fraction of jointly-called CpGs with |Δβ| > cutoff
    """
    valid = ~(np.isnan(beta_a) | np.isnan(beta_b))
    n_cpgs = int(valid.sum())

    if n_cpgs < 2:
        return {
            "n_cpgs": n_cpgs,
            "pearson_r": float("nan"),
            "pearson_p": float("nan"),
            "median_abs_delta": float("nan"),
            "mean_abs_delta": float("nan"),
            "frac_large_delta": float("nan"),
        }

    a = beta_a[valid]
    b = beta_b[valid]
    delta = np.abs(a - b)

    r, p = pearsonr(a, b)

    return {
        "n_cpgs": n_cpgs,
        "pearson_r": float(r),
        "pearson_p": float(p),
        "median_abs_delta": float(np.median(delta)),
        "mean_abs_delta": float(np.mean(delta)),
        "frac_large_delta": float((delta > large_delta_cutoff).sum() / n_cpgs),
    }


def classify_pair(
    metrics: Dict,
    r_threshold: float,
    delta_threshold: float,
    large_delta_frac: float,
) -> str:
    """Return a verdict string for a flagged pair based on methylation metrics."""
    r   = metrics["pearson_r"]
    med = metrics["median_abs_delta"]
    fld = metrics["frac_large_delta"]

    if any(np.isnan(v) for v in [r, med, fld]):
        return "INSUFFICIENT_DATA"

    if r >= r_threshold and med <= delta_threshold and fld <= large_delta_frac:
        return "CONFIRMED_DUPLICATE"
    elif r >= r_threshold - 0.01 or med <= delta_threshold * 2:
        # High genotype concordance but borderline methylation — worth flagging
        return "BORDERLINE"
    else:
        return "METHYLATION_MISMATCH"


# ---------------------------------------------------------------------------
# Main analysis
# ---------------------------------------------------------------------------

def run_analysis(
    flagged_path: str,
    gds_path: str,
    out_dir: str,
    prefix: str,
    r_threshold: float,
    delta_threshold: float,
    large_delta_frac: float,
    chunk_size: int,
) -> pd.DataFrame:
    """
    Core logic: for each flagged pair, pull beta values from GDS and compute
    methylation similarity metrics.

    Betas are loaded in chunks of chunk_size sites to avoid memory spikes.
    Only the unique set of samples involved in flagged pairs is loaded.
    """
    # ── Read flagged pairs ────────────────────────────────────────────────
    flagged = pd.read_csv(flagged_path, sep="\t")
    if flagged.empty:
        log.info("No flagged duplicate pairs — nothing to do.")
        return pd.DataFrame()

    log.info("Flagged pairs to assess: %d", len(flagged))

    # Unique samples needed
    samples_needed = sorted(
        set(flagged["sample_a"].tolist() + flagged["sample_b"].tolist())
    )
    log.info("Unique samples to load betas for: %d", len(samples_needed))

    # ── Get all CpG site names from GDS ──────────────────────────────────
    log.info("Reading site list from GDS: %s", gds_path)
    all_sites = list_gds_sites(gds_path)
    n_sites = len(all_sites)
    log.info("Total CpG sites in GDS: %d", n_sites)

    # ── Accumulate metrics per pair in chunks ─────────────────────────────
    # We initialise accumulators as lists and reduce at the end to avoid
    # holding the full (n_sites × n_samples) matrix in memory at once.

    # Accumulators: per pair, lists of per-chunk arrays
    pair_keys = list(zip(flagged["sample_a"], flagged["sample_b"]))
    # {(a,b): {"deltas": [], "a_vals": [], "b_vals": []}}
    pair_acc: Dict[Tuple, Dict] = {
        k: {"deltas": [], "a_vals": [], "b_vals": []}
        for k in pair_keys
    }

    n_chunks = (n_sites + chunk_size - 1) // chunk_size
    log.info(
        "Loading betas in %d chunks of ≤%d sites …", n_chunks, chunk_size
    )

    for chunk_idx in range(n_chunks):
        start = chunk_idx * chunk_size
        end   = min(start + chunk_size, n_sites)
        chunk_sites = all_sites[start:end]

        if (chunk_idx + 1) % 10 == 0 or chunk_idx == n_chunks - 1:
            log.info(
                "  Chunk %d/%d  (sites %d–%d)",
                chunk_idx + 1, n_chunks, start + 1, end,
            )

        beta_chunk = fetch_betas_for_samples(gds_path, chunk_sites, samples_needed)

        for sample_a, sample_b in pair_keys:
            if sample_a not in beta_chunk.columns or sample_b not in beta_chunk.columns:
                continue
            a_vals = beta_chunk[sample_a].values
            b_vals = beta_chunk[sample_b].values
            pair_acc[(sample_a, sample_b)]["a_vals"].append(a_vals)
            pair_acc[(sample_a, sample_b)]["b_vals"].append(b_vals)

    # ── Compute final metrics per pair ────────────────────────────────────
    log.info("Computing methylation metrics for %d pairs …", len(pair_keys))
    results = []

    for _, row in flagged.iterrows():
        sample_a = row["sample_a"]
        sample_b = row["sample_b"]
        key = (sample_a, sample_b)

        acc = pair_acc.get(key, {})
        a_all = np.concatenate(acc.get("a_vals", [])) if acc.get("a_vals") else np.array([])
        b_all = np.concatenate(acc.get("b_vals", [])) if acc.get("b_vals") else np.array([])

        if len(a_all) == 0 or len(b_all) == 0:
            metrics = {
                "n_cpgs": 0,
                "pearson_r": float("nan"),
                "pearson_p": float("nan"),
                "median_abs_delta": float("nan"),
                "mean_abs_delta": float("nan"),
                "frac_large_delta": float("nan"),
            }
        else:
            metrics = methylation_metrics(
                a_all, b_all, large_delta_cutoff=0.2
            )

        verdict = classify_pair(
            metrics, r_threshold, delta_threshold, large_delta_frac
        )

        results.append(
            {
                "sample_a":             sample_a,
                "sample_b":             sample_b,
                "genotype_concordance": row["concordance"],
                "genotype_n_snps":      row["n_compared"],
                **metrics,
                "verdict":              verdict,
            }
        )

    return pd.DataFrame(results)


# ---------------------------------------------------------------------------
# Output writers
# ---------------------------------------------------------------------------

def write_evidence_table(df: pd.DataFrame, path: str) -> None:
    float_cols = [
        "genotype_concordance", "pearson_r", "pearson_p",
        "median_abs_delta", "mean_abs_delta", "frac_large_delta",
    ]
    out = df.copy()
    for col in float_cols:
        if col in out.columns:
            out[col] = out[col].apply(
                lambda v: f"{v:.6f}" if pd.notna(v) else "NA"
            )
    out.to_csv(path, sep="\t", index=False)
    log.info("Evidence table written: %s", path)


def write_summary(df: pd.DataFrame, thresholds: Dict, path: str) -> None:
    lines = [
        "=" * 65,
        "METHYLATION CONCORDANCE SUMMARY  (genotype-flagged pairs)",
        "=" * 65,
        f"  Total pairs assessed           : {len(df)}",
    ]

    if len(df) > 0:
        for verdict in ["CONFIRMED_DUPLICATE", "BORDERLINE",
                        "METHYLATION_MISMATCH", "INSUFFICIENT_DATA"]:
            n = (df["verdict"] == verdict).sum()
            lines.append(f"  {verdict:<30}: {n}")

        confirmed = df[df["verdict"] == "CONFIRMED_DUPLICATE"]
        mismatch  = df[df["verdict"] == "METHYLATION_MISMATCH"]

        lines += [
            "",
            "  Methylation thresholds used:",
            f"    Pearson r     ≥ {thresholds['r_threshold']:.4f}",
            f"    Median |Δβ|   ≤ {thresholds['delta_threshold']:.4f}",
            f"    Frac |Δβ|>0.2 ≤ {thresholds['large_delta_frac']:.4f}",
        ]

        if len(confirmed) > 0:
            lines += [
                "",
                "  CONFIRMED DUPLICATES (same sample, same methylation):",
            ]
            for _, r in confirmed.iterrows():
                lines.append(
                    f"    {r['sample_a']}  ↔  {r['sample_b']}"
                    f"  r={r['pearson_r']:.4f}  median|Δβ|={r['median_abs_delta']:.4f}"
                )

        if len(mismatch) > 0:
            lines += [
                "",
                "  METHYLATION MISMATCHES (same genotype, different methylation):",
                "  → Possible sample swap or monozygotic twins",
            ]
            for _, r in mismatch.iterrows():
                lines.append(
                    f"    {r['sample_a']}  ↔  {r['sample_b']}"
                    f"  r={r['pearson_r']:.4f}  median|Δβ|={r['median_abs_delta']:.4f}"
                )

    lines.append("=" * 65)
    text = "\n".join(lines)
    with open(path, "w") as fh:
        fh.write(text + "\n")
    print("\n" + text)
    log.info("Summary written: %s", path)


def write_html_report(df: pd.DataFrame, thresholds: Dict, path: str) -> None:
    """Generate a self-contained HTML report with an interactive pair table."""

    def fmt(v, decimals=4):
        if pd.isna(v):
            return "—"
        try:
            return f"{float(v):.{decimals}f}"
        except (TypeError, ValueError):
            return str(v)

    verdict_colour = {
        "CONFIRMED_DUPLICATE":   "#2ecc71",   # green
        "BORDERLINE":            "#f39c12",   # amber
        "METHYLATION_MISMATCH":  "#e74c3c",   # red
        "INSUFFICIENT_DATA":     "#95a5a6",   # grey
    }
    verdict_label = {
        "CONFIRMED_DUPLICATE":   "✔ Confirmed duplicate",
        "BORDERLINE":            "⚠ Borderline",
        "METHYLATION_MISMATCH":  "✘ Methylation mismatch",
        "INSUFFICIENT_DATA":     "? Insufficient data",
    }

    rows_html = ""
    for _, r in df.iterrows():
        v = r.get("verdict", "INSUFFICIENT_DATA")
        colour = verdict_colour.get(v, "#95a5a6")
        label  = verdict_label.get(v, v)
        rows_html += f"""
        <tr>
          <td>{r['sample_a']}</td>
          <td>{r['sample_b']}</td>
          <td>{fmt(r.get('genotype_concordance'), 4)}</td>
          <td>{r.get('genotype_n_snps', '—')}</td>
          <td>{r.get('n_cpgs', '—')}</td>
          <td>{fmt(r.get('pearson_r'), 6)}</td>
          <td>{fmt(r.get('median_abs_delta'), 6)}</td>
          <td>{fmt(r.get('mean_abs_delta'), 6)}</td>
          <td>{fmt(r.get('frac_large_delta'), 4)}</td>
          <td style="color:{colour};font-weight:bold">{label}</td>
        </tr>"""

    n_total     = len(df)
    n_confirmed = (df["verdict"] == "CONFIRMED_DUPLICATE").sum() if n_total else 0
    n_mismatch  = (df["verdict"] == "METHYLATION_MISMATCH").sum() if n_total else 0
    n_border    = (df["verdict"] == "BORDERLINE").sum() if n_total else 0

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>Methylation Duplicate Evidence Report</title>
<style>
  body      {{ font-family: Arial, sans-serif; margin: 30px; color: #333; }}
  h1        {{ color: #2c3e50; }}
  h2        {{ color: #34495e; border-bottom: 1px solid #ccc; padding-bottom: 4px; }}
  .cards    {{ display: flex; gap: 20px; margin: 20px 0; flex-wrap: wrap; }}
  .card     {{ background: #f8f9fa; border: 1px solid #dee2e6; border-radius: 8px;
               padding: 16px 24px; text-align: center; min-width: 140px; }}
  .card .n  {{ font-size: 2.2em; font-weight: bold; }}
  .card .l  {{ font-size: 0.9em; color: #666; margin-top: 4px; }}
  .green    {{ color: #2ecc71; }}
  .amber    {{ color: #f39c12; }}
  .red      {{ color: #e74c3c; }}
  table     {{ border-collapse: collapse; width: 100%; margin-top: 16px;
               font-size: 0.88em; }}
  th        {{ background: #2c3e50; color: white; padding: 8px 10px;
               text-align: left; }}
  td        {{ padding: 7px 10px; border-bottom: 1px solid #eee; }}
  tr:hover  {{ background: #f5f5f5; }}
  .thresh   {{ background: #f0f4f8; border-left: 4px solid #2c3e50;
               padding: 12px 18px; margin: 12px 0; font-size: 0.92em; }}
  .thresh b {{ display: inline-block; width: 220px; }}
  .note     {{ font-size: 0.85em; color: #555; margin-top: 6px; }}
</style>
</head>
<body>
<h1>Methylation Duplicate Evidence Report</h1>

<p>
Pairs below were flagged by genotype concordance (≥{thresholds['r_threshold']} on rs-probe SNPs).
This report adds genome-wide methylation profile comparison to confirm or refute each flag.
</p>

<div class="cards">
  <div class="card">
    <div class="n">{n_total}</div>
    <div class="l">Pairs assessed</div>
  </div>
  <div class="card">
    <div class="n green">{n_confirmed}</div>
    <div class="l">Confirmed duplicates</div>
  </div>
  <div class="card">
    <div class="n amber">{n_border}</div>
    <div class="l">Borderline</div>
  </div>
  <div class="card">
    <div class="n red">{n_mismatch}</div>
    <div class="l">Methylation mismatches<br><span class="note">(possible swap / MZ twins)</span></div>
  </div>
</div>

<h2>Thresholds applied</h2>
<div class="thresh">
  <b>Pearson r ≥</b> {thresholds['r_threshold']}<br>
  <b>Median |Δβ| ≤</b> {thresholds['delta_threshold']}<br>
  <b>Fraction |Δβ| &gt; 0.2 ≤</b> {thresholds['large_delta_frac']}<br>
</div>
<p class="note">
  A pair is CONFIRMED_DUPLICATE only when <em>all three</em> methylation criteria are met
  in addition to the genotype concordance flag.  METHYLATION_MISMATCH means identical
  genotypes but divergent methylation — consistent with a sample swap or monozygotic twins.
</p>

<h2>Pair-level evidence</h2>
<table>
  <thead>
    <tr>
      <th>Sample A</th><th>Sample B</th>
      <th>Geno concordance</th><th>Geno SNPs</th>
      <th>CpGs compared</th>
      <th>Pearson r</th>
      <th>Median |Δβ|</th>
      <th>Mean |Δβ|</th>
      <th>Frac |Δβ|>0.2</th>
      <th>Verdict</th>
    </tr>
  </thead>
  <tbody>
    {rows_html if rows_html else '<tr><td colspan="10" style="text-align:center;color:#888">No pairs to display</td></tr>'}
  </tbody>
</table>

</body>
</html>"""

    with open(path, "w") as fh:
        fh.write(html)
    log.info("HTML report written: %s", path)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Methylation profile concordance for genotype-flagged duplicate pairs"
    )
    p.add_argument("--flagged",   required=True,
                   help="flagged_duplicates.tsv from pairwise_concordance.py")
    p.add_argument("--gds",       required=True,
                   help="Normalised beta GDS file (beta_normalised.gds)")
    p.add_argument("--out",       required=True, help="Output directory")
    p.add_argument("--prefix",    required=True, help="Output file prefix")
    p.add_argument("--r_threshold",      type=float, default=0.999,
                   help="Pearson r threshold for CONFIRMED_DUPLICATE (default: 0.999)")
    p.add_argument("--delta_threshold",  type=float, default=0.01,
                   help="Max median |Δβ| for CONFIRMED_DUPLICATE (default: 0.01)")
    p.add_argument("--large_delta_frac", type=float, default=0.01,
                   help="Max fraction of |Δβ|>0.2 for CONFIRMED_DUPLICATE (default: 0.01)")
    p.add_argument("--chunk_size", type=int, default=5000,
                   help="CpGs per GDS read chunk (default: 5000; lower = less memory)")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    os.makedirs(args.out, exist_ok=True)

    # Check inputs
    for path, label in [(args.flagged, "--flagged"), (args.gds, "--gds")]:
        if not os.path.exists(path):
            log.error("%s path not found: %s", label, path)
            sys.exit(1)

    thresholds = {
        "r_threshold":      args.r_threshold,
        "delta_threshold":  args.delta_threshold,
        "large_delta_frac": args.large_delta_frac,
    }

    def out(suffix: str) -> str:
        return os.path.join(args.out, f"{args.prefix}_{suffix}")

    df = run_analysis(
        flagged_path     = args.flagged,
        gds_path         = args.gds,
        out_dir          = args.out,
        prefix           = args.prefix,
        r_threshold      = args.r_threshold,
        delta_threshold  = args.delta_threshold,
        large_delta_frac = args.large_delta_frac,
        chunk_size       = args.chunk_size,
    )

    if df.empty:
        log.info("No pairs to report.")
        # Write empty placeholder files so the process always has outputs
        for suffix in [
            "methylation_duplicate_evidence.tsv",
            "methylation_duplicate_summary.txt",
            "methylation_duplicate_report.html",
        ]:
            with open(out(suffix), "w") as fh:
                fh.write("# No flagged pairs\n")
        return

    write_evidence_table(df,             out("methylation_duplicate_evidence.tsv"))
    write_summary(df, thresholds,        out("methylation_duplicate_summary.txt"))
    write_html_report(df, thresholds,    out("methylation_duplicate_report.html"))


if __name__ == "__main__":
    main()