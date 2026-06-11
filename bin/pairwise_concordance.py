#!/usr/bin/env python3
"""
pairwise_concordance.py
=======================
All-vs-all pairwise genotype concordance from a PLINK-style .traw file.

Adapted for the meffil methylation pipeline: input is genotype calls (0/1/2)
derived from the 65 rs-prefixed SNP probes on the EPIC/EPICv2 array.

Logic matches H3Aflow's concordance module:
  - For each sample pair, count positions where BOTH samples have a valid call.
  - Concordance = (matching calls) / (jointly-called positions).
  - Pairs with fewer than MIN_SNPS jointly called positions are flagged as
    unreliable and reported but excluded from the duplicate-flag threshold.
  - Pairs at or above CONCORDANCE_THRESHOLD are reported as likely duplicates
    or sample swaps.

Usage:
    python3 pairwise_concordance.py \\
        --traw   <prefix>_snp_calls.traw \\
        --out    <out_dir> \\
        --prefix <run_prefix> \\
        [--threshold 0.99] \\
        [--min_snps  50]

Outputs:
    <prefix>_concordance_matrix.csv   – N×N matrix of concordance values
    <prefix>_concordance_pairs.tsv    – long-form pairwise results
    <prefix>_flagged_duplicates.tsv   – pairs above threshold
    <prefix>_concordance_summary.txt  – plain-text run summary
"""

import argparse
import csv
import itertools
import logging
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


# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------

def parse_traw(traw_path: str) -> Tuple[List[str], List[str], np.ndarray]:
    """
    Parse a PLINK transposed genotype file (.traw).

    Columns: CHR  SNP  CM  POS  COUNTED  ALT  sample1  sample2 ...
    Genotype encoding: 0=AA, 1=AB, 2=BB, NA=missing

    Returns
    -------
    snp_ids    : list of SNP identifiers (one per row)
    sample_ids : list of sample identifiers
    gt_matrix  : int8 array, shape (n_snps, n_samples); -1 = missing
    """
    log.info("Parsing .traw file: %s", traw_path)

    snp_ids: List[str] = []
    sample_ids: List[str] = []
    rows: List[List[int]] = []

    with open(traw_path, "r", newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader)
        # Columns 0-5: CHR SNP CM POS COUNTED ALT  → data starts at col 6
        sample_ids = header[6:]

        for line in reader:
            if not line:
                continue
            snp_ids.append(line[1])
            calls = []
            for val in line[6:]:
                val = val.strip()
                if val in ("NA", "nan", "", "."):
                    calls.append(-1)
                else:
                    try:
                        calls.append(int(float(val)))
                    except ValueError:
                        calls.append(-1)
            rows.append(calls)

    gt_matrix = np.array(rows, dtype=np.int8)
    log.info(
        "Loaded %d SNPs × %d samples", gt_matrix.shape[0], gt_matrix.shape[1]
    )
    return snp_ids, sample_ids, gt_matrix


# ---------------------------------------------------------------------------
# Concordance calculation
# ---------------------------------------------------------------------------

def pairwise_concordance(
    gt_matrix: np.ndarray,
    sample_ids: List[str],
    min_snps: int = 50,
) -> List[Dict]:
    """
    Compute all-vs-all concordance.

    Returns a list of dicts with keys:
        sample_a, sample_b, n_compared, n_concordant, concordance, reliable
    """
    n_samples = gt_matrix.shape[1]
    results = []

    total_pairs = n_samples * (n_samples - 1) // 2
    log.info("Computing concordance for %d sample pairs …", total_pairs)

    for i, j in itertools.combinations(range(n_samples), 2):
        col_i = gt_matrix[:, i]
        col_j = gt_matrix[:, j]

        # Positions where both samples have a valid call
        valid_mask = (col_i >= 0) & (col_j >= 0)
        n_compared = int(valid_mask.sum())

        if n_compared == 0:
            concordance = float("nan")
            n_concordant = 0
        else:
            n_concordant = int((col_i[valid_mask] == col_j[valid_mask]).sum())
            concordance = n_concordant / n_compared

        results.append(
            {
                "sample_a": sample_ids[i],
                "sample_b": sample_ids[j],
                "n_compared": n_compared,
                "n_concordant": n_concordant,
                "concordance": concordance,
                "reliable": n_compared >= min_snps,
            }
        )

    return results


def build_matrix(
    results: List[Dict], sample_ids: List[str]
) -> np.ndarray:
    """Build symmetric N×N concordance matrix from pair list."""
    n = len(sample_ids)
    idx = {s: i for i, s in enumerate(sample_ids)}
    mat = np.full((n, n), np.nan)
    np.fill_diagonal(mat, 1.0)

    for r in results:
        i = idx[r["sample_a"]]
        j = idx[r["sample_b"]]
        val = r["concordance"]
        mat[i, j] = val
        mat[j, i] = val

    return mat


# ---------------------------------------------------------------------------
# Output writers
# ---------------------------------------------------------------------------

def write_matrix(mat: np.ndarray, sample_ids: List[str], path: str) -> None:
    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow([""] + sample_ids)
        for i, sid in enumerate(sample_ids):
            row_vals = [
                "" if np.isnan(mat[i, j]) else f"{mat[i, j]:.6f}"
                for j in range(len(sample_ids))
            ]
            writer.writerow([sid] + row_vals)
    log.info("Concordance matrix written: %s", path)


def write_pairs(results: List[Dict], path: str) -> None:
    fieldnames = [
        "sample_a", "sample_b",
        "n_compared", "n_concordant",
        "concordance", "reliable",
    ]
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for r in results:
            row = dict(r)
            if isinstance(row["concordance"], float) and not np.isnan(row["concordance"]):
                row["concordance"] = f"{row['concordance']:.6f}"
            writer.writerow(row)
    log.info("Pairwise results written: %s  (%d pairs)", path, len(results))


def write_flagged(
    results: List[Dict], threshold: float, path: str
) -> int:
    flagged = [
        r for r in results
        if r["reliable"]
        and isinstance(r["concordance"], float)
        and not np.isnan(r["concordance"])
        and r["concordance"] >= threshold
    ]
    flagged.sort(key=lambda x: x["concordance"], reverse=True)

    fieldnames = [
        "sample_a", "sample_b",
        "n_compared", "n_concordant",
        "concordance", "flag",
    ]
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for r in flagged:
            writer.writerow(
                {
                    "sample_a": r["sample_a"],
                    "sample_b": r["sample_b"],
                    "n_compared": r["n_compared"],
                    "n_concordant": r["n_concordant"],
                    "concordance": f"{r['concordance']:.6f}",
                    "flag": "LIKELY_DUPLICATE",
                }
            )

    if flagged:
        log.warning(
            "⚠  %d pair(s) flagged as likely duplicates (concordance ≥ %.3f):",
            len(flagged), threshold,
        )
        for r in flagged:
            log.warning(
                "   %s  ↔  %s  (concordance=%.4f, n_snps=%d)",
                r["sample_a"], r["sample_b"],
                r["concordance"], r["n_compared"],
            )
    else:
        log.info("No duplicate pairs detected above threshold %.3f", threshold)

    log.info("Flagged duplicates written: %s", path)
    return len(flagged)


def write_summary(
    results: List[Dict],
    n_samples: int,
    n_snps: int,
    n_flagged: int,
    threshold: float,
    min_snps: int,
    path: str,
) -> None:
    concordances = [
        r["concordance"]
        for r in results
        if r["reliable"]
        and isinstance(r["concordance"], float)
        and not np.isnan(r["concordance"])
    ]

    unreliable = sum(1 for r in results if not r["reliable"])

    lines = [
        "=" * 60,
        "GENOTYPE CONCORDANCE SUMMARY",
        "=" * 60,
        f"  Samples analysed      : {n_samples}",
        f"  SNP probes used       : {n_snps}",
        f"  Total pairs           : {len(results)}",
        f"  Reliable pairs        : {len(results) - unreliable}",
        f"    (≥ {min_snps} jointly called SNPs)",
        f"  Unreliable pairs      : {unreliable}",
        "",
        "  Concordance distribution (reliable pairs):",
    ]

    if concordances:
        arr = np.array(concordances)
        lines += [
            f"    Min    : {arr.min():.4f}",
            f"    Median : {np.median(arr):.4f}",
            f"    Mean   : {arr.mean():.4f}",
            f"    Max    : {arr.max():.4f}",
            f"    Std    : {arr.std():.4f}",
        ]
        # Distribution buckets
        buckets = [
            (0.00, 0.70, "unrelated   (0.00–0.70)"),
            (0.70, 0.85, "1st-degree  (0.70–0.85)"),
            (0.85, 0.99, "ambiguous   (0.85–0.99)"),
            (0.99, 1.01, "duplicates  (≥ 0.99)   "),
        ]
        lines.append("")
        lines.append("  Concordance buckets:")
        for lo, hi, label in buckets:
            count = int(((arr >= lo) & (arr < hi)).sum())
            lines.append(f"    {label}: {count}")
    else:
        lines.append("    (no reliable pairs found)")

    lines += [
        "",
        f"  Duplicate threshold   : ≥ {threshold:.3f}",
        f"  Flagged duplicates    : {n_flagged}",
        "=" * 60,
    ]

    text = "\n".join(lines)
    with open(path, "w") as fh:
        fh.write(text + "\n")

    # Also echo to stdout
    print("\n" + text)
    log.info("Summary written: %s", path)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="All-vs-all pairwise concordance from EPIC array rs-probe calls"
    )
    p.add_argument("--traw",      required=True, help=".traw genotype file")
    p.add_argument("--out",       required=True, help="Output directory")
    p.add_argument("--prefix",    required=True, help="Output file prefix")
    p.add_argument(
        "--threshold", type=float, default=0.99,
        help="Concordance threshold for flagging duplicates (default: 0.99)",
    )
    p.add_argument(
        "--min_snps", type=int, default=50,
        help="Minimum jointly-called SNPs required for a reliable comparison "
             "(default: 50; use lower values if many probes are missing)",
    )
    return p.parse_args()


def main() -> None:
    args = parse_args()
    os.makedirs(args.out, exist_ok=True)

    # Parse input
    snp_ids, sample_ids, gt_matrix = parse_traw(args.traw)
    n_snps   = gt_matrix.shape[0]
    n_samples = gt_matrix.shape[1]

    if n_samples < 2:
        log.error("Need at least 2 samples for concordance analysis.")
        sys.exit(1)

    log.info("Running concordance: %d samples, %d SNP probes", n_samples, n_snps)
    log.info("Duplicate threshold : %.3f", args.threshold)
    log.info("Min SNPs (reliable) : %d",   args.min_snps)

    # Compute
    results = pairwise_concordance(gt_matrix, sample_ids, min_snps=args.min_snps)

    # Build output paths
    def out(suffix: str) -> str:
        return os.path.join(args.out, f"{args.prefix}_{suffix}")

    # Write outputs
    mat = build_matrix(results, sample_ids)
    write_matrix(mat,     sample_ids, out("concordance_matrix.csv"))
    write_pairs(results,              out("concordance_pairs.tsv"))
    n_flagged = write_flagged(results, args.threshold, out("flagged_duplicates.tsv"))
    write_summary(
        results,
        n_samples=n_samples,
        n_snps=n_snps,
        n_flagged=n_flagged,
        threshold=args.threshold,
        min_snps=args.min_snps,
        path=out("concordance_summary.txt"),
    )


if __name__ == "__main__":
    main()
