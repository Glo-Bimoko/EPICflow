#!/usr/bin/env python3
"""
generate_id_map.py

Produces ID_map.csv from the H3A plates manifest and the EPIC manifest.

Output columns:
  H3A_sample_id  — Sample_ID from all_6_plates.csv
  EPIC_sample_id — Sample_ID from all_6_epic.csv

Only samples present in BOTH files are written to the map.
Samples exclusive to one file are reported to stdout but not included.

Usage:
  python generate_id_map.py \
      --plates all_6_plates.csv \
      --epic   all_6_epic.csv \
      --out    ID_map.csv
"""

import argparse
import sys
import pandas as pd


def build_id_map(plates_path: str, epic_path: str, out_path: str) -> None:
    plates = pd.read_csv(plates_path, dtype=str)
    epic   = pd.read_csv(epic_path,   dtype=str)

    # Normalise column name (handle 'Sample ID' or 'Sample_ID')
    def get_sample_id_col(df: pd.DataFrame, source: str) -> str:
        for candidate in ("Sample_ID", "Sample ID"):
            if candidate in df.columns:
                return candidate
        raise KeyError(
            f"Could not find a 'Sample_ID' or 'Sample ID' column in {source}. "
            f"Available columns: {df.columns.tolist()}"
        )

    plates_col = get_sample_id_col(plates, plates_path)
    epic_col   = get_sample_id_col(epic,   epic_path)

    plates_ids = set(plates[plates_col].dropna().str.strip())
    epic_ids   = set(epic[epic_col].dropna().str.strip())

    both        = sorted(plates_ids & epic_ids)
    only_plates = sorted(plates_ids - epic_ids)
    only_epic   = sorted(epic_ids   - plates_ids)

    print(f"H3A (plates) samples : {len(plates_ids)}")
    print(f"EPIC samples         : {len(epic_ids)}")
    print(f"Matched (map rows)   : {len(both)}")
    if only_plates:
        print(f"Only in plates (excluded): {len(only_plates)} — e.g. {only_plates[:5]}")
    if only_epic:
        print(f"Only in EPIC   (excluded): {len(only_epic)} — e.g. {only_epic[:5]}")

    id_map = pd.DataFrame({
        "H3A_sample_id":  both,
        "EPIC_sample_id": both,
    })

    id_map.to_csv(out_path, index=False)
    print(f"\nWrote {len(id_map)} rows → {out_path}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--plates", default="all_6_plates.csv",
                        help="Path to H3A plates CSV (default: all_6_plates.csv)")
    parser.add_argument("--epic",   default="all_6_epic.csv",
                        help="Path to EPIC CSV (default: all_6_epic.csv)")
    parser.add_argument("--out",    default="ID_map.csv",
                        help="Output path (default: ID_map.csv)")
    args = parser.parse_args()

    try:
        build_id_map(args.plates, args.epic, args.out)
    except (FileNotFoundError, KeyError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
