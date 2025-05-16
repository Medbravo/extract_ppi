#!/usr/bin/env python3
"""collect_scores.py

Combine all Rosetta *.sc score files found under a given directory into a
single tab-separated table and sort the entries by ddG (ascending, most
negative first) and shape complementarity (descending, highest first).

Usage:
    python collect_scores.py /path/to/output_dir [--pattern "*.sc"]

The script writes `combined_scores.tsv` inside the provided directory.
"""
from __future__ import annotations

import argparse
import csv
import glob
import os
from pathlib import Path
from typing import List, Dict

import pandas as pd

# -----------------------------------------------------------------------------

def parse_scorefile(path: Path) -> List[Dict[str, str]]:
    """Parse a Rosetta score file (.sc) and return list of row dicts."""
    rows: List[Dict[str, str]] = []
    header = []
    with path.open() as fh:
        for line in fh:
            if not line.startswith("SCORE:"):
                continue
            parts = line.strip().split()
            # Remove the leading tag
            tag = parts.pop(0)  # SCORE:

            # This is the header line (contains alphabetic column names)
            if all(not p.replace('.', '', 1).lstrip('-').isdigit() for p in parts[:-1]):
                header = parts
                continue

            # Data line (values + description at end)
            if not header:
                # file is malformed (header after data)
                continue
            # Align header with parts; there may be extra description column at end
            row_dict = {h: v for h, v in zip(header, parts)}
            # Description (pose tag) is usually last entry, under column 'description'
            if len(parts) == len(header) + 1:
                row_dict["description"] = parts[-1]
            rows.append(row_dict)
    return rows


def main():
    parser = argparse.ArgumentParser(description="Consolidate Rosetta .sc files")
    parser.add_argument("directory", type=str, help="Root directory containing .sc files (searched recursively)")
    parser.add_argument("--pattern", default="*.sc", help="Glob pattern for .sc files (default: *.sc)")
    args = parser.parse_args()

    root = Path(args.directory).expanduser().resolve()
    if not root.is_dir():
        raise SystemExit(f"Provided path {root} is not a directory")

    files = [Path(p) for p in glob.glob(str(root / "**" / args.pattern), recursive=True)]
    if not files:
        raise SystemExit(f"No score files matching {args.pattern} found under {root}")

    all_rows: List[Dict[str, str]] = []
    for sc_file in files:
        rows = parse_scorefile(sc_file)
        for row in rows:
            row["scorefile"] = sc_file.name
            all_rows.append(row)

    if not all_rows:
        raise SystemExit("No valid SCORE lines parsed from score files.")

    df = pd.DataFrame(all_rows)

    # Attempt to cast numeric columns
    for col in df.columns:
        df[col] = pd.to_numeric(df[col], errors="ignore")

    # Determine columns representing ddg & shape complementarity
    ddg_cols = [c for c in df.columns if c.lower() in {"ddg", "interface_delta_x", "delta_g"}]
    sc_cols = [c for c in df.columns if c.lower().startswith("shape_complementarity") or c.lower() == "sc_value"]

    if not ddg_cols:
        print("WARNING: Could not find a ddG column (ddg/interface_delta_X). Sorting by total_score if present.")
    if not sc_cols:
        print("WARNING: Could not find shape complementarity column. Sorting by ddG only.")

    # Build sort order list
    sort_cols = []
    ascending = []
    if ddg_cols:
        sort_cols.append(ddg_cols[0])
        ascending.append(True)   # lower ddG is better
    if sc_cols:
        sort_cols.append(sc_cols[0])
        ascending.append(False)  # higher SC is better
    if not sort_cols:
        # fallback to total_score
        if "total_score" in df.columns:
            sort_cols.append("total_score")
            ascending.append(True)
        else:
            print("No suitable columns for sorting; leaving as-is.")

    df_sorted = df.sort_values(by=sort_cols, ascending=ascending)

    out_file = root / "combined_scores.tsv"
    df_sorted.to_csv(out_file, sep="\t", index=False, quoting=csv.QUOTE_NONE)
    print(f"Wrote {len(df_sorted)} rows to {out_file}")


if __name__ == "__main__":
    main() 