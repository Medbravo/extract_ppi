#!/usr/bin/env python3
"""copy_success_pdbs.py

Copy all PDB files generated for a given design job that passed the SUCCESS
filter (as indicated in the *_confidence_summary.tsv produced by
`evaluate_predictions.py`).

The script expects the evaluation pipeline directory layout used in the
README, i.e.

- input Alphafold predictions are located under `in/<job_id>/`.
- The TSV summary is located in the *current* directory (evaluation_pipeline)
  and is named `<job_id>_confidence_summary.tsv`.
- Destination directory will be created alongside the TSV and will be named
  `<job_id>_success_pdbs`.

Usage:
    python copy_success_pdbs.py <job_id>

Example:
    python copy_success_pdbs.py 1747410238309-7728-b4dab45c_job
"""
from __future__ import annotations

import argparse
import shutil
from pathlib import Path

import pandas as pd


def copy_pdbs(job_id: str) -> None:
    root = Path(__file__).resolve().parent  # evaluation_pipeline directory
    tsv_path = root / f"{job_id}_confidence_summary.tsv"
    if not tsv_path.exists():
        raise SystemExit(f"Could not find TSV summary: {tsv_path}")

    # Read the TSV and filter successful entries
    df = pd.read_csv(tsv_path, sep="\t")
    if "SUCCESS" not in df.columns:
        raise SystemExit("The TSV file does not contain a 'SUCCESS' column.")

    success_df = df[df["SUCCESS"].astype(str).str.upper() == "YES"]
    if success_df.empty:
        print("No successful designs found; nothing to copy.")
        return

    # Input & destination paths: the PDBs are stored under the project-level
    #   ../in/<job_id>/ directory (sibling of evaluation_pipeline).
    job_input_dir = root.parent / "in" / job_id
    if not job_input_dir.is_dir():
        raise SystemExit(f"Input directory not found: {job_input_dir}")

    dest_dir = root / f"{job_id}_success_pdbs"
    dest_dir.mkdir(exist_ok=True)

    num_copied = 0
    for rel_dir in success_df["directory"]:
        # The final component of the directory string corresponds to the task
        # ID (e.g. 1747410239544-9037-05d4947d).  Alphafold writes a PDB named
        # exactly "<task_id>_model.pdb" in the job directory root.  Build that
        # specific filename directly.
        task_id = Path(rel_dir).name
        pdb_path = job_input_dir / f"{task_id}_model.pdb"
        matches = [pdb_path] if pdb_path.exists() else []
        if not matches:
            print(f"WARNING: no PDB matching '{pdb_path}' found in {job_input_dir}")
            continue

        for pdb_path in matches:
            target_path = dest_dir / pdb_path.name
            # Avoid unintentional overwrite (multiple tasks could share the same
            # filename if there are relaxed+unrelaxed variants).
            if target_path.exists():
                # Append a counter to preserve all files
                stem = target_path.stem
                suffix = target_path.suffix
                idx = 1
                while (dest_dir / f"{stem}_{idx}{suffix}").exists():
                    idx += 1
                target_path = dest_dir / f"{stem}_{idx}{suffix}"

            shutil.copy2(pdb_path, target_path)
            num_copied += 1

    print(f"Copied {num_copied} PDB file(s) to {dest_dir}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Copy successful PDBs for a job")
    parser.add_argument("job_id", help="Job identifier (e.g. 1747410..._job)")
    args = parser.parse_args()

    copy_pdbs(args.job_id)


if __name__ == "__main__":
    main() 