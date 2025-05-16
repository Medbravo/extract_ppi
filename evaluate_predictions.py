import argparse
import os
import csv
import json

# Import helper functions from existing module
try:
    # When the script is run from inside the extract_ppi folder
    from extract_structure_confidence_scores import calculate_complex_confidence
    from extract_cif_sequence import extract_sequence_from_cif
except ImportError:
    # Fallback when extract_ppi is installed as a package or cwd is different
    from extract_ppi.extract_structure_confidence_scores import calculate_complex_confidence
    from extract_ppi.extract_cif_sequence import extract_sequence_from_cif


def process_prediction(cif_file: str, json_file: str):
    """Compute average pLDDT and inter-chain pAE for a single prediction.

    Parameters
    ----------
    cif_file : str
        Path to the CIF structure file ending with ``_model.cif``.
    json_file : str
        Path to the confidence JSON file ending with ``_confidences.json``.

    Returns
    -------
    tuple[float | None, float | None]
        ``(avg_plddt, avg_inter_pae)`` – ``None`` if inputs are missing or errors occur.
    """

    # Basic validation
    if not (os.path.isfile(cif_file) and os.path.isfile(json_file)):
        return None, None

    # Extract chain lengths from structure
    chain_sequences = extract_sequence_from_cif(cif_file)
    if not chain_sequences:
        return None, None
    chain_lengths = [length for _, (_, length) in sorted(chain_sequences.items())]

    # Compute confidence metrics
    avg_plddt, avg_inter_pae = calculate_complex_confidence(json_file, chain_lengths)
    return avg_plddt, avg_inter_pae


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Iterate over all sub-directories in a given folder, compute confidence "
            "metrics using extract_structure_confidence_scores.py, and save a summary table."
        )
    )
    parser.add_argument(
        "parent_dir",
        help="Path to folder containing prediction sub-directories (e.g. in/1747032813539-7091-ae5c995d_job)",
    )
    parser.add_argument(
        "-o",
        "--output",
        default="confidence_summary.tsv",
        help="Output TSV file for results (default: confidence_summary.tsv)",
    )
    parser.add_argument(
        "--plddt_threshold",
        type=float,
        default=70.0,
        help="Minimum average pLDDT for SUCCESS (default: 70.0). If your pLDDTs are 0‒1, use 0.7).",
    )
    parser.add_argument(
        "--pae_threshold",
        type=float,
        default=10.0,
        help="Maximum average inter-chain pAE for SUCCESS (default: 10.0 Å)",
    )

    args = parser.parse_args()

    parent_dir = os.path.abspath(args.parent_dir)
    if not os.path.isdir(parent_dir):
        parser.error(f"{parent_dir} is not a directory or does not exist")

    results = []
    # Walk recursively to locate output folders containing both required files
    for root, dirs, files in os.walk(parent_dir):
        if not files:
            continue

        # Get relevant confidence/structure files in the current directory
        cif_file = json_file = summary_file = None
        for file in files:
            if file.endswith('_model.cif'):
                cif_file = os.path.join(root, file)
            elif file.endswith('_confidences.json') and not file.endswith('summary_confidences.json'):
                json_file = os.path.join(root, file)
            elif file.endswith('summary_confidences.json'):
                summary_file = os.path.join(root, file)

        avg_plddt, avg_inter_pae = process_prediction(cif_file, json_file)

        # Extract unique chain_iptm value
        chain_iptm_val = None
        if summary_file and os.path.isfile(summary_file):
            try:
                with open(summary_file) as f:
                    data = json.load(f)
                iptm_list = data.get("chain_iptm", [])
                unique_vals = {v for v in iptm_list if isinstance(v, (int, float))}
                if len(unique_vals) == 1:
                    chain_iptm_val = unique_vals.pop()
            except Exception:
                pass

        # Build FASTA string for all chains in the model
        fasta = None
        chain_sequences = extract_sequence_from_cif(cif_file) if cif_file and os.path.isfile(cif_file) else {}
        if chain_sequences:
            fasta = " ; ".join(
                f">{cid}:{seq}" for cid, (seq, _) in sorted(chain_sequences.items())
            )

        # Determine success status
        success = (
            avg_plddt is not None
            and avg_inter_pae is not None
            and avg_plddt > args.plddt_threshold
            and avg_inter_pae < args.pae_threshold
        )

        rel_path = os.path.relpath(root, parent_dir)

        results.append(
            {
                "directory": rel_path,
                "avg_pLDDT": None if avg_plddt is None else f"{avg_plddt:.2f}",
                "avg_inter_pAE": None if avg_inter_pae is None else f"{avg_inter_pae:.2f}",
                "chain_iptm": None if chain_iptm_val is None else f"{chain_iptm_val:.2f}",
                "SUCCESS": "YES" if success else "NO",
                "fasta": fasta,
            }
        )

        # Print status to console
        status_msg = (
            "SUCCESS" if success else "FAIL" if (avg_plddt is not None and avg_inter_pae is not None) else "ERROR"
        )
        print(f"{rel_path}: {status_msg} (pLDDT={avg_plddt}, pAE={avg_inter_pae})")

    # Write table
    out_path = os.path.abspath(args.output)
    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=["directory", "avg_pLDDT", "avg_inter_pAE", "chain_iptm", "SUCCESS", "fasta"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(results)

    print(f"\nSaved results to {out_path}")


if __name__ == "__main__":
    main() 