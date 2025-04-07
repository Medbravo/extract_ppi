#!/usr/bin/env python
"""
Pipeline Complex Retrieval

This script retrieves protein binders from PDB for a given target protein
and extracts protein-protein interfaces.
"""

import os
import sys
import argparse
from pathlib import Path
from extract_ppi import retrieve_binders_fromPDB, confirm_experimental_binding_and_add_interface

def run_pipeline(target_seq_file, target_name, pdb_seqres_db, pdb_dir='./pdbs', 
                min_seq_id=0.5, distance_threshold=2.5):
    """
    Run the complete pipeline for retrieving protein complexes and extracting interfaces.
    
    Args:
        target_seq_file: Path to the cropped target sequence file (fasta)
        target_name: Name of the target protein (e.g., 'egfr')
        pdb_seqres_db: Path to the mmseqs database of PDB sequences
        pdb_dir: Directory to store PDB files
        min_seq_id: Minimum sequence identity threshold for mmseqs
        distance_threshold: Distance threshold for interface detection in Angstroms
        
    Returns:
        Path to the final fasta file containing binders with interface information
    """
    print(f"Step 1: Retrieving binders from PDB for {target_name}...")
    binders_file = retrieve_binders_fromPDB(
        target_seq_file,
        target_name,
        pdb_seqres_db,
        min_seq_id=min_seq_id
    )
    
    if not binders_file:
        print("Failed to retrieve binders. Exiting.")
        return None
    
    print(f"Step 2: Confirming experimental binding and extracting interfaces...")
    result_file = confirm_experimental_binding_and_add_interface(
        target_name,
        pdb_dir=pdb_dir,
        distance_threshold=distance_threshold
    )
    
    print(f"Pipeline completed successfully. Results saved to: {result_file}")
    return result_file

def main():
    parser = argparse.ArgumentParser(description='Pipeline for retrieving protein complexes and extracting interfaces')
    parser.add_argument('--target_seq', type=str, required=True,
                        help='Path to the cropped target sequence file (fasta)')
    parser.add_argument('--target_name', type=str, required=True,
                        help='Name of the target protein (e.g., egfr)')
    parser.add_argument('--pdb_seqres_db', type=str, required=True,
                        help='Path to the mmseqs database of PDB sequences')
    parser.add_argument('--pdb_dir', type=str, default='./pdbs',
                        help='Directory to store PDB files (default: ./pdbs)')
    parser.add_argument('--min_seq_id', type=float, default=0.5,
                        help='Minimum sequence identity for mmseqs (default: 0.5)')
    parser.add_argument('--distance', type=float, default=2.5,
                        help='Distance threshold for interface detection in Angstroms (default: 2.5)')
    
    args = parser.parse_args()
    
    run_pipeline(
        args.target_seq,
        args.target_name,
        args.pdb_seqres_db,
        pdb_dir=args.pdb_dir,
        min_seq_id=args.min_seq_id,
        distance_threshold=args.distance
    )

if __name__ == "__main__":
    main()