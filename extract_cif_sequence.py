#!/usr/bin/env python3

import os
import sys
import argparse
from Bio.PDB import MMCIFParser

def extract_sequence_from_cif(cif_file):
    """
    Extract sequence and length of each chain from a CIF file.
    
    Args:
        cif_file (str): Path to the CIF file
    
    Returns:
        dict: Dictionary with chain IDs as keys and tuples of (sequence, length) as values
    """
    # Create MMCIF parser
    parser = MMCIFParser()
    
    try:
        # Parse the CIF file
        structure_id = os.path.basename(cif_file).split('.')[0]
        structure = parser.get_structure(structure_id, cif_file)
        
        # Dictionary to store sequences
        chain_sequences = {}
        
        # Map from 3-letter to 1-letter amino acid codes
        three_to_one = {
            'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
            'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
            'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
            'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
            'UNK': 'X'  # Unknown amino acid
        }
        
        # Extract sequence for each chain
        for model in structure:
            for chain in model:
                chain_id = chain.get_id()
                sequence = ""
                
                # Get sequence from residues
                for residue in chain:
                    # Skip hetero-atoms and water
                    if residue.id[0] != " ":
                        continue
                    
                    resname = residue.get_resname()
                    if resname in three_to_one:
                        sequence += three_to_one[resname]
                    else:
                        sequence += "X"  # For non-standard amino acids
                
                # Store sequence and its length
                chain_sequences[chain_id] = (sequence, len(sequence))
        
        return chain_sequences
    
    except Exception as e:
        print(f"Error processing {cif_file}: {str(e)}", file=sys.stderr)
        return {}

def main():
    parser = argparse.ArgumentParser(description='Extract sequence and length from CIF files')
    parser.add_argument('cif_file', help='Path to the CIF file')
    args = parser.parse_args()
    
    # Check if file exists
    if not os.path.exists(args.cif_file):
        print(f"Error: File {args.cif_file} does not exist", file=sys.stderr)
        sys.exit(1)
    
    # Extract sequences
    chain_sequences = extract_sequence_from_cif(args.cif_file)
    
    # Print results
    if chain_sequences:
        print(f"Chain information for {os.path.basename(args.cif_file)}:")
        print("-" * 50)
        for chain_id, (sequence, length) in sorted(chain_sequences.items()):
            print(f"Chain ID: {chain_id}")
            print(f"Length: {length}")
            print(f"Sequence: {sequence}")
            print("-" * 50)
    else:
        print("No chains found or there was an error processing the file.")

if __name__ == "__main__":
    main() 