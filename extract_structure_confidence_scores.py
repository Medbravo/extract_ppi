import json
import numpy as np
import sys
import os
from extract_cif_sequence import extract_sequence_from_cif

def calculate_complex_confidence(json_path, chain_lengths):
    """
    Calculates average pLDDT and average inter-chain pAE for a protein complex.

    Args:
        json_path (str): Path to the JSON file containing prediction scores.
        chain_lengths (list): A list of integers representing the length of
                              each chain.

    Returns:
        tuple: A tuple containing:
               - average_plddt (float): The mean pLDDT score across all residues.
               - average_inter_pae (float): The mean pAE score for inter-chain
                                            residue pairs.
    """
    total_residues = sum(chain_lengths)
    
    try:
        with open(json_path, 'r') as f:
            data = json.load(f)
    except Exception as e:
        print(f"Error reading JSON file: {e}")
        return None, None
    
    plddt_scores = data['atom_plddts']
    
    
    if plddt_scores is None:
        print("Could not find pLDDT scores in the JSON file")
        return None, None
    
    # Calculate average pLDDT
    try:
        average_plddt = np.mean(plddt_scores)
    except Exception as e:
        print(f"Error calculating average pLDDT: {e}")
        return None, None
    
    # Get pAE matrix
    pae_matrix = None
    
   
    pae_matrix = np.array(data['pae'])
    
    if pae_matrix is None:
        print("Could not find pAE matrix in the JSON file")
        return average_plddt, None
    
    # Calculate inter-chain pAE
    try:
        inter_chain_pae_values = []
        current_offset_row = 0
        
        for i, len_chain_i in enumerate(chain_lengths):
            start_row = current_offset_row
            end_row = current_offset_row + len_chain_i
            
            current_offset_col = 0
            for j, len_chain_j in enumerate(chain_lengths):
                start_col = current_offset_col
                end_col = current_offset_col + len_chain_j
                
                # Only consider pairs from different chains
                if i != j and end_row <= pae_matrix.shape[0] and end_col <= pae_matrix.shape[1]:
                    inter_block = pae_matrix[start_row:end_row, start_col:end_col]
                    inter_chain_pae_values.extend(inter_block.flatten().tolist())
                
                current_offset_col += len_chain_j
            current_offset_row += len_chain_i
        
        # Calculate average inter-chain pAE
        if not inter_chain_pae_values:
            if len(chain_lengths) <= 1:
                print("Note: Only one chain detected. Inter-chain pAE is not applicable.")
                average_inter_pae = np.nan
            else:
                print("Warning: No inter-chain pAE values collected.")
                average_inter_pae = np.nan
        else:
            average_inter_pae = np.mean(inter_chain_pae_values)
        
        return average_plddt, average_inter_pae
    
    except Exception as e:
        print(f"Error calculating inter-chain pAE: {e}")
        return average_plddt, None

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python extract_structure_confidence_scores.py <input_directory>")
        print("The directory should contain *_model.cif and *_confidences.json files.")
        sys.exit(1)
    
    # Get input directory
    input_dir = sys.argv[1]
    
    try:
        files = os.listdir(input_dir)
    except Exception as e:
        print(f"Error accessing directory {input_dir}: {e}")
        sys.exit(1)
    
    # Find specific files with expected naming pattern
    cif_file = None
    json_file = None
    
    # Get files ending with _model.cif and _confidences.json
    for file in files:
        if file.endswith('_model.cif'):
            cif_file = os.path.join(input_dir, file)
        elif file.endswith('_confidences.json'):
            json_file = os.path.join(input_dir, file)
    
    # Check if files were found
    if cif_file is None:
        print(f"Error: No *_model.cif file found in {input_dir}")
        sys.exit(1)
    
    if json_file is None:
        print(f"Error: No *_confidences.json file found in {input_dir}")
        sys.exit(1)
    
    print(f"Processing directory: {input_dir}")
    print(f"CIF file: {os.path.basename(cif_file)}")
    print(f"JSON file: {os.path.basename(json_file)}")
    
    # Extract chain information from CIF file
    chain_sequences = extract_sequence_from_cif(cif_file)
    
    if not chain_sequences:
        print("Error: Could not extract chain information from CIF file.")
        sys.exit(1)
    
    # Get chain lengths in sorted order by chain ID
    chain_lengths = [length for _, (_, length) in sorted(chain_sequences.items())]
    chain_ids = sorted(chain_sequences.keys())
    
    print(f"Detected chains: {chain_ids}")
    print(f"Chain lengths: {chain_lengths} (Total residues: {sum(chain_lengths)})")
    
    # Calculate confidence scores
    avg_plddt, avg_inter_pae = calculate_complex_confidence(json_file, chain_lengths)
    
    # Print results
    if avg_plddt is not None:
        print(f"\nAverage pLDDT for the complex: {avg_plddt:.2f}")
    else:
        print("\nCould not calculate Average pLDDT.")
    
    if avg_inter_pae is not None:
        if np.isnan(avg_inter_pae):
            print("Average Inter-Chain pAE: Not Applicable (only one chain)")
        else:
            print(f"Average Inter-Chain pAE: {avg_inter_pae:.2f} Angstroms")
    else:
        print("Could not calculate Average Inter-Chain pAE.")