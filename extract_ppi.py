#!/usr/bin/env python
import os
import sys
import argparse
from pathlib import Path
import numpy as np
from collections import defaultdict, Counter
from scipy.spatial import cKDTree
from Bio.PDB import PDBParser, Entity, MMCIFParser
from Bio.PDB.PDBList import PDBList
from Bio import SeqIO
import subprocess
import re
import datetime
import json

def download_pdb(pdb_id, pdb_dir):
    """Download PDB file from the PDB database if not already present."""
    # Check if file already exists with any of the common extensions
    possible_filenames = [
        f"{pdb_id}.pdb",
        f"{pdb_id}.ent",
        f"{pdb_id}.cif",
        f"{pdb_id.lower()}.ent",
        f"{pdb_id.lower()}.pdb",
        f"{pdb_id.lower()}.cif"
    ]
    
    for filename in possible_filenames:
        file_path = Path(pdb_dir) / filename
        if file_path.exists():
            print(f"Found existing PDB file: {file_path}")
            return str(file_path)
    
    # Download if not exists
    print(f"Downloading PDB structure '{pdb_id}'...")
    pdb_list = PDBList()
    filename = pdb_list.retrieve_pdb_file(pdb_id, file_format="pdb", pdir=pdb_dir)
    return filename

def get_experimental_method(structure):
    """Get the experimental method from the PDB file."""
    header = structure.header
    if 'structure_method' in header:
        return header['structure_method']
    return "UNKNOWN"

def get_deposition_date(structure):
    """Get the deposition date from the PDB file."""
    header = structure.header
    if 'deposition_date' in header:
        return header['deposition_date']
    return "UNKNOWN"

def get_quality_score(structure):
    """Get quality score from the PDB file (resolution for X-ray, etc.)"""
    header = structure.header
    if 'resolution' in header and header['resolution'] is not None:
        return f"{header['resolution']:.2f}"
    elif 'resolution' in header:
        return "NA"
    return "NA"

def extract_ppi(pdb_file, distance_threshold=5, contact_probs_json=None, threshold_probability_contact=0.7):
    """Extract protein-protein interface residues based on distance threshold and optionally filter by contact probabilities."""
    # Use the appropriate parser based on file extension
    if pdb_file.lower().endswith('.cif'):
        parser = MMCIFParser(QUIET=True)
    else:
        parser = PDBParser(QUIET=True)
        
    # Load contact probability data if provided
    contact_prob_data = None
    if contact_probs_json:
        try:
            with open(contact_probs_json, 'r') as f:
                contact_prob_data = json.load(f)
        except Exception as e:
            print(f"Warning: Could not load contact probabilities from {contact_probs_json}: {e}")
            contact_prob_data = None
        
    try:
        structure = parser.get_structure("complex", pdb_file)
        model = structure[0]
        
        # Get experimental method
        exp_method = get_experimental_method(structure)
        
        # Get deposition date and quality score
        deposition_date = get_deposition_date(structure)
        quality_score = get_quality_score(structure)
        
        # Create a dictionary to hold chains and their atoms
        chain_atoms = {}
        chain_residues = {}
        
        # Collect atoms and residues for each chain
        for chain in model:
            chain_id = chain.get_id()
            chain_atoms[chain_id] = []
            chain_residues[chain_id] = {}
            
            for residue in chain:
                # Skip non-amino acid residues (like water or ligands)
                if residue.get_id()[0] != " ":
                    continue
                
                res_id = residue.get_id()[1]
                chain_residues[chain_id][res_id] = residue
                
                for atom in residue:
                    chain_atoms[chain_id].append((atom, residue))
        
        # Find interfaces between all chains
        chain_interfaces = {}
        
        for chain1_id in chain_atoms:
            for chain2_id in chain_atoms:
                if chain1_id >= chain2_id:  # Skip self-interactions and duplicates
                    continue
                    
                interface_key = f"{chain1_id}_{chain2_id}"
                chain_interfaces[interface_key] = {"chain1": set(), "chain2": set()}
                
                # Create arrays of atom coordinates for each chain
                atoms1 = [(atom.get_coord(), i) for i, (atom, _) in enumerate(chain_atoms[chain1_id])]
                atoms2 = [(atom.get_coord(), i) for i, (atom, _) in enumerate(chain_atoms[chain2_id])]
                
                if not atoms1 or not atoms2:
                    continue
                    
                coords1 = np.array([a[0] for a in atoms1])
                coords2 = np.array([a[0] for a in atoms2])
                
                # Use KDTree to find all pairs of atoms within the distance threshold
                tree1 = cKDTree(coords1)
                tree2 = cKDTree(coords2)
                
                pairs = tree1.query_ball_tree(tree2, distance_threshold)
                
                # Process each interaction
                for i, close_indices in enumerate(pairs):
                    if close_indices:  # If there are interactions
                        atom1_idx = atoms1[i][1]
                        residue1 = chain_atoms[chain1_id][atom1_idx][1]
                        res1_id = residue1.get_id()[1]
                        chain_interfaces[interface_key]["chain1"].add(res1_id)
                        
                        for j in close_indices:
                            atom2_idx = atoms2[j][1]
                            residue2 = chain_atoms[chain2_id][atom2_idx][1]
                            res2_id = residue2.get_id()[1]
                            chain_interfaces[interface_key]["chain2"].add(res2_id)
        
        # Apply contact probability filtering if data is available
        if contact_prob_data:
            chain_interfaces = filter_interfaces_by_contact_probs(
                chain_interfaces, contact_prob_data, threshold_probability_contact
            )
        
        # Format the interfaces for output
        formatted_interfaces = []
        
        for interface_key, interface in chain_interfaces.items():
            chain1_id, chain2_id = interface_key.split("_")
            
            # Format the interface residues for chain 1
            chain1_hotspots = []
            for res_id in sorted(interface["chain1"]):
                residue = chain_residues[chain1_id][res_id]
                res_name = residue.get_resname()
                chain1_hotspots.append(f"{res_name}{res_id}")
            
            # Format the interface residues for chain 2
            chain2_hotspots = []
            for res_id in sorted(interface["chain2"]):
                residue = chain_residues[chain2_id][res_id]
                res_name = residue.get_resname()
                chain2_hotspots.append(f"{res_name}{res_id}")
            
            # Add deposition date and quality score to the output - don't add prefixes here, they'll be added later
            formatted_interface = f"{chain1_id}_{chain2_id}|{pdb_file}|{deposition_date}|{quality_score}|METHOD_{exp_method}|{' '.join(chain1_hotspots)}|{' '.join(chain2_hotspots)}"
            formatted_interfaces.append(formatted_interface)
        
        return formatted_interfaces
    except Exception as e:
        print(f"Error parsing file {pdb_file}: {e}")
        return []

def filter_interfaces_by_contact_probs(chain_interfaces, contact_prob_data, threshold_probability_contact):
    """Filter interface residues based on contact probabilities from JSON data."""
    try:
        contact_probs = np.array(contact_prob_data['contact_probs'])
        token_chain_ids = contact_prob_data['token_chain_ids']
        token_res_ids = contact_prob_data['token_res_ids']
        
        # Find unique chains in the contact probability data
        unique_chains = sorted(set(token_chain_ids))
        if len(unique_chains) < 2:
            print("Warning: Less than 2 chains found in contact probability data")
            return chain_interfaces
            
        # Create masks for each chain
        chain_masks = {}
        chain_res_ids = {}
        for chain_id in unique_chains:
            mask = np.array(token_chain_ids) == chain_id
            chain_masks[chain_id] = mask
            chain_res_ids[chain_id] = np.array(token_res_ids)[mask]
        
        # Filter each interface
        filtered_interfaces = {}
        
        for interface_key, interface in chain_interfaces.items():
            chain1_id, chain2_id = interface_key.split("_")
            
            # Map PDB chain IDs to contact probability chain IDs
            # This assumes the chain IDs match, but you might need to adjust this mapping
            prob_chain1_id = None
            prob_chain2_id = None
            
            # Try to find matching chains in contact probability data
            for prob_chain_id in unique_chains:
                if prob_chain_id == chain1_id or (len(unique_chains) == 2 and prob_chain1_id is None):
                    prob_chain1_id = prob_chain_id
                elif prob_chain_id == chain2_id or (len(unique_chains) == 2 and prob_chain2_id is None):
                    prob_chain2_id = prob_chain_id
            
            # If we have exactly 2 chains and haven't assigned them yet, assign them
            if len(unique_chains) == 2 and prob_chain1_id is None:
                prob_chain1_id = unique_chains[0]
                prob_chain2_id = unique_chains[1]
            elif len(unique_chains) == 2 and prob_chain2_id is None:
                prob_chain2_id = unique_chains[1] if prob_chain1_id == unique_chains[0] else unique_chains[0]
            
            if prob_chain1_id is None or prob_chain2_id is None:
                print(f"Warning: Could not map chains {chain1_id}, {chain2_id} to contact probability data")
                filtered_interfaces[interface_key] = interface
                continue
            
            # Get masks and residue IDs for both chains
            mask1 = chain_masks[prob_chain1_id]
            mask2 = chain_masks[prob_chain2_id]
            res_ids1 = chain_res_ids[prob_chain1_id]
            res_ids2 = chain_res_ids[prob_chain2_id]
            
            # Calculate maximum contact probabilities for each chain
            max_probs1 = np.max(contact_probs[mask1][:, mask2], axis=1)
            max_probs2 = np.max(contact_probs[mask2][:, mask1], axis=1)
            
            # Find residues above threshold
            interface_indices1 = np.where(max_probs1 > threshold_probability_contact)[0]
            interface_indices2 = np.where(max_probs2 > threshold_probability_contact)[0]
            
            # Get the actual residue IDs that pass the threshold
            filtered_res_ids1 = set(res_ids1[interface_indices1])
            filtered_res_ids2 = set(res_ids2[interface_indices2])
            
            # Apply window expansion (similar to extract_interface_indices)
            window_size = 0 #restrict to no window expansion
            expanded_res_ids1 = set()
            expanded_res_ids2 = set()
            
            for res_id in filtered_res_ids1:
                for w in range(res_id - window_size, res_id + window_size + 1):
                    if min(res_ids1) <= w <= max(res_ids1):
                        expanded_res_ids1.add(w)
            
            for res_id in filtered_res_ids2:
                for w in range(res_id - window_size, res_id + window_size + 1):
                    if min(res_ids2) <= w <= max(res_ids2):
                        expanded_res_ids2.add(w)
            
            # Filter the original interface residues
            filtered_chain1 = interface["chain1"].intersection(expanded_res_ids1)
            filtered_chain2 = interface["chain2"].intersection(expanded_res_ids2)
            
            filtered_interfaces[interface_key] = {
                "chain1": filtered_chain1,
                "chain2": filtered_chain2
            }
        
        return filtered_interfaces
        
    except Exception as e:
        print(f"Error filtering interfaces by contact probabilities: {e}")
        return chain_interfaces

def retrieve_binders_fromPDB(target_seq_file, target_name, pdb_seqres_db, min_seq_id=0.5, cov_mode=2, cov_threshold=0):
    """
    Retrieve all protein binders for a given target from PDB.
    
    Args:
        target_seq_file: Path to the cropped target sequence file (fasta)
        target_name: Name of the target protein (e.g., 'egfr')
        pdb_seqres_db: Path to the mmseqs database of PDB sequences
        min_seq_id: Minimum sequence identity threshold (default: 0.5)
        cov_mode: Coverage mode for mmseqs (default: 2)
        cov_threshold: Coverage threshold (default: 0)
        
    Returns:
        Path to the output fasta file containing all potential binder sequences
    """
    # Create output directory if it doesn't exist
    out_dir = Path('out')
    out_dir.mkdir(exist_ok=True)
    
    # Create temp directory if it doesn't exist
    tmp_dir = Path('tmp')
    tmp_dir.mkdir(exist_ok=True)
    
    # Output file paths
    result_file = out_dir / 'result.m8'
    output_fasta = out_dir / f'{target_name}_protein_complexes.fasta'
    
    # Run mmseqs easy-search to find similar sequences
    cmd = [
        'mmseqs', 'easy-search', 
        str(target_seq_file), 
        str(pdb_seqres_db), 
        str(result_file), 
        str(tmp_dir),
        '--min-seq-id', str(min_seq_id),
        '--cov-mode', str(cov_mode),
        '-c', str(cov_threshold)
    ]
    
    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running mmseqs: {e}")
        return None
    
    # Parse the result file to get target chains
    target_chains = set()
    with open(result_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                target_chains.add(parts[1])  # Format: 1IVO_A
    
    # Parse the PDB sequence database to extract binders
    pdb_seqres_file = str(pdb_seqres_db).replace('_db', '.txt')
    
    # Dictionary to keep track of unique sequences and their header information
    unique_sequences = {}
    chain_headers = {}
    
    # First pass: collect all chains, their sequences, and header information
    with open(pdb_seqres_file, 'r') as f:
        current_id = None
        current_header = ""
        current_seq = ""
        
        for line in f:
            if line.startswith('>'):
                if current_id and current_seq:
                    pdb_chain_id = current_id.split()[0]  # Get just the chain ID part (e.g., "1IVO_A")
                    unique_sequences[pdb_chain_id] = current_seq
                    chain_headers[pdb_chain_id] = current_header
                
                current_header = line.strip()[1:]  # Store the full header without '>'
                current_id = current_header.split()[0]  # Extract PDB ID and chain (e.g., "1IVO_A")
                current_seq = ""
            else:
                current_seq += line.strip()
        
        # Add the last sequence
        if current_id and current_seq:
            pdb_chain_id = current_id.split()[0]
            unique_sequences[pdb_chain_id] = current_seq
            chain_headers[pdb_chain_id] = current_header
    
    # Extract PDB IDs that contain target chains
    pdb_ids_with_target = set()
    for chain_id in target_chains:
        pdb_id = chain_id.split('_')[0]
        pdb_ids_with_target.add(pdb_id)
    
    # Dictionary to track unique sequences for output
    output_sequences = {}
    
    # Second pass: write binder sequences to output file
    with open(output_fasta, 'w') as outfile:
        for chain_id, sequence in unique_sequences.items():
            pdb_id, chain_num = chain_id.split('_')
            
            # Check if this PDB has at least one chain matching the target
            if pdb_id in pdb_ids_with_target:
                # Skip chains that match the target
                if chain_id in target_chains:
                    continue
                
                # Check if sequence is already in our output collection
                if sequence in output_sequences:
                    continue
                
                # Add to output sequences
                output_sequences[sequence] = chain_id
                
                # Get the original header information
                header_info = chain_headers[chain_id]
                
                # Extract description part (everything after the first space)
                if ' ' in header_info:
                    description = header_info[header_info.index(' ')+1:]
                else:
                    description = ""
                
                # Write to output file with description preserved
                if description:
                    outfile.write(f">{pdb_id}_{chain_num}|{pdb_id}|{chain_num} {description}\n")
                else:
                    outfile.write(f">{pdb_id}_{chain_num}|{pdb_id}|{chain_num}\n")
                outfile.write(f"{sequence}\n")
    
    return str(output_fasta)

def confirm_experimental_binding_and_add_interface(target_name, pdb_dir='./pdbs', distance_threshold=2.5, contact_probs_json=None, threshold_probability_contact=0.7):
    """
    Confirm experimental binding and add interface information to the protein complexes file.
    
    Args:
        target_name: Name of the target protein (e.g., 'egfr')
        pdb_dir: Directory to store PDB files
        distance_threshold: Distance threshold for interface detection in Angstroms
        contact_probs_json: Path to JSON file containing contact probabilities for filtering
        threshold_probability_contact: Threshold for contact probability filtering
        
    Returns:
        Path to the updated fasta file
    """
    # Create directories if they don't exist
    out_dir = Path('out')
    out_dir.mkdir(exist_ok=True)
    
    pdb_dir_path = Path(pdb_dir)
    pdb_dir_path.mkdir(exist_ok=True)
    
    # Input and output file paths
    input_fasta = out_dir / f'{target_name}_protein_complexes.fasta'
    
    # Dictionary to store updated entries
    updated_entries = {}
    
    # Read the input fasta file and store descriptions
    records = []
    descriptions = {}
    protein_names = {}
    
    with open(input_fasta, 'r') as f:
        for line in f:
            if line.startswith('>'):
                header = line.strip()[1:]  # Remove the '>' character
                record_id = header.split()[0]  # Get the ID part before the first space
                
                # Store description (everything after the ID)
                if ' ' in header:
                    full_desc = header[header.index(' ')+1:]
                    descriptions[record_id] = full_desc
                    
                    # Extract protein name (last part after "length:XX")
                    if "length:" in full_desc and "  " in full_desc:
                        parts = full_desc.split("  ")
                        if len(parts) > 1:
                            protein_names[record_id] = parts[-1].strip()
                        else:
                            protein_names[record_id] = ""
                    else:
                        protein_names[record_id] = ""
                else:
                    descriptions[record_id] = ""
                    protein_names[record_id] = ""
                
                record = {'id': record_id, 'seq': ''}
                records.append(record)
            else:
                records[-1]['seq'] += line.strip()
    
    # Process each record
    for record in records:
        record_id = record['id']
        sequence = record['seq']
        description = descriptions.get(record_id, "")
        protein_name = protein_names.get(record_id, "")
        
        # Parse the record ID to get PDB ID
        parts = record_id.split('|')
        pdb_chain = parts[0]
        pdb_id = parts[1]
        chain_id = parts[2]
        
        try:
            # Download the PDB file
            pdb_file = download_pdb(pdb_id, str(pdb_dir_path))
            
            # Extract interfaces
            interfaces = extract_ppi(pdb_file, distance_threshold=distance_threshold, 
                                    contact_probs_json=contact_probs_json, 
                                    threshold_probability_contact=threshold_probability_contact)
            
            # Find interfaces involving this chain
            interface_found = False
            for interface in interfaces:
                interface_parts = interface.split('|')
                chains = interface_parts[0].split('_')
                
                # Check if this chain is in the interface
                if chain_id in chains:
                    # Get the partner chain
                    partner_chain = chains[1] if chains[0] == chain_id else chains[0]
                    
                    # Get the deposition date and quality score - don't remove prefixes here
                    deposition_date = interface_parts[2]
                    quality_score = interface_parts[3]
                    
                    # Get the experimental method - remove METHOD_ prefix
                    exp_method = interface_parts[4].replace('METHOD_', '')
                    
                    # Get the hotspot residues
                    if chains[0] == chain_id:
                        chain_hotspots = interface_parts[5]
                        partner_hotspots = interface_parts[6]
                    else:
                        chain_hotspots = interface_parts[6]
                        partner_hotspots = interface_parts[5]
                    
                    # Create updated entry with protein name followed by interface information - add prefixes
                    if protein_name:
                        updated_id = f"{pdb_chain}|{pdb_id}|{protein_name}|DATE_{deposition_date}|SCORE_{quality_score}|METHOD_{exp_method}|{partner_hotspots}|{chain_hotspots}"
                    else:
                        updated_id = f"{pdb_chain}|{pdb_id}|DATE_{deposition_date}|SCORE_{quality_score}|METHOD_{exp_method}|{partner_hotspots}|{chain_hotspots}"
                    
                    updated_entries[record_id] = (updated_id, sequence)
                    interface_found = True
                    break
            
            # If no interface was found, keep the original entry with protein name
            if not interface_found:
                if protein_name:
                    updated_entries[record_id] = (f"{pdb_chain}|{pdb_id}|{protein_name}|DATE_UNKNOWN|SCORE_NA", sequence)
                else:
                    updated_entries[record_id] = (f"{pdb_chain}|{pdb_id}|DATE_UNKNOWN|SCORE_NA", sequence)
                
        except Exception as e:
            print(f"Error processing PDB {pdb_id}: {str(e)}")
            # Keep the original entry if there's an error
            if protein_name:
                updated_entries[record_id] = (f"{pdb_chain}|{pdb_id}|{protein_name}|DATE_UNKNOWN|SCORE_NA", sequence)
            else:
                updated_entries[record_id] = (f"{pdb_chain}|{pdb_id}|DATE_UNKNOWN|SCORE_NA", sequence)
    
    # Write updated entries to output file
    with open(input_fasta, 'w') as outfile:
        for original_id, (updated_id, sequence) in updated_entries.items():
            outfile.write(f">{updated_id}\n")
            outfile.write(f"{sequence}\n")
    
    return str(input_fasta)

def main():
    parser = argparse.ArgumentParser(description='Extract protein-protein interfaces from PDB structures')
    
    # Mutually exclusive group for input source
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument('--input', type=str, 
                             help='PDB identifier, path to a PDB/CIF file, or directory containing CIF and optionally confidences.json')
    input_group.add_argument('--pdb_id', type=str, 
                             help='PDB identifier (for backwards compatibility)')

    parser.add_argument('--pdb_dir', type=str, default='./pdbs', 
                        help='Directory to store PDB files if downloading (default: ./pdbs)')
    parser.add_argument('--distance', type=float, default=2.5,
                        help='Distance threshold for interface detection in Angstroms (default: 2.5)')
    parser.add_argument('--contact_probs_json', type=str, default=None,
                        help='Path to JSON file containing contact probabilities for filtering interface residues')
    parser.add_argument('--threshold_probability_contact', type=float, default=0.9,
                        help='Threshold for contact probability filtering (default: 0.9)')
    
    args = parser.parse_args()
    
    # Create PDB directory if it doesn't exist (used for downloads)
    pdb_dir = Path(args.pdb_dir)
    pdb_dir.mkdir(exist_ok=True)
    
    pdb_file = None
    pdb_id_for_reporting = "UNKNOWN"
    input_source_description = ""

    try:
        if args.pdb_id:
            # Handle --pdb_id (backwards compatibility)
            pdb_id = args.pdb_id
            pdb_id_for_reporting = pdb_id
            input_source_description = f"PDB ID: {pdb_id}"
            print(f"Processing PDB ID (using --pdb_id): {pdb_id}")
            # Pass the pdb_dir to the download function
            pdb_file = download_pdb(pdb_id, str(pdb_dir))
        
        elif args.input:
            # Handle --input (file path, directory, or PDB ID)
            input_source_description = f"input: {args.input}"
            input_path = Path(args.input)
            
            if input_path.is_file():
                print(f"Using local PDB file: {input_path}")
                pdb_file = str(input_path)
                # Use the filename stem as the ID for reporting
                pdb_id_for_reporting = input_path.stem
                
            elif input_path.is_dir():
                print(f"Using directory: {input_path}")
                # Look for CIF files in the directory
                cif_files = list(input_path.glob("*.cif"))
                if not cif_files:
                    raise FileNotFoundError(f"No CIF files found in directory: {input_path}")
                
                # Use the first CIF file found
                pdb_file = str(cif_files[0])
                pdb_id_for_reporting = cif_files[0].stem
                print(f"Found CIF file: {pdb_file}")
                
                # Check for confidences.json file in the same directory
                # Try different naming patterns for confidences files
                confidences_patterns = [
                    "confidences.json",
                    "*confidences.json",
                    "*_confidences.json"
                ]
                
                confidences_file = None
                for pattern in confidences_patterns:
                    matching_files = list(input_path.glob(pattern))
                    if matching_files:
                        confidences_file = matching_files[0]
                        break
                
                if confidences_file and args.contact_probs_json is None:
                    args.contact_probs_json = str(confidences_file)
                    print(f"Found confidences file: {confidences_file}")
                    print("Automatically using contact probabilities for filtering")
                
            else:
                # Assume it's a PDB ID and try to download
                pdb_id = args.input
                pdb_id_for_reporting = pdb_id
                print(f"Input '{args.input}' not found as file or directory, assuming PDB ID: {pdb_id}")
                # Pass the pdb_dir to the download function
                pdb_file = download_pdb(pdb_id, str(pdb_dir))
        else:
             # This case should not be reached due to the required mutually exclusive group
             parser.error("Either --input or --pdb_id must be provided.")

        if not pdb_file or not Path(pdb_file).exists():
             raise FileNotFoundError(f"Could not find or download PDB from source: {input_source_description}")

        interfaces = extract_ppi(pdb_file, distance_threshold=args.distance, 
                                 contact_probs_json=args.contact_probs_json, 
                                 threshold_probability_contact=args.threshold_probability_contact)
        
        if interfaces:
            for interface in interfaces:
                # Add DATE_ and SCORE_ prefixes to the fields
                parts = interface.split('|')
                date = parts[2]
                score = parts[3]
                parts[2] = f"DATE_{date}"
                parts[3] = f"SCORE_{score}"
                # Use the determined pdb_id_for_reporting in the output
                parts[1] = pdb_id_for_reporting 
                formatted_interface = "|".join(parts)
                print(f">{formatted_interface}")
        else:
            print(f"No interfaces found in PDB source: {input_source_description}")
            
    except Exception as e:
        print(f"Error processing PDB source '{input_source_description}': {str(e)}")
        sys.exit(1)
        
    
        
if __name__ == "__main__":
    main()