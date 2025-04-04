#!/usr/bin/env python
import os
import sys
import argparse
from pathlib import Path
import numpy as np
from collections import defaultdict
from scipy.spatial import cKDTree
from Bio.PDB import PDBParser, Entity
from Bio.PDB.PDBList import PDBList

def download_pdb(pdb_id, pdb_dir):
    """Download PDB file from the PDB database if not already present."""
    # Check if file already exists
    expected_filename = f"pdb{pdb_id.lower()}.ent"
    file_path = Path(pdb_dir) / expected_filename
    
    if file_path.exists():
        return str(file_path)
    
    # Download if not exists
    pdb_list = PDBList()
    filename = pdb_list.retrieve_pdb_file(pdb_id, file_format="pdb", pdir=pdb_dir)
    return filename

def get_experimental_method(structure):
    """Get the experimental method from the PDB file."""
    header = structure.header
    if 'structure_method' in header:
        return header['structure_method']
    return "UNKNOWN"

def extract_ppi(pdb_file, distance_threshold=2.5):
    """Extract protein-protein interface residues based on distance threshold."""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("complex", pdb_file)
    model = structure[0]
    
    # Get experimental method
    exp_method = get_experimental_method(structure)
    
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
        
        formatted_interface = f"{chain1_id}_{chain2_id}|{pdb_file}|EXP_{exp_method}|{' '.join(chain1_hotspots)}|{' '.join(chain2_hotspots)}"
        formatted_interfaces.append(formatted_interface)
    
    return formatted_interfaces

def main():
    parser = argparse.ArgumentParser(description='Extract protein-protein interfaces from PDB structures')
    parser.add_argument('pdb_id', help='PDB identifier')
    parser.add_argument('--pdb_dir', type=str, default='./pdbs', 
                        help='Directory to store PDB files (default: ./pdbs)')
    parser.add_argument('--distance', type=float, default=4.0,
                        help='Distance threshold for interface detection (default: 4.0)')
    
    args = parser.parse_args()
    
    # Create PDB directory if it doesn't exist
    pdb_dir = Path(args.pdb_dir)
    pdb_dir.mkdir(exist_ok=True)
    
    pdb_id = args.pdb_id.upper()
    
    try:
        # Pass the pdb_dir to the download function
        pdb_file = download_pdb(pdb_id, str(pdb_dir))
        interfaces = extract_ppi(pdb_file, distance_threshold=args.distance)
        
        for interface in interfaces:
            print(f">{interface}")
            
    except Exception as e:
        print(f"Error processing PDB {pdb_id}: {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    main()