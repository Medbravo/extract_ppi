from Bio import PDB
import os
import glob
import argparse

def convert_cif_to_pdb(cif_path, pdb_path):
    """Convert a .cif file to .pdb format"""
    parser = PDB.MMCIFParser()
    io = PDB.PDBIO()
    
    # Parse the CIF file
    structure = parser.get_structure('structure', cif_path)
    
    # Save as PDB
    io.set_structure(structure)
    io.save(pdb_path)

def batch_convert_folder(input_folder, output_folder=None):
    """Recursively convert all .cif files found under *input_folder* to .pdb format.
    The resulting .pdb files are saved in the root of *input_folder* (or *output_folder* if provided)."""
    # If no output folder specified, use input folder
    if output_folder is None:
        output_folder = input_folder
        
    # Create output folder if it doesn't exist
    os.makedirs(output_folder, exist_ok=True)
    
    # Find all .cif files recursively
    cif_files = glob.glob(os.path.join(input_folder, "**", "*.cif"), recursive=True)
    
    for cif_file in cif_files:
        print(cif_file)
        # Generate output PDB filename
        base_name = os.path.splitext(os.path.basename(cif_file))[0]
        pdb_file = os.path.join(output_folder, f"{base_name}.pdb")
        
        try:
            convert_cif_to_pdb(cif_file, pdb_file)
            print(f"Converted {cif_file} to {pdb_file}")
        except Exception as e:
            print(f"Error converting {cif_file}: {str(e)}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Recursively convert all .cif files under the given input folder(s) to .pdb format."
    )
    parser.add_argument(
        "input_folders",
        nargs="+",
        help="One or more directories to search recursively for .cif files."
    )
    parser.add_argument(
        "-o",
        "--output_folder",
        default=None,
        help="Destination directory for converted .pdb files. If omitted, each input folder is used as its own output directory."
    )

    args = parser.parse_args()

    for input_folder in args.input_folders:
        batch_convert_folder(input_folder, args.output_folder)
