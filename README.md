# Extract PPI - Protein-Protein Interface Extractor

A Python tool to extract and analyze protein-protein interfaces from PDB structures.

## Installation

```bash
# Clone this repository
git clone https://github.com/Medbravo/extract_ppi.git
cd extract_ppi

# Install dependencies
pip install -r requirements.txt
```

## Usage

```bash
python extract_ppi.py PDB_ID [options]
```

### Arguments

- `PDB_ID`: The 4-character PDB identifier (required)

### Options

- `--pdb_dir PATH`: Directory to store PDB files (default: ./pdbs)
- `--distance FLOAT`: Distance threshold for interface detection in Angstroms (default: 4.0)

### Examples

```bash
# Basic usage with default parameters
python extract_ppi.py 1NQL

# Specify a custom PDB directory and distance threshold
python extract_ppi.py 1NQL --pdb_dir ./my_pdbs --distance 3.5
```

## Output Format

The output is in a simple text format with one line per interface:

```
>CHAIN1_CHAIN2|PDB_FILE|EXP_METHOD|CHAIN1_RESIDUES|CHAIN2_RESIDUES
```

Where:
- `CHAIN1_CHAIN2`: The chain identifiers forming the interface
- `PDB_FILE`: Path to the PDB file
- `EXP_METHOD`: Experimental method used to determine the structure
- `CHAIN1_RESIDUES`: Space-separated list of residues from chain 1 at the interface
- `CHAIN2_RESIDUES`: Space-separated list of residues from chain 2 at the interface

## Dependencies

- numpy
- scipy
- biopython

You can install the dependencies using `uv`:

```bash
# Install with uv (recommended)
uv pip install -r requirements.txt

# Or using standard pip
pip install -r requirements.txt
```