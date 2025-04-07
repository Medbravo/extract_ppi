# Extract PPI - Protein-Protein Interface Extractor

A Python tool to extract and analyze protein-protein interfaces from PDB structures and retrieve protein binders for a given target.

## Installation

```bash
# Clone this repository
git clone https://github.com/Medbravo/extract_ppi.git
cd extract_ppi

# Install dependencies
pip install -r requirements.txt
```

## Basic PPI Extraction

Extract protein-protein interfaces from a PDB structure:

```bash
python extract_ppi.py --mode ppi --pdb_id PDB_ID [options]
```

### Arguments

- `--mode ppi`: Run in PPI extraction mode (required)
- `--pdb_id`: The 4-character PDB identifier (required)

### Options

- `--pdb_dir PATH`: Directory to store PDB files (default: ./pdbs)
- `--distance FLOAT`: Distance threshold for interface detection in Angstroms (default: 2.5)

### Example

```bash
# Extract interfaces from PDB ID 1NQL
python extract_ppi.py --mode ppi --pdb_id 1NQL --distance 2.5
```

## Pipeline Complex Retrieval

The pipeline performs the following tasks:
1. Retrieve all protein binders for a given target from PDB
2. Extract protein-protein interfaces for experimentally confirmed complexes

### Running the Pipeline

```bash
python pipeline_complex_retrieval.py --target_seq PATH --target_name NAME --pdb_seqres_db PATH [options]
```

### Arguments

- `--target_seq`: Path to the cropped target sequence (FASTA file) (required)
- `--target_name`: Name of the target protein (e.g., egfr) (required)
- `--pdb_seqres_db`: Path to the mmseqs database of PDB sequences (required)

### Options

- `--pdb_dir PATH`: Directory to store PDB files (default: ./pdbs)
- `--min_seq_id FLOAT`: Minimum sequence identity for mmseqs (default: 0.5)
- `--distance FLOAT`: Distance threshold for interface detection in Angstroms (default: 5)

### Example

```bash
# Run the pipeline for EGFR
python pipeline_complex_retrieval.py \
  --target_seq in/egfr_cropped.fasta \
  --target_name egfr \
  --pdb_seqres_db in/pdb_seqres_db
```

### Running Individual Steps

You can also run individual steps of the pipeline:

#### 1. Retrieve Binders from PDB

```bash
python extract_ppi.py --mode retrieve_binders \
  --target_seq in/egfr_cropped.fasta \
  --target_name egfr \
  --pdb_seqres_db in/pdb_seqres_db
```

#### 2. Confirm Experimental Binding and Extract Interfaces

```bash
python extract_ppi.py --mode confirm_binding \
  --target_name egfr \
  --pdb_dir ./pdbs \
  --distance 2.5
```

## Output Format

The output is a FASTA file named `{target_name}_protein_complexes.fasta` with entries in the format:

```
>{PDB_ID}_{CHAIN}|{PDB_ID}|{CHAIN}|EXP_{METHOD}|{TARGET_HOTSPOTS}|{BINDER_HOTSPOTS}
{PROTEIN_SEQUENCE}
```

Where:
- `PDB_ID`: The 4-character PDB identifier
- `CHAIN`: The chain identifier in the PDB file
- `METHOD`: Experimental method used to determine the structure (X-RAY, NMR, etc.)
- `TARGET_HOTSPOTS`: Space-separated list of residues from the target at the interface
- `BINDER_HOTSPOTS`: Space-separated list of residues from the binder at the interface
- `PROTEIN_SEQUENCE`: The amino acid sequence of the binder

## Dependencies

- numpy
- scipy
- biopython
- mmseqs2 (must be installed and available in PATH)

You can install the Python dependencies using `pip`:

```bash
# Install with pip
pip install -r requirements.txt
```

For mmseqs2 installation, follow the instructions at [MMseqs2 documentation](https://github.com/soedinglab/MMseqs2#installation).