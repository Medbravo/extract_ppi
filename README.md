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

Extract protein-protein interfaces from a PDB structure. If not in the pdb_dir, the script will download it.

Example:

```bash
python extract_ppi.py --pdb_id 2gs6 --pdb_dir in/ --distance 5
# Extract interfaces from PDB ID 1NQL
python extract_ppi.py --pdb_id 1NQL --distance 2.5

# With contact probability filtering
python extract_ppi.py --input 1ABC --distance 5.0 \
    --contact_probs_json path/to/confidences.json \
    --threshold_probability_contact 0.8

# Directory input with automatic confidences detection
python extract_ppi.py --input /path/to/alphafold3_output/ \
    --distance 5.0 --threshold_probability_contact 0.8
```

### Arguments

- `--pdb_id`: The 4-character PDB identifier (required)
- `--input`: PDB identifier, path to a PDB/CIF file, or directory containing CIF and optionally confidences.json

### Options

- `--pdb_dir PATH`: Directory to store PDB files or where it is expected to have the PDB (default: ./pdbs)
- `--distance FLOAT`: Distance threshold for interface detection in Angstroms (default: 2.5)
- `--contact_probs_json PATH`: Path to JSON file containing contact probabilities for filtering interface residues
- `--threshold_probability_contact FLOAT`: Threshold for contact probability filtering (default: 0.9)

## Enhanced Contact Probability Filtering

The tool now supports filtering interface residues based on contact probabilities from AlphaFold3 or similar predictions:

### Key Features
- **Automatic detection**: When using directory input, automatically detects and loads `confidences.json` files
- **Precise filtering**: Only residues with contact probability above threshold are selected (no window expansion)
- **Multiple input formats**: Supports PDB IDs, file paths, or directories containing CIF and confidences files

### Directory Input
When providing a directory path, the tool automatically:
1. Searches for CIF files in the directory
2. Looks for confidences files (patterns: `confidences.json`, `*confidences.json`, `*_confidences.json`)
3. Enables contact probability filtering if confidences file is found

### Contact Probability Logic
- Calculates maximum contact probability for each residue with partner chain residues
- Filters residues where `max_contact_prob > threshold`
- **No window expansion** (window_size = 0) for precise interface detection
- Combines distance-based and probability-based filtering for accurate results

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

## Output Format

The `extract_ppi` tool outputs interface information in the following format:

```
>{CHAIN1}_{CHAIN2}|{PDB_ID}|DATE_{DATE}|SCORE_{SCORE}|METHOD_{METHOD}|{CHAIN1_HOTSPOTS}|{CHAIN2_HOTSPOTS}|{CHAIN1}_{COUNT1}|{CHAIN2}_{COUNT2}
```

### Example Output:
```
>A_B|1747690647190-4494-65a0dd34_model|DATE_|SCORE_NA|METHOD_|ILE13 ARG17 ASN18 PRO20 ASN23 GLY24 TRP25 MET27 THR28 ILE31 ALA32 PHE35|LEU349 HIS370 LEU372 PRO373 VAL374 LEU406 GLN408 ALA409 HIS433 GLY434 PHE436 ALA439 VAL441 SER442 ILE462 SER464 LYS489 ILE491 SER492|A_12|B_19
```

### Field Descriptions:
- `CHAIN1_CHAIN2`: Interface identifier (e.g., A_B)
- `PDB_ID`: The structure identifier or filename
- `DATE_`: Deposition date (if available from PDB header)
- `SCORE_`: Quality score/resolution (if available from PDB header)
- `METHOD_`: Experimental method used to determine the structure (X-RAY, NMR, etc.)
- `CHAIN1_HOTSPOTS`: Space-separated list of interface residues from chain 1 (format: RES_NAME + RES_NUMBER)
- `CHAIN2_HOTSPOTS`: Space-separated list of interface residues from chain 2 (format: RES_NAME + RES_NUMBER)
- `CHAIN1_COUNT`: Chain 1 identifier and number of hotspot residues (e.g., A_12)
- `CHAIN2_COUNT`: Chain 2 identifier and number of hotspot residues (e.g., B_19)

### Pipeline Output Format

For the pipeline complex retrieval, the output is a FASTA file named `{target_name}_protein_complexes.fasta` with entries in the format:

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
