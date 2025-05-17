#!/bin/bash
#SBATCH --job-name=relax_metrics
#SBATCH --nodes=1
#SBATCH --time=1:00:00                 # adjust wall-clock limit as needed
#SBATCH --output=exelogs/out/relax_metrics.%A_%a.out
#SBATCH --error=exelogs/err/relax_metrics.%A_%a.err
# >>> IMPORTANT <<<
# Set the array range when you submit, e.g.:
#   sbatch --array=1-$(ls /path/to/pdbs/*.pdb | wc -l) relax_metrics_array.slurm
# Or hard-code:   #SBATCH --array=1-100
# Or 
# NUM=$(find /home/auri/extract_ppi/in/EGFR_ECD/1746470428303-2758-312e7774_job_5XWD_24_affibodies_any_hotspot \
#             -maxdepth 1 -type f \( -name "*.pdb" -o -name "*.pdb.gz" \) | wc -l)
# sbatch --array=1-"$NUM" relax_metrics_array.slurm
#-------------------------------------------------------------------------------

set -euo pipefail

###########################################################################
# USER CONFIG (edit paths below)
###########################################################################
ROSETTA_BIN="/home/auri/Rosetta/rosetta.source.release-371/main/source/bin/rosetta_scripts.linuxgccrelease"
ROSETTA_DB="/home/auri/Rosetta/rosetta.source.release-371/main/database"
DALPHABALL_EXEC="/home/auri/Rosetta/rosetta.source.release-371/main/source/external/DAlphaBall/DAlphaBall.gcc"
XML_SCRIPT="/home/auri/Rosetta/relax_and_metrics.xml"           # full path so worker node can find it

# PDB directory is now provided as the first positional argument when you
# submit the job, e.g.:
#   sbatch --array=1-$(ls my_pdbs/*.pdb | wc -l) relax_metrics_array.slurm my_pdbs
#
# Exit with a helpful message if missing.
if [[ $# -lt 1 ]]; then
    echo "Usage: sbatch [--array=...] relax_metrics_array.slurm <PDB_DIR>" >&2
    exit 1
fi

PDB_DIR="$1"                                 # directory with input complexes
OUTPUT_ROOT="$PDB_DIR/relax_metrics_$(date +%Y%m%d)"  # all runs for this submission
NSTRUCT=1                 # number of output models per input structure

###########################################################################

# Build an ordered list of PDB (and gzipped PDB) files ---------------------
mapfile -t PDB_LIST < <(find "$PDB_DIR" -maxdepth 1 -type f \( -name "*.pdb" -o -name "*.pdb.gz" \) | sort)
TOTAL_FILES=${#PDB_LIST[@]}

if (( TOTAL_FILES == 0 )); then
    echo "[ERROR] No PDB files found in $PDB_DIR" >&2
    exit 1
fi

TASK_ID=${SLURM_ARRAY_TASK_ID:-1}
if (( TASK_ID < 1 || TASK_ID > TOTAL_FILES )); then
    echo "[ERROR] SLURM_ARRAY_TASK_ID ($TASK_ID) out of range 1-$TOTAL_FILES" >&2
    exit 1
fi

PDB_FILE="${PDB_LIST[$((TASK_ID-1))]}"
BASE_NAME=$(basename "$PDB_FILE")
BASE_NAME=${BASE_NAME%%.pdb}
BASE_NAME=${BASE_NAME%%.pdb.gz}

CURRENT_OUT="${OUTPUT_ROOT}/${BASE_NAME}"
mkdir -p "$CURRENT_OUT"

echo "[INFO] Task $TASK_ID processing $PDB_FILE → $CURRENT_OUT";

# Run Rosetta ----------------------------------------------------------------
"$ROSETTA_BIN" \
    -database "$ROSETTA_DB" \
    -parser:protocol "$XML_SCRIPT" \
    -s "$PDB_FILE" \
    -nstruct "$NSTRUCT" \
    -out:path:all "$CURRENT_OUT" \
    -out:prefix "${BASE_NAME}_" \
    -out:file:scorefile "${BASE_NAME}.sc" \
    -holes:dalphaball "$DALPHABALL_EXEC" \
    -run:preserve_header \
    -overwrite \
    -beta_nov16  # remove or replace with other options you need

# End of task
echo "[DONE] $PDB_FILE" 