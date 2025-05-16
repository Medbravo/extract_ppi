#!/usr/bin/env bash

################################################################################
# run_relax_metrics.sh
#
# Run Rosetta FastRelax followed by interface analysis (SC) and ddG evaluation
# for every PDB file located in the specified input directory.
#
# Requirements:
#   * Rosetta (compiled) – rosetta_scripts.linuxgccrelease executable
#   * A recent Rosetta database
#   * DAlphaBall executable (only if shape-complementarity requires it)
#
# Edit the USER CONFIG section below to match your installation.
################################################################################

set -euo pipefail

# ------------------- USER CONFIG ---------------------------------------------
ROSETTA_BIN="/home/auri/Rosetta/rosetta.source.release-371/main/source/bin/rosetta_scripts.linuxgccrelease"
ROSETTA_DB="/home/auri/Rosetta/rosetta.source.release-371/main/database"
DALPHABALL_EXEC="/home/auri/Rosetta/rosetta.source.release-371/main/source/external/DAlphaBall/DAlphaBall.gcc"

# PDB_DIR="/home/auri/Rosetta/inputs/1747032813539-7091-ae5c995d_pdbs"   # directory containing input complexes
# PDB_DIR="/home/auri/Rosetta/inputs/1744796811477-1950-0f55c863"   # directory containing input complexes
# PDB_DIR="/home/auri/Rosetta/inputs/task_1746537047984-8610-94b7aba0_results/1746537047984-8610-94b7aba0"
# PDB_DIR="/home/auri/extract_ppi/in/bce_var"
PDB_DIR="/home/auri/extract_ppi/in/EGFR_ECD/1746470428303-2758-312e7774_job_5XWD_24_affibodies_any_hotspot"
XML_SCRIPT="$(dirname "$0")/relax_and_metrics.xml"                  # XML protocol shipped alongside this script
OUTPUT_DIR="/home/auri/Rosetta/output/relax_metrics_$(date +%Y%m%d_%H%M%S)"  # timestamped output folder
NSTRUCT=1   # number of output models per input structure

# -----------------------------------------------------------------------------

# Check prerequisites ---------------------------------------------------------
[[ -x "${ROSETTA_BIN}" ]] || { echo "ERROR: Rosetta executable not found at ${ROSETTA_BIN}"; exit 1; }
[[ -d "${ROSETTA_DB}" ]] || { echo "ERROR: Rosetta database not found at ${ROSETTA_DB}"; exit 1; }
[[ -f "${XML_SCRIPT}" ]] || { echo "ERROR: XML protocol not found at ${XML_SCRIPT}"; exit 1; }
[[ -d "${PDB_DIR}" ]] || { echo "ERROR: PDB directory not found at ${PDB_DIR}"; exit 1; }

mkdir -p "${OUTPUT_DIR}"

# Main loop -------------------------------------------------------------------
find "${PDB_DIR}" -maxdepth 1 -type f \( -name "*.pdb" -o -name "*.pdb.gz" \) | while read -r pdb_file; do
    # Derive a clean base name (remove .pdb or .pdb.gz)
    base_name=$(basename "${pdb_file}")
    base_name=${base_name%%.pdb}
    base_name=${base_name%%.gz}

    echo "[INFO] Processing ${base_name}"

    current_output_dir="${OUTPUT_DIR}/${base_name}"
    mkdir -p "${current_output_dir}"

    # Construct Rosetta command ------------------------------------------------
    "$ROSETTA_BIN" \
        -database "$ROSETTA_DB" \
        -parser:protocol "$XML_SCRIPT" \
        -s "$pdb_file" \
        -nstruct "$NSTRUCT" \
        -out:path:all "$current_output_dir" \
        -out:prefix "${base_name}_" \
        -out:file:scorefile "${base_name}.sc" \
        -holes:dalphaball "$DALPHABALL_EXEC" \
        -run:preserve_header \
        -overwrite 2>&1 | tee "${current_output_dir}/${base_name}.log"

    echo "[DONE] ${base_name} → results in ${current_output_dir}"
    echo "--------------------------------------------------------------------"

done

echo "All structures processed. Combined output directory: ${OUTPUT_DIR}"