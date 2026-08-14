#!/bin/bash
# runscript.sh: Run a single row for HPC job array submission
#
# Usage (SLURM example):
#   sbatch --array=1-N runscript.sh params_file.dat
#   export SPECLOC_RUN_ID=my_batch_20260814   # optional, shared output folder
#
# Usage (manual):
#   ./runscript.sh <row_index> [params_file]
#
# Examples:
#   ./runscript.sh 1 params_file.dat
#   ./runscript.sh 5

set -e

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$SCRIPT_DIR"

# Determine row index
ROW_INDEX="${SLURM_ARRAY_TASK_ID:-$1}"
if [ -z "$ROW_INDEX" ]; then
    echo "Usage: $0 <row_index> [params_file]"
    exit 1
fi

# Determine params file
if [ -n "$SLURM_ARRAY_TASK_ID" ]; then
    PARAMS_ARG="$1"
else
    PARAMS_ARG="$2"
fi

if [ -z "$PARAMS_ARG" ]; then
    # Use most recent .dat file
    PARAMS_FILE=$(ls -t param_sets/*.dat 2>/dev/null | head -n 1)
    if [ -z "$PARAMS_FILE" ]; then
        echo "Error: No .dat files found in param_sets/"
        exit 1
    fi
else
    PARAMS_FILE="$PARAMS_ARG"
fi

if [ ! -f "$PARAMS_FILE" ]; then
    echo "Error: Parameter file not found: $PARAMS_FILE"
    exit 1
fi

echo "=========================================="
echo "HPC Job: Row $ROW_INDEX"
echo "=========================================="
echo "Params file: $PARAMS_FILE"
echo "Row index:   $ROW_INDEX"
echo "=========================================="
echo ""

# Run the single row
julia --startup-file=no main.jl "$ROW_INDEX" "$PARAMS_FILE"

echo ""
echo "=========================================="
echo "Row $ROW_INDEX completed!"
echo "=========================================="
