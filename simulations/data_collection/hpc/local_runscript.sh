#!/bin/bash
# local_runscript.sh: Run all rows from a parameter .dat file in parallel on local machine
#
# Usage:
#   ./local_runscript.sh [params_file] [num_cores]
#
# Examples:
#   ./local_runscript.sh                  # Uses newest .dat and auto-detected cores
#   ./local_runscript.sh params_file.dat  # Specific .dat file
#   ./local_runscript.sh params_file.dat 4 # Specific .dat file, 4 parallel jobs

set -e

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$SCRIPT_DIR"

detect_cores() {
    if command -v nproc >/dev/null 2>&1; then
        nproc
    elif command -v sysctl >/dev/null 2>&1; then
        sysctl -n hw.ncpu
    else
        echo 4
    fi
}

# Determine params file
if [ -z "$1" ]; then
    # Use most recent .dat file
    PARAMS_FILE=$(ls -t param_sets/*.dat 2>/dev/null | head -n 1)
    if [ -z "$PARAMS_FILE" ]; then
        echo "Error: No .dat files found in param_sets/"
        exit 1
    fi
else
    PARAMS_FILE="$1"
fi

if [ ! -f "$PARAMS_FILE" ]; then
    echo "Error: Parameter file not found: $PARAMS_FILE"
    exit 1
fi

# Determine number of parallel jobs
if [ -z "$2" ]; then
    NUM_JOBS=$(detect_cores)
else
    NUM_JOBS="$2"
fi

# Count rows in .dat file (header + data rows)
NUM_ROWS=$(tail -n +2 "$PARAMS_FILE" | wc -l | tr -d ' ')

RUN_ID="local_$(date +%Y%m%d_%H%M%S)_rows${NUM_ROWS}"
RESULTS_DIR="results/${RUN_ID}"

echo "=========================================="
echo "Local Parallel Batch Runner"
echo "=========================================="
echo "Params file: $PARAMS_FILE"
echo "Total rows:  $NUM_ROWS"
echo "Parallel jobs: $NUM_JOBS"
echo "Run ID: $RUN_ID"
echo "Output dir: $RESULTS_DIR"
echo "=========================================="
echo ""

if [ "$NUM_ROWS" -eq 0 ]; then
    echo "Error: No data rows found in $PARAMS_FILE"
    exit 1
fi

# Check for xargs or fall back to sequential
if command -v xargs >/dev/null 2>&1; then
    echo "Using xargs for job distribution"
    seq 1 "$NUM_ROWS" | xargs -P "$NUM_JOBS" -I {} \
        bash -c 'row="$1"; echo "[Row ${row}/'"$NUM_ROWS"'] Starting..."; SPECLOC_RUN_ID="'"$RUN_ID"'" julia --startup-file=no main.jl "$row" "'"$PARAMS_FILE"'" 2>&1 | sed "s/^/[Row ${row}] /"' _ {}
else
    echo "Warning: xargs not found, running sequentially (slow)"
    for row in $(seq 1 "$NUM_ROWS"); do
        echo "[Row $row/$NUM_ROWS] Starting..."
        SPECLOC_RUN_ID="$RUN_ID" julia --startup-file=no main.jl "$row" "$PARAMS_FILE"
    done
fi

echo ""
echo "=========================================="
echo "All rows completed!"
echo "=========================================="
