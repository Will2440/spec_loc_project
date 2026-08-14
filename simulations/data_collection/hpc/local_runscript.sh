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

echo "Output: $RESULTS_DIR"
echo ""

if [ "$NUM_ROWS" -eq 0 ]; then
    echo "Error: No data rows found in $PARAMS_FILE"
    exit 1
fi

# Create temp directory for progress tracking
PROGRESS_DIR=$(mktemp -d)
trap "rm -rf $PROGRESS_DIR" EXIT

# Progress bar function
show_progress() {
    local completed=$1
    local total=$2
    local start_time=$3
    
    local percent=$((completed * 100 / total))
    local bar_size=20
    local filled=$((percent * bar_size / 100))
    local empty=$((bar_size - filled))
    
    # Build progress bar
    local bar="["
    for ((i=0; i<filled; i++)); do bar+="█"; done
    for ((i=0; i<empty; i++)); do bar+="░"; done
    bar+="]"
    
    # Calculate ETA
    local elapsed=$(($(date +%s) - start_time))
    local rate=0
    if [ $completed -gt 0 ]; then
        rate=$((elapsed / completed))
    fi
    
    local remaining=$((total - completed))
    local eta_seconds=$((remaining * rate))
    
    local eta_str=""
    if [ $eta_seconds -gt 0 ]; then
        local eta_mins=$((eta_seconds / 60))
        local eta_secs=$((eta_seconds % 60))
        if [ $eta_mins -gt 0 ]; then
            eta_str=$(printf "ETA: %dm%02ds" $eta_mins $eta_secs)
        else
            eta_str=$(printf "ETA: %ds" $eta_secs)
        fi
    fi
    
    printf "\r%-50s %3d%% [%d/%d] %s" "$bar" "$percent" "$completed" "$total" "$eta_str"
}

# Start time
START_TIME=$(date +%s)

# Run jobs with progress tracking
if command -v xargs >/dev/null 2>&1; then
    seq 1 "$NUM_ROWS" | xargs -P "$NUM_JOBS" -I {} \
        bash -c 'row="$1"; SPECLOC_RUN_ID="'"$RUN_ID"'" julia --startup-file=no main.jl "$row" "'"$PARAMS_FILE"'" >/dev/null 2>&1; touch "'"$PROGRESS_DIR"'/done_${row}"' _ {}
    
    # Monitor progress
    completed=0
    while [ $completed -lt $NUM_ROWS ]; do
        completed=$(ls -1 "$PROGRESS_DIR" 2>/dev/null | wc -l)
        show_progress $completed $NUM_ROWS $START_TIME
        if [ $completed -lt $NUM_ROWS ]; then
            sleep 0.5
        fi
    done
    show_progress $NUM_ROWS $NUM_ROWS $START_TIME
    echo ""
else
    echo "Warning: xargs not found, running sequentially (slow)"
    for row in $(seq 1 "$NUM_ROWS"); do
        SPECLOC_RUN_ID="$RUN_ID" julia --startup-file=no main.jl "$row" "$PARAMS_FILE"
        show_progress $row $NUM_ROWS $START_TIME
    done
    echo ""
fi
