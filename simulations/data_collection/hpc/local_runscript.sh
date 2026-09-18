#!/bin/bash
# local_runscript.sh: Run all rows from a parameter .dat file in parallel on local machine.
#
# Usage:
#   ./local_runscript.sh                           # newest .dat, auto-detected cores
#   ./local_runscript.sh params_file.dat           # specific file
#   ./local_runscript.sh params_file.dat 4         # specific file, 4 parallel jobs

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

# Resolve params file
if [ -z "${1:-}" ]; then
    PARAMS_FILE=$(ls -t param_sets/*.dat 2>/dev/null | head -n 1)
    if [ -z "$PARAMS_FILE" ]; then
        echo "Error: No .dat files found in param_sets/. Run param_prep.jl first."
        exit 1
    fi
else
    PARAMS_FILE="$1"
fi

[ -f "$PARAMS_FILE" ] || { echo "Error: Parameter file not found: $PARAMS_FILE"; exit 1; }

NUM_JOBS="${2:-$(detect_cores)}"
NUM_ROWS=$(tail -n +2 "$PARAMS_FILE" | wc -l | tr -d ' ')

[ "$NUM_ROWS" -eq 0 ] && { echo "Error: No data rows in $PARAMS_FILE"; exit 1; }

RUN_ID="local_$(date +%Y%m%d_%H%M%S)_rows${NUM_ROWS}"
echo "Params:  $PARAMS_FILE"
echo "Rows:    $NUM_ROWS  |  Jobs: $NUM_JOBS  |  Run ID: $RUN_ID"
echo ""

PROGRESS_DIR=$(mktemp -d)
trap "rm -rf $PROGRESS_DIR" EXIT

show_progress() {
    local completed=$1 total=$2 start=$3
    local percent=$((completed * 100 / total))
    local bar_size=20
    local filled=$((percent * bar_size / 100))
    local empty=$((bar_size - filled))
    local bar="["
    for ((i=0; i<filled; i++)); do bar+="█"; done
    for ((i=0; i<empty; i++)); do bar+="░"; done
    bar+="]"
    local elapsed=$(($(date +%s) - start))
    local eta_str=""
    if [ $completed -gt 0 ]; then
        local rate=$((elapsed / completed))
        local eta_s=$(((total - completed) * rate))
        local em=$((eta_s / 60)); local es=$((eta_s % 60))
        [ $em -gt 0 ] && eta_str=$(printf "ETA: %dm%02ds" $em $es) || eta_str=$(printf "ETA: %ds" $es)
    fi
    printf "\r%-50s %3d%% [%d/%d] %s" "$bar" "$percent" "$completed" "$total" "$eta_str"
}

START_TIME=$(date +%s)

seq 1 "$NUM_ROWS" | xargs -P "$NUM_JOBS" -I {} \
    bash -c 'row="$1"
             SPECLOC_RUN_ID="'"$RUN_ID"'" julia --startup-file=no main.jl "$row" "'"$PARAMS_FILE"'" \
             >/dev/null 2>&1
             touch "'"$PROGRESS_DIR"'/done_${row}"' _ {}

# Final progress update
completed=$(ls -1 "$PROGRESS_DIR" 2>/dev/null | wc -l)
show_progress $completed $NUM_ROWS $START_TIME
echo ""

ELAPSED=$(($(date +%s) - START_TIME))
printf "Done: %d rows in %dm%02ds — results in results/%s/\n" \
    "$NUM_ROWS" "$((ELAPSED/60))" "$((ELAPSED%60))" "$RUN_ID"
