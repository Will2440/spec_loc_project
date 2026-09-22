#!/bin/bash
#SBATCH --job-name=SpecLoc_symcosdiff_gsampled
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=general
#SBATCH --qos=regular
#SBATCH --cpus-per-task=4
#SBATCH --array=1-1
#SBATCH --time=24:00:00
#SBATCH --mem=24G

set -euo pipefail

ROW_INDEX="${SLURM_ARRAY_TASK_ID:-${1:-}}"
PARAMS_FILE="${2:-${1:-}}"

if [ -z "${ROW_INDEX:-}" ]; then
    echo "Usage: $0 <row_index> [params_file]" >&2
    exit 1
fi

# Get the actual hpc directory from PARAMS_FILE path
PARAMS_FILE_ABS="$(cd "$(dirname "$PARAMS_FILE")" 2>/dev/null && pwd)/$(basename "$PARAMS_FILE")" || true
if [ -z "$PARAMS_FILE_ABS" ]; then
    SCRIPT_PATH="$(readlink -f "${BASH_SOURCE[0]}")"
    SCRIPT_DIR="$(dirname "$SCRIPT_PATH")"
    PARENT_DIR="$(dirname "$SCRIPT_DIR")"
    PARAMS_FILE_ABS="$(ls -t "$PARENT_DIR/param_sets"/*.dat 2>/dev/null | head -n 1)"
else
    PARAMS_FILE="$PARAMS_FILE_ABS"
fi

if [ -z "${PARAMS_FILE_ABS:-}" ] || [ ! -f "${PARAMS_FILE_ABS}" ]; then
    echo "Error: Parameter file not found: ${PARAMS_FILE_ABS:-<none>}" >&2
    exit 1
fi

PARAMS_FILE="$PARAMS_FILE_ABS"
PARENT_DIR="$(dirname "$(dirname "$PARAMS_FILE")")"

# Get script directory and cd there
SCRIPT_PATH="$(readlink -f "${BASH_SOURCE[0]}")"
SCRIPT_DIR="$(dirname "$SCRIPT_PATH")"

cd "$SCRIPT_DIR"

# Load Julia and OpenBLAS modules for BLAS threading support
module load Julia/1.11.6-linux-x86_64 || module load Julia || true
module load OpenBLAS/0.3.29-GCC-14.2.0 || true

# Set Julia depot path based on cluster
if [[ "$HOSTNAME" == *"hyperion"* ]]; then
    export JULIA_DEPOT_PATH="$HOME/.julia"
elif [[ "$HOSTNAME" == *"bp"* ]] || [[ "$HOSTNAME" == *"bluecrystal"* ]]; then
    export JULIA_DEPOT_PATH="/user/work/hb21877/.julia"
else
    # Default fallback
    export JULIA_DEPOT_PATH="$HOME/.julia"
fi

# Configure multithreading for both Julia and BLAS
# Set to 4 threads (benchmark-optimal on Hyperion)
export JULIA_NUM_THREADS=4
export OPENBLAS_NUM_THREADS=4
export OPENBLAS_VERBOSE=0  # Suppress BLAS startup messages

# Auto-install required packages if missing
julia --startup-file=no -e 'import Pkg; Pkg.add(["JLD2", "KrylovKit"])' 2>/dev/null

LOGS_DIR="$PARENT_DIR/logs"
mkdir -p "$LOGS_DIR"

LOG_FILE="$LOGS_DIR/job_${SLURM_ARRAY_JOB_ID:-noID}_row_${ROW_INDEX}_$(date +%Y%m%d_%H%M%S).log"

julia --startup-file=no "$PARENT_DIR/main.jl" "$ROW_INDEX" "$PARAMS_FILE" > "$LOG_FILE" 2>&1
