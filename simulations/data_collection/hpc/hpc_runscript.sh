#!/bin/bash
#SBATCH --job-name=SpecLoc_QWZ_noW
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=compute
#SBATCH --cpus-per-task=1
#SBATCH --account=phys030424
#SBATCH --array=1-200
#SBATCH --time=05:00:00
#SBATCH --mem=1G

set -euo pipefail

SCRIPT_DIR="/user/work/hb21877/SpecLoc_project/spec_loc_project/simulations/data_collection/hpc"
cd "$SCRIPT_DIR"

ROW_INDEX="${SLURM_ARRAY_TASK_ID:-${1:-}}"
PARAMS_FILE="${2:-${1:-}}"

if [ -z "${ROW_INDEX:-}" ]; then
    echo "Usage: $0 <row_index> [params_file]" >&2
    exit 1
fi

if [ -z "${PARAMS_FILE:-}" ]; then
    PARAMS_FILE="$(ls -t "$SCRIPT_DIR/param_sets"/*.dat 2>/dev/null | head -n 1)"
fi

if [ -z "${PARAMS_FILE:-}" ] || [ ! -f "${PARAMS_FILE}" ]; then
    echo "Error: Parameter file not found: ${PARAMS_FILE:-<none>}" >&2
    exit 1
fi

module add languages/julia || true
export JULIA_DEPOT_PATH="/user/work/hb21877/.julia"
export JULIA_NUM_THREADS=1

julia --startup-file=no "$SCRIPT_DIR/main.jl" "$ROW_INDEX" "$PARAMS_FILE"