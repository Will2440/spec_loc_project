# HPC Batch Runner Scripts

## local_runscript.sh (Local Testing)

Run **all rows** from a parameter .dat file in **parallel** on your laptop/local machine.

### Usage

```bash
cd simulations/data_collection/hpc

# Auto-detect newest .dat file and number of cores
./local_runscript.sh

# Specific .dat file, auto-detect cores
./local_runscript.sh param_sets/params_20260814_152612_explicit_Eauto_rows16.dat

# Specific .dat file and number of parallel jobs
./local_runscript.sh param_sets/params_20260814_152612_explicit_Eauto_rows16.dat 4
```

### Features

- **Parallel execution**: Automatically detects CPU cores and runs multiple rows simultaneously
- **Auto job distribution**: Uses xargs `-P` for parallel row dispatch
- **Fallback support**: Falls back to sequential if parallel tools unavailable
- **Output tagging**: Each row's output is prefixed with `[Row N]` for clarity
- **Progress tracking**: Prints total rows and parallel job count before starting
- **Shared run folder**: One local launch writes all row outputs into a single `results/<run_id>/` folder

### Example

```bash
# Run all 16 rows with 4 parallel jobs
./local_runscript.sh param_sets/params_20260814_152612_explicit_Eauto_rows16.dat 4
```

Output shows rows running concurrently, with row/chunk-tagged files written into one shared run directory.

---

## runscript.sh (HPC Job Array)

Run a **single row** for use in HPC job arrays (SLURM, PBS, etc.)

### Usage

```bash
# Direct invocation
cd simulations/data_collection/hpc
./runscript.sh 1 params_file.dat
./runscript.sh 5

# SLURM job array
sbatch --array=1-16 runscript.sh params_file.dat
```

### Features

- **Single-row execution**: Each invocation runs exactly one row
- **SLURM integration**: Automatically detects `SLURM_ARRAY_TASK_ID` for job arrays
- **Auto params detection**: Falls back to newest .dat file if not specified
- **Lightweight**: No threading overhead; one core per job

### Example SLURM Script

Create `submit_batch.sh`:

```bash
#!/bin/bash
#SBATCH --job-name=specloc_batch
#SBATCH --array=1-16
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=results/slurm-%A_%a.log

cd simulations/data_collection/hpc
./runscript.sh params_file.dat
```

Submit with:
```bash
sbatch submit_batch.sh
```

---

## Key Differences

| Aspect | local_runscript.sh | runscript.sh |
|--------|-------------------|--------------|
| **Execution** | All rows in parallel | Single row per call |
| **Use case** | Local testing on laptop | HPC batch submission |
| **Parallelism** | Uses all available cores | One core per job |
| **Output** | Interleaved with row tags | Sequential, single result |
| **Job control** | Shell script | SLURM/PBS array-friendly |

---

## Workflow

1. **Local testing**:
   ```bash
   ./local_runscript.sh param_sets/params_small_test.dat
   ```

2. **After validation, submit to HPC**:
   ```bash
   sbatch --array=1-256 runscript.sh param_sets/params_large_run.dat
   ```

3. **Process results**:
   ```bash
   cd ../../../simulations/data_processing
   julia main.jl  # Uses newest results dir by default
   ```

4. **Visualize**:
   ```bash
   python3 useful_scripts/qwz_plt_viewer.py
   ```
