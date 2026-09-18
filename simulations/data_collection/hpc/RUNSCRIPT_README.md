# HPC Batch Runner Scripts

## Overview

This folder contains the updated HPC-compatible spectral localiser pipeline for the QWZ model.
Key improvements over `hpc_legacy/`:
- **New perturbation types**: `:sym_cos_sum/diff/add/sub`, `:asym_sin_sum/diff/add/sub`, `:tilt` (legacy)
- **Orbital sublattice embedding**: `orbital_displacement` and `phi` parameters shift the position operators used in the spectral localiser
- **Disorder averaging**: `n_disorder_realisations` averages spectral localiser over independent disorder realisations
- **Kappa scaling**: `scale_kappa_to_L` sets κ = scale / L_x for finite-size robustness
- **Fast sparse solver**: LDLt (Sylvester's law) for signature; KrylovKit shift-and-invert for gap/spectrum slices
- **COO Hamiltonian assembly**: sparse build with correct diagonal hopping terms for `_add/_sub` types

---

## Files

| File | Purpose |
|------|---------|
| `param_prep.jl` | Generates `.dat` parameter files with full new parameter set |
| `main.jl` | Runs one row: `julia main.jl <row_index> [params_file]` |
| `solver.jl` | `SpecLocSolver` module — all physics computations |
| `local_runscript.sh` | Parallel local testing (all rows via xargs) |
| `hpc_runscript.sh` | SLURM job array runner (one row per task) |

---

## New Parameters (columns 31–35 in .dat)

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `orbital_displacement` | Float64 | `0.0` | Sublattice position offset magnitude `d` for X, Y operators |
| `phi` | Float64 | `0.0` | Sublattice displacement angle φ (radians) |
| `n_disorder_realisations` | Int | `1` | Number of disorder samples to average over (ignored when W=0) |
| `scale_kappa_to_L` | Bool | `false` | If true, κ = `kappa_scales[i] / Lx_obc` |
| `kappa_scales` | Vector{Float64} | `[0.0004]` | Scale factors when `scale_kappa_to_L=true` |

---

## Perturbation Types

| Symbol | k-space scalar term |
|--------|-------------------|
| `:none` | 0 |
| `:sym_cos_sum` | γ(cos kx + cos ky) |
| `:sym_cos_diff` | γ(cos kx − cos ky) |
| `:sym_cos_add` | γ cos(kx + ky) — uses diagonal hoppings |
| `:sym_cos_sub` | γ cos(kx − ky) — uses diagonal hoppings |
| `:asym_sin_sum` | γ(sin kx + sin ky) |
| `:asym_sin_diff` | γ(sin kx − sin ky) |
| `:asym_sin_add` | γ sin(kx + ky) — uses diagonal hoppings |
| `:asym_sin_sub` | γ sin(kx − ky) — uses diagonal hoppings |
| `:tilt` | −γ sin kx (legacy) |
| `:symmetric` | alias for `:sym_cos_sum` |

---

## Quick Start

### Local testing
```bash
cd simulations/data_collection/hpc

# 1. Generate parameter file
julia --startup-file=no param_prep.jl

# 2. Run all rows in parallel (auto-detects cores)
./local_runscript.sh

# 3. Or specify file and parallelism
./local_runscript.sh param_sets/params_*.dat 4
```

### HPC (SLURM)
```bash
# Edit SCRIPT_DIR and account in hpc_runscript.sh, then:
sbatch --array=1-<N_rows> hpc_runscript.sh param_sets/params_*.dat
```

### Process results
```bash
cd simulations/data_processing
julia main.jl   # uses newest results/ dir by default
```

---

## Julia Dependencies

`solver.jl` requires **KrylovKit** in addition to stdlib packages.
Install once in the target Julia environment:
```julia
import Pkg; Pkg.add("KrylovKit")
```
All other dependencies (`LinearAlgebra`, `SparseArrays`, `Statistics`, `Random`) are stdlib.

---

## Legacy Reference

The original scripts (without new perturbation types or sparse solver) are preserved unchanged in `hpc_legacy/`.
