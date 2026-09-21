# Julia Package Setup for HPC Pipeline

## Quick Setup (Do This First!)

Run this once per login session (packages may be wiped between logins):

```bash
cd /scratch/caiger/spec_loc_project/simulations/data_collection/hpc
module load Julia
julia setup_julia_packages.jl
```

This takes ~30-60 seconds and installs:
- **JLD2** - HDF5-based data file format for results
- **KrylovKit** - Krylov subspace methods for eigenvalue solving

## Required Packages

| Package | Purpose | Required For |
|---------|---------|--------------|
| JLD2 | Save/load results | main.jl (results output) |
| KrylovKit | Eigenvalue solver | solver.jl (spectral computation) |
| LinearAlgebra | Standard library | Built-in |
| SparseArrays | Standard library | Built-in |
| Random | Standard library | Built-in |
| Dates | Standard library | Built-in |
| DelimitedFiles | Standard library | Built-in |
| Printf | Standard library | Built-in |

## Manual Installation (if script fails)

```bash
julia -e 'import Pkg; Pkg.add("JLD2")'
julia -e 'import Pkg; Pkg.add("KrylovKit")'
```

## Verification

After setup, test that packages load:

```bash
julia -e 'using JLD2; using KrylovKit; println("✓ All packages loaded successfully!")'
```

## Cluster-Specific Notes

- **Hyperion**: Julia depot stored in `$SCRATCH/.julia`
- **Blue Pebble**: Julia depot stored in `/user/work/$USER/.julia`

The runscripts auto-detect the cluster and set appropriate paths.

## Troubleshooting

**If packages disappear on next login:**
This is normal on HPC systems. Re-run `julia setup_julia_packages.jl` at the start of your session.

**If you need persistent installations:**
Ask the HPC admins to install packages system-wide, or use a Julia environment file (Project.toml/Manifest.toml) - contact if you need this set up.
