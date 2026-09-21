# Spectral Localiser κ-Scan Analysis

This directory contains scripts to analyze and visualize spectral localiser data from JLD2 result files, specifically for parameter scans over κ (kappa) values.

## Files

| File | Purpose |
|------|---------|
| `load_data.jl` | Core utility functions for loading JLD2 results and extracting data |
| `plot_kappa_scan.jl` | Main plotting routine that creates three plots from κ-scan data |
| `run_analysis.jl` | Convenience runner script (recommended entry point) |
| `README.md` | This file |

## Usage

### Quick Start

From this directory, run:

```bash
julia run_analysis.jl /path/to/result.jld2
```

For example:
```bash
julia run_analysis.jl ../../../data_collection/hpc/results/20260921_144508_r1/row1_chunk1_A1.0_B1.0_m-1.0_g0.0_ptsym_cos_diff_dtnone_d0.0_phi0.0_sx50_sy50.jld2
```

### Output Location & Organization

By default, plots are automatically organized in a parameter-based directory structure:

```
plots/
└── A{A}_B{B}_m{m}_g{gamma}_W{W}_pt{pert_type}_dt{disorder_type}_d{orbital_d}_phi{phi}/
    ├── 01_gap_vs_kappa_pt{pert_type}_A{A}_B{B}_m{m}_g{gamma}_W{W}_d{orbital_d}.png
    ├── 02_index_vs_kappa_pt{pert_type}_A{A}_B{B}_m{m}_g{gamma}_W{W}_d{orbital_d}.png
    └── 03_signature_vs_kappa_pt{pert_type}_A{A}_B{B}_m{m}_g{gamma}_W{W}_d{orbital_d}.png
```

The directory and filenames include all relevant system parameters for easy identification.

### What the Script Does

1. **Load** the JLD2 file and extract gap, index, and signature data
2. **Average** over the energy dimension for each κ value
3. **Create 3 plots** with clear filenames:
   - `01_gap_vs_kappa_...` — Spectral localiser gap vs κ
   - `02_index_vs_kappa_...` — Spectral localiser invariant vs κ  
   - `03_signature_vs_kappa_...` — Spectral localiser signature vs κ
4. **Print summary statistics** to the terminal

### Custom Output Directory

To override the automatic organization and save to a custom location:

```bash
julia run_analysis.jl /path/to/result.jld2 /custom/output/dir
```

When a custom directory is provided, plots still use parameter-based filenames but are saved directly to that path (no automatic subdirectory created).

## Module Functions

### `load_data.jl`

#### `load_specloc_result(filepath::String)`
Load a JLD2 result file and return a NamedTuple containing:
- `metadata`: Dict of simulation parameters (A, B, m, γ, W, etc.)
- `axes`: Dict with vectors (gammas, Ws, kappas, Es)
- `gap`: 4D array (ng, nW, nk, nE)
- `signature`: 4D array (ng, nW, nk, nE)
- `index`: 4D array (ng, nW, nk, nE) [= 0.5 × signature]

#### `extract_gap_vs_kappa(result::NamedTuple)`
Extract gap vs κ for single γ, single W case.
Returns: `(kappas::Vector, gap_vs_kappa::Vector)`

#### `extract_index_vs_kappa(result::NamedTuple)`
Extract spectral localiser invariant vs κ for single γ, single W case.
Returns: `(kappas::Vector, index_vs_kappa::Vector)`

#### `extract_signature_vs_kappa(result::NamedTuple)`
Extract spectral localiser signature vs κ for single γ, single W case.
Returns: `(kappas::Vector, signature_vs_kappa::Vector)`

### `plot_kappa_scan.jl`

#### `plot_kappa_scan(result_file::String; output_dir::String="")`
Main entry point. Loads result file, extracts κ-scan data, creates plots, and prints summary statistics.

- If `output_dir=""` (default): plots are saved to `./plots/<param_summary>/` where `<param_summary>` encodes all system parameters
- If `output_dir` is provided: plots are saved directly to that directory (with parameter-based filenames)

## Example Output

Running the analysis on your κ-scan data will produce output like:

```
📂 Loading: /path/to/result.jld2
✓ Saved: /Users/Will/Documents/spec_loc_project/simulations/data_processing/kappa_scan/plots/A1.00_B1.00_m-1.00_g0.00_W0.00_ptsym_cos_diff_dtnone_d0.00_phi0.00/01_gap_vs_kappa_ptsym_cos_diff_A1.00_B1.00_m-1.00_g0.00_W0.00_d0.00.png
✓ Saved: /Users/Will/Documents/spec_loc_project/simulations/data_processing/kappa_scan/plots/A1.00_B1.00_m-1.00_g0.00_W0.00_ptsym_cos_diff_dtnone_d0.00_phi0.00/02_index_vs_kappa_ptsym_cos_diff_A1.00_B1.00_m-1.00_g0.00_W0.00_d0.00.png
✓ Saved: /Users/Will/Documents/spec_loc_project/simulations/data_processing/kappa_scan/plots/A1.00_B1.00_m-1.00_g0.00_W0.00_ptsym_cos_diff_dtnone_d0.00_phi0.00/03_signature_vs_kappa_ptsym_cos_diff_A1.00_B1.00_m-1.00_g0.00_W0.00_d0.00.png

============================================================
Spectral Localiser κ-scan Summary
============================================================
System parameters:
  A=1.0000  B=1.0000  m=-1.0000
  γ=0.0000  W=0.0000
  Perturbation type: sym_cos_diff
  Disorder type: none
  Orbital displacement: 0.0
  Orbital angle φ: 0.0

Data shapes:
  Number of κ values: 50
  κ range: [0.000001, 0.100000]

Invariant statistics:
  Min index:  0.000000
  Max index:  1.000000
  Mean index: 0.580000

Gap statistics (E-averaged):
  Min gap:  1.003656e-04
  Max gap:  9.785678e-01
  Mean gap: 2.251735e-01

Plots saved to: /Users/Will/Documents/spec_loc_project/simulations/data_processing/kappa_scan/plots/A1.00_B1.00_m-1.00_g0.00_W0.00_ptsym_cos_diff_dtnone_d0.00_phi0.00
============================================================
```

## Data Format Notes

The JLD2 files saved by `solver.jl` have the following structure:

- **Dimensions**: (n_gammas, n_Ws, n_kappas, n_Es)
- For κ-scan data, typically: n_gammas=1, n_Ws=1, n_kappas=variable, n_Es=variable
- The `gap`, `signature`, and `index` arrays follow this structure
- All energy-dependent quantities are averaged over the E dimension

## Requirements

- Julia 1.6+
- JLD2.jl (for loading)
- Plots.jl + backend (e.g., GR, PyPlot, PlotlyJS)
- Statistics (stdlib)

## Installation

If needed, install dependencies from Julia REPL:

```julia
using Pkg
Pkg.add(["JLD2", "Plots"])
```

## Notes

- The analysis assumes a **single γ and single W** case (as typical for κ-scan data)
- Energy averaging is performed via `mean(data; dims=last_dimension)`
- Plot titles and directory/filenames are auto-generated from metadata
- **Directory structure** encodes all system parameters, enabling easy navigation across different parameter sets
- Filenames are prefixed with `01_`, `02_`, `03_` for convenient sorting and identification
- All physical units are as defined in the simulation (`solver.jl`)
