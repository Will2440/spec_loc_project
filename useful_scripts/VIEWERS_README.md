# SpecLoc Lazy Viewers

Three focused viewers for exploring spectral localiser simulation results.

## Overview

The viewer system is split into three specialized applications, each with a clear data model and purpose:

### 1. Ribbon Spectrum + Band3D Viewer (`qwz_ribbon_viewer.py`)

**Purpose**: Inspect per-gamma ribbon spectrum and 3D band structure

**Plots** (3×3 grid):
- Row 1: IPR-fixed, dC/dE-fixed, Chern-acc-fixed
- Row 2: IPR-auto, dC/dE-auto, Chern-acc-auto  
- Row 3: Band3D, (empty), (empty)

**Controls**:
- Sliders: A, B, m, gamma
- Toggles: Perturbation type, Range mode (fixed/auto)

**Launch**:
```bash
python3 useful_scripts/qwz_ribbon_viewer.py /path/to/lazy_packet
```

---

### 2. OBC Density of States Viewer (`qwz_obc_viewer.py`)

**Purpose**: Inspect DOS and LDOS for open boundary conditions

**Plots** (1×3 grid):
- Total DOS vs Energy
- LDOS at Target Energy
- LDOS at Lowest |E| States

**Controls**:
- Sliders: A, B, m, gamma, E (target energy for LDOS)
- Toggle: Perturbation type

**Auto-detection**: Automatically detects whether gamma-resolved or legacy records exist per plot type

**Launch**:
```bash
python3 useful_scripts/qwz_obc_viewer.py /path/to/lazy_packet
```

---

### 3. Spectral Localiser Viewer (`qwz_specloc_viewer.py`)

**Purpose**: Inspect aggregated spectral localiser analysis

**Core Plots** (always shown):
- Signature vs E
- log(gap) vs E
- Spectrum vs gamma
- Gap: gamma×E heatmap
- Signature: gamma×E heatmap

**Extended Plots** (toggle on/off):
- Gap/Signature: gamma×W, gamma×kappa, W×kappa
- Spectrum vs W, Spectrum vs kappa

**Controls**:
- Sliders: A, B, m, gamma, W, kappa, E
- Toggles: Perturbation type, Extended plots
- Spectrum y-limits: Custom y-axis range for spectrum plots

**Launch**:
```bash
python3 useful_scripts/qwz_specloc_viewer.py /path/to/lazy_packet
```

---

## Launch All Viewers

Start all three viewers simultaneously:

```bash
python3 useful_scripts/qwz_launch_all.py /path/to/lazy_packet
```

Press `Ctrl+C` to close all viewers.

---

## Workflow

### 1. Build Lazy Packet

```bash
julia simulations/data_processing/lazy_packet_builder.jl \
    "<input_root>" \
    "<output_root>" \
    "<run_id>"
```

This scans all `.jld2` files in `<input_root>` and creates:
- `lazy_packet_records.tsv` - Individual plot records
- `lazy_packet_groups.tsv` - Aggregated specloc groups
- `lazy_packet_stats.tsv` - Global statistics (e.g., dC/dE range)

### 2. Launch Viewers

```bash
python3 useful_scripts/qwz_launch_all.py "<output_root>/<run_id>/lazy_packet"
```

### 3. Navigate Data

Use sliders to navigate parameter space. Plots render on-demand and are cached locally in `lazy_cache/`.

---

## Features

**All Viewers**:
- ✅ Lazy on-demand rendering (no pre-rendering wait)
- ✅ Local disk cache with LRU management
- ✅ Per-plot metadata captions showing record ID and parameters
- ✅ Perturbation type toggle
- ✅ Debounced slider updates for smooth navigation

**Ribbon/OBC Viewers**:
- ✅ Direct case record matching (simple, fast)
- ✅ Auto-detects gamma-resolved vs legacy data

**Specloc Viewer**:
- ✅ Group record matching (aggregated across multiple cases)
- ✅ Extended plot toggle for large packets
- ✅ Custom spectrum y-limits
- ✅ E-reference marker on signature/loggap plots

---

## Architecture

**Data Pipeline**:
```
JLD2 files → lazy_packet_builder.jl → {records,groups,stats}.tsv
                                            ↓
                        Viewer → lazy_render.jl → cached PNGs
```

**Record Types**:
- `case`: Per-case plots (ribbon, band3d, dos, ldos)
- `specloc_group`: Aggregated specloc plots spanning multiple cases

**Matching Logic**:
- **Ribbon/OBC**: Direct lookup by (A, B, m, perturbation_type, gamma, ...)
- **Specloc**: Group lookup by (A, B, m, perturbation_type, ...) + parameter-specific axes

---

## Legacy Viewer

The old monolithic viewer is preserved as:
```
useful_scripts/qwz_lazy_plt_viewer_legacy.py
```

---

## Troubleshooting

**"No matching record"**: 
- Check that perturbation type matches data (toggle to cycle through available types)
- For OBC: Legacy packets may not have gamma-resolved records (caption shows `[legacy]` vs `[gamma-resolved]`)
- For Specloc: Ensure packet was built with updated `lazy_packet_builder.jl` that emits all plot types

**Blank/stale images**:
- Force cache rebuild: `SPECLOC_LAZY_FORCE_REBUILD=true` when launching viewer
- Or manually delete `lazy_cache/` directory

**Slow rendering**:
- Reduce prefetch radius: `SPECLOC_LAZY_PREFETCH_RADIUS=1`
- Or disable prefetch: `SPECLOC_LAZY_PREFETCH=0`

---

## Development Notes

**Line counts** (approx):
- Ribbon viewer: ~370 lines
- OBC viewer: ~390 lines  
- Specloc viewer: ~480 lines
- Launcher: ~60 lines
- **Total**: ~1300 lines (vs ~850 in legacy monolithic viewer)

**Complexity**: Each viewer is self-contained with simple, direct logic. No fallbacks, no virtual types, no gamma-aware branching across unrelated plot families.
