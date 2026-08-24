# SpecLoc Viewer Refactor - Migration Guide

## What Changed (2026-08-17)

Split the monolithic viewer into **three focused viewers** + a launcher.

### Before (Monolithic)

**Single viewer** (`qwz_lazy_plt_viewer.py`):
- Mixed 3 incompatible data models (case per-gamma, case legacy, group aggregated)
- ~850 lines with complex fallback logic
- Virtual plot types for missing heatmaps
- Gamma-aware branching across all OBC plots
- Selection cache mixing unrelated plot families
- "No matching record" errors from lookup mismatches

### After (Focused Split)

**Three specialized viewers** (~1300 total lines, simpler per-viewer):

1. **Ribbon Viewer** (`qwz_ribbon_viewer.py`, ~370 lines)
   - Only ribbon spectrum + band3D
   - Direct case lookup: (A, B, m, gamma, perturbation)
   - No fallbacks, no virtual types

2. **OBC Viewer** (`qwz_obc_viewer.py`, ~390 lines)
   - Only DOS/LDOS plots
   - Auto-detects gamma vs legacy per plot type independently
   - Clear caption showing `[gamma-resolved]` vs `[legacy]`

3. **Specloc Viewer** (`qwz_specloc_viewer.py`, ~480 lines)
   - Only spectral localiser plots
   - Group record lookup (aggregated data)
   - All heatmaps are real records (no virtual fallback)
   - Extended plots toggle for large packets

**Launcher** (`qwz_launch_all.py`, ~60 lines):
- Opens all three viewers simultaneously
- Single Ctrl+C closes all

---

## Benefits

| Aspect | Before | After |
|--------|--------|-------|
| **Complexity** | High (one viewer handles everything) | Low (each viewer single-purpose) |
| **Matching logic** | Complex fallbacks + virtual types | Direct record lookup |
| **Debugging** | Hard (pointing errors across types) | Easy (isolated scope) |
| **Cache** | One cache for all types | Per-viewer cache (independent) |
| **Gamma handling** | All-or-nothing OBC gating | Auto-detect per plot type |
| **Heatmaps** | Virtual fallback | Real records emitted by builder |
| **Line count** | ~850 lines | ~370 + ~390 + ~480 = ~1240 lines (but simpler) |
| **Selection cache** | One cache mixing unrelated types | Per-viewer, plot-relevant keys only |

---

## Quick Start

### 1. Rebuild Your Packet (Recommended)

This ensures all plot types (including heatmaps) are explicit records:

```bash
julia simulations/data_processing/lazy_packet_builder.jl \
    "<input_root>" \
    "<output_root>" \
    "20260817_refactor"
```

### 2. Launch All Viewers

```bash
python3 useful_scripts/qwz_launch_all.py \
    "<output_root>/20260817_refactor/lazy_packet"
```

### 3. Or Launch Individually

```bash
# Ribbon + Band3D only
python3 useful_scripts/qwz_ribbon_viewer.py /path/to/lazy_packet

# DOS/LDOS only
python3 useful_scripts/qwz_obc_viewer.py /path/to/lazy_packet

# Spectral localiser only
python3 useful_scripts/qwz_specloc_viewer.py /path/to/lazy_packet
```

---

## Backward Compatibility

### Old Packets (Pre-Refactor)

- **Ribbon/OBC viewers**: Will work with old packets
- **Specloc viewer**: May show "no matching record" for extended heatmaps if packet was built with old builder (missing gamma×E, gamma×W, etc.)
  - **Fix**: Rebuild packet with new `lazy_packet_builder.jl`

### Legacy Viewer

Preserved as `useful_scripts/qwz_lazy_plt_viewer_legacy.py` for reference. Not recommended for new work.

---

## Code Changes Summary

### New Files
- `useful_scripts/qwz_ribbon_viewer.py` (new)
- `useful_scripts/qwz_obc_viewer.py` (new)
- `useful_scripts/qwz_specloc_viewer.py` (new)
- `useful_scripts/qwz_launch_all.py` (new)
- `useful_scripts/VIEWERS_README.md` (new)
- `useful_scripts/VIEWERS_MIGRATION.md` (this file)

### Renamed Files
- `qwz_lazy_plt_viewer.py` → `qwz_lazy_plt_viewer_legacy.py`

### Modified Files
- `simulations/data_processing/lazy_packet_builder.jl`
  - Now emits all 13 specloc plot types explicitly (no more missing heatmaps)
  - Independent gamma detection for dos/ldos_target/ldos_lowest

- `simulations/data_processing/plotting.jl`
  - Independent gamma handling for DOS/LDOS (was all-or-nothing)

---

## Testing Checklist

After rebuilding your packet and launching viewers:

### Ribbon Viewer
- [ ] Ribbon spectrum IPR shows edge states
- [ ] dC/dE colorbar range matches packet stats
- [ ] Chern accumulation crosses zero correctly
- [ ] Band3D shows gap at expected k-points
- [ ] Perturbation toggle cycles through available types
- [ ] Range toggle switches between fixed/auto limits

### OBC Viewer
- [ ] DOS shows gap/peak structure
- [ ] LDOS at target energy shows localized states
- [ ] LDOS at lowest |E| shows bulk vs edge
- [ ] Caption shows `[gamma-resolved]` for new data or `[legacy]` for old data
- [ ] Gamma slider changes which gamma slice is shown (if gamma-resolved)

### Specloc Viewer
- [ ] Signature vs E shows transitions (steps)
- [ ] log(gap) vs E shows minima at phase boundaries
- [ ] Spectrum vs gamma shows connected eigenvalue evolution
- [ ] Gap gamma×E heatmap shows phase diagram
- [ ] Signature gamma×E heatmap shows topological index
- [ ] Extended plots toggle adds/removes W/kappa families
- [ ] Spectrum y-limits input constrains vertical range

### All Viewers
- [ ] Per-plot metadata caption shows record ID and parameters
- [ ] Images render on-demand (not blank)
- [ ] Slider navigation is smooth and responsive
- [ ] Cache grows in `lazy_cache/` directory

---

## Troubleshooting

**"No matching record" in Specloc Viewer for heatmaps:**
- Your packet was built with old builder that only emitted 5 plot types
- **Fix**: Rebuild packet with updated `lazy_packet_builder.jl`

**OBC shows `[legacy]` but you expect gamma-resolved:**
- Your simulation data lacks `dos_by_gamma` / `ldos_target_by_gamma` / `ldos_lowest_by_gamma` payloads
- **Fix**: Re-run simulations with updated `solver.jl` that computes gamma-resolved OBC

**Viewers crash on launch:**
- Check PyQt5 is installed: `pip install PyQt5`
- Check packet path exists and contains `lazy_packet_records.tsv`

---

## Performance Notes

**Compared to legacy monolithic viewer:**
- ✅ Faster record matching (smaller candidate sets per viewer)
- ✅ Simpler cache management (no mixing of unrelated types)
- ✅ Lower memory usage per viewer (each loads only relevant records)
- ✅ Easier debugging (isolated failure modes)

**Recommended for large packets:**
- Launch only the viewer(s) you need (not all three)
- Disable prefetch: `SPECLOC_LAZY_PREFETCH=0`
- Keep extended plots off in specloc viewer (default for >5000 records)

---

## Further Reading

See `useful_scripts/VIEWERS_README.md` for detailed usage documentation.
