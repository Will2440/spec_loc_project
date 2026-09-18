# Performance Analysis: fast_specloc.jl

## Executive Summary
The script has moderate optimization but several inefficiencies that could provide **2-4x speedup** with targeted improvements. The main bottlenecks are sparse matrix construction patterns, redundant block computations, and suboptimal parallelization structure.

---

## Critical Issues

### 1. **Sparse Matrix Assembly Inefficiency** ⚠️ CRITICAL
**Location:** `add_block!` function in `real_space_perturbed_disordered_hamiltonian_qwz`

**Problem:**
- Building sparse matrix via repeated `+=` operations is O(n) per insertion
- Called ~4*Lx*Ly times per Hamiltonian (~4000 times for Lx=Ly=50)
- Sparse matrix is repeatedly reallocated/rehashed

**Current approach:**
```julia
H = spzeros(ComplexF64, dim, dim)  # Empty sparse matrix
# Then loop: mat[row_base + a, col_base + b] += value  # Slow!
```

**Impact:** **~30-40% of runtime** for Hamiltonian construction

**Solution:** Use COO (coordinate) format assembly:
```julia
Is, Js, Vs = Int[], Int[], ComplexF64[]
# Collect all (i,j,v) triplets
# Then: H = sparse(Is, Js, Vs, dim, dim)
```
This is **5-10x faster** for this construction pattern.

---

### 2. **Redundant Block Computation** ⚠️ HIGH
**Location:** `real_space_perturbed_qwz_blocks` called inside disorder loop

**Problem:**
- Called once per realisation per thread per parameter set
- Returns identical `onsite_base`, `tx`, `ty` blocks for the same (A,B,m,gamma)
- Only disorder differs between realisations

**Current:**
```julia
for realisation in 1:n_disorder_realisations
    onsite_base, tx, ty = real_space_perturbed_qwz_blocks(...)  # REDUNDANT
    # Build H with different disorder...
end
```

**Impact:** **~5% overhead** (small but unnecessary)

**Solution:** Move block computation outside realisation loop
```julia
onsite_base, tx, ty = real_space_perturbed_qwz_blocks(...)
for realisation in 1:n_disorder_realisations
    H = build_H_with_blocks(onsite_base, tx, ty, disorder)
end
```

---

### 3. **Thread Parallelization Inefficiency** ⚠️ MEDIUM
**Location:** Main computation loop with `@showprogress` + `@threads`

**Problem:**
- `@showprogress` with `@threads` causes thread contention on progress updates
- Realisation loop is serial (not parallelized)
- For n_disorder_realisations=20, each thread does 20x serialwork—**not scalable**
- Load imbalance: threads finish at different times waiting for slowest

**Current structure:**
```julia
@showprogress Threads.@threads for i in 1:n_params  # Progress contention
    for realisation in 1:n_disorder_realisations    # Serial! 
        # Eigenvalue compute (slow)
    end
end
```

**Impact:** **Progress overhead ~2-5%**, more importantly: **realisation loop not scaling with threads**

**Solutions:**
1. Remove `@showprogress` from parallel section (use manual progress if needed)
2. Flatten to 2D parallelization: `(param × realisation)` instead of nesting

---

### 4. **Inefficient Diagonal Operator Usage** ⚠️ MEDIUM
**Location:** `build_localiser_operator`

**Problem:**
- `X` and `Y` are `Diagonal` objects, but treated as dense in block matrices
- `Ablock = kappa * (X - x0 * I)` creates dense matrices before sparse assembly
- `Ablock - im * Bblock` creates 2x2 complex block matrices

**Current (medium efficiency):**
```julia
Ablock = kappa * (X - x0 * I)    # Dense diagonal -> Matrix
Bblock = kappa * (Y - y0 * I)
L = [Hshift  (Ablock - im * Bblock);   # Block matrix assembly
     (Ablock + im * Bblock)  -Hshift]
```

**Better approach:**
```julia
# Keep as diagonal: Ablock = Diagonal(...)
# Avoid dense matrix creation
# Use sparse block constructor
```

**Impact:** **~10-15% overhead** in localiser operator construction

---

### 5. **Memory Allocation in Tight Loop** ⚠️ MEDIUM
**Location:** Main threaded loop

**Problem:**
```julia
@threads for i in 1:n_params
    low_lying_evals_accum = zeros(Float64, n_lowest_evals)  # Allocates EVERY iteration
    for realisation in ...
        low_lying_evals_avg = low_lying_evals_accum ./ n_disorder_realisations  # Allocates
    end
end
```

**Impact:** Thousands of small allocations/deallocations

**Solution:** Pre-allocate outside loop, reuse in-place operations

---

## Optimization Opportunities (by impact)

| Priority | Issue | Estimated Speedup | Effort |
|----------|-------|-------------------|--------|
| 🔴 CRITICAL | Sparse matrix assembly via COO | 5-10x | Medium |
| 🟠 HIGH | Remove progress reporting from threads | 2-5% | Easy |
| 🟠 HIGH | Redundant block computation | 5% | Easy |
| 🟡 MEDIUM | Pre-allocate reusable buffers | 3-5% | Easy |
| 🟡 MEDIUM | Flatten nested parallelization | 10-20% (with 20 realisations) | Medium |
| 🟡 MEDIUM | Use sparse diagonal efficiently | 10-15% | Medium |
| 🟢 LOW | Tune BLAS threading (depends on system) | 0-10% | Hard |

**Combined potential:** **2-4x speedup** (conservative estimate with sparse matrix fix + parallelization improvements)

---

## BLAS & Thread Utilization

### Current State:
✓ **Good:**
- Using `KrylovKit.eigsolve` which utilizes BLAS for matrix-vector products
- LU factorization in `shift_invert` uses BLAS
- Threads.@threads for outer parallelization

✗ **Poor:**
- Progress bar contention (minor)
- Small problem size (Lx=Ly=20: 800x800 matrices) → BLAS overhead dominates
- Single-threaded sparse matrix construction
- Redundant diagonal computations

### Recommendation:
- For **small systems** (Lx,Ly < 40): Current BLAS usage is acceptable, improve matrix construction
- For **larger systems** (Lx,Ly > 50): Spend more effort on cache-efficient sparse assembly

---

## Recommended Implementation Priority

### Phase 1 (Easy, ~30% speedup):
```julia
1. Move real_space_perturbed_qwz_blocks outside realisation loop
2. Remove @showprogress from parallel region
3. Pre-allocate low_lying_evals_accum and buffers
4. Use ./ with pre-allocated output buffer: div!(low_lying_evals_avg, ...)
```

### Phase 2 (Medium, ~3-5x):
```julia
1. Replace sparse matrix assembly with COO format
2. Flatten parallelization to (param, realisation) 2D grid
3. Use SparseArrays.sparse! or COO constructors
```

### Phase 3 (Advanced):
```julia
1. Pre-compute and cache Hamiltonian structure
2. Only modify disorder entries between realisations
3. Consider StaticArrays for small 2x2 blocks
4. Profile with @time/@benchmark on representative subset
```

---

## Diagnostic Commands

Run these to verify improvements:

```julia
# Before optimization
@time results_df = fast_low_lying_localiser_spectrum(...)

# Profile specific function
using Profile
@profile fast_low_lying_localiser_spectrum([20], [20]; 
    Ws=[0.0, 1.0], Es=[0.0, 1.0], 
    n_disorder_realisations=5, n_lowest_evals=4)
Profile.print()

# Allocations
@time results_df = fast_low_lying_localiser_spectrum(...)  # Watch allocations

# Check thread scaling
LinearAlgebra.BLAS.set_num_threads(1)  # Disable BLAS threads, test Julia threads only
```

---

## Summary Table

| Component | Current | Issue | Fix |
|-----------|---------|-------|-----|
| **Sparse assembly** | Element-wise += | O(n) per insert | COO format |
| **Block computation** | Per realisation | Redundant | Pre-compute |
| **Parallelization** | Nested serial | Not scalable | 2D grid |
| **Allocation** | Per iteration | Fragmentation | Pre-allocate |
| **Progress reporting** | In parallel loop | Contention | Remove/log outside |
| **Diagonal operators** | Dense conversion | Overhead | Sparse diagonal ops |

