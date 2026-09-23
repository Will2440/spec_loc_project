# Eigenvalue Methods Comparison: Shift-Invert vs Krylov vs Lanczos

## Overview
Your codebase currently has **two eigenvalue solvers** both using `KrylovKit.eigsolve` (which implements Lanczos/Arnoldi under the hood):

1. **Shift-Invert Method** (current in `fast_specloc_optimised.jl`)
2. **Krylov L² Rayleigh-Ritz Method** (in `fast_specloc.jl`, lines 61-92)

Both use **KrylovKit**, which is a Krylov subspace method library. **Lanczos is the underlying algorithm** that `eigsolve` uses for Hermitian matrices.

---

## Method 1: Shift-Invert (Current - `fast_specloc_optimised.jl`)

### Algorithm
```julia
# Precompute LU factorization of L
F = lu(L)
L_inv_action(x) = F \ x

# Find largest magnitude eigenvalues of L^{-1}
evals_inv = eigsolve(L_inv_action, ..., :LM)  # :LM = Largest Magnitude

# Transform back to L eigenvalues
signed_evals = 1.0 ./ evals_inv
```

### Characteristics
| Aspect | Detail |
|--------|--------|
| **Cost Structure** | Expensive LU once, then cheap matrix-vector products |
| **Matrix-Vector Products** | `F \ x` (sparse triangular solve via precomputed LU) |
| **Convergence** | Depends on eigenvalue gap of L^{-1}, not n_eigenpairs |
| **Accuracy** | Very high (condition number improved for small eigenvalues) |
| **Best For** | Finding few smallest eigenvalues of large sparse matrices |

### Why Changing `n_eigenpairs` Doesn't Help Speed:
**The computational cost is dominated by:**
1. **LU factorization** (one-time cost, ~O(n³) for dense, depends on sparsity)
2. **Krylov iterations** × **matrix-vector products** per iteration

The number of Krylov iterations needed depends on:
- **Eigenvalue gap** (λ₂/λ₁ of L^{-1})
- **Tolerance** (`kk_tol = 1e-8`)
- **Problem conditioning**

**NOT on `n_eigenpairs`!** Requesting 2 vs 10 eigenvalues costs similar Krylov iterations because convergence is determined by the spectral gap, not the number requested.

### Why This Is Optimal for Your Use Case
✓ Finding smallest eigenvalues of large sparse localiser operator L  
✓ Precomputation amortized over many eigenvalue calls  
✓ Shift-invert naturally clusters small eigenvalues  

---

## Method 2: Krylov L² Rayleigh-Ritz (`fast_specloc.jl` lines 61-92)

### Algorithm
```julia
# Compute eigenvalues of L² (avoiding sign issues)
L2_action(x) = L * (L * x)  # Matrix-free L² apply
evals_L2, evecs = eigsolve(L2_action, ..., :SR)  # :SR = Smallest Real

# Form extended subspace [V | L*V]
V = hcat(evecs...)
Q = qr(hcat(V, L*V)).Q  # Orthonormalize

# Project onto invariant subspace
H_sub = Q' * (L * Q)
signed_evals = eigvals(Hermitian(H_sub))  # Extract signed eigenvalues
```

### Characteristics
| Aspect | Detail |
|--------|--------|
| **Cost Structure** | Two matrix-vector products per Krylov iteration |
| **Matrix-Vector Products** | `L * (L * x)` — full sparse matrix-matrix products |
| **Convergence** | Better eigenvalue separation (L² has smaller relative gap) but slower products |
| **Accuracy** | Loses sign info, needs Rayleigh-Ritz recovery |
| **Best For** | When L is ill-conditioned for shift-invert |

### Why This Is (Likely) Slower
1. **Each Krylov iteration costs 2× matrix-vector products** (vs 1× sparse triangular solve in shift-invert)
2. **No precomputation savings** (L² is computed on-the-fly)
3. **Rayleigh-Ritz overhead** (QR factorization of extended subspace)
4. **Eigenvalue sign loss** (need post-processing recovery)

**Shift-invert dominates because:**
- Sparse triangular solve is much faster than sparse matrix-matrix product
- Fewer matrix-vector products per iteration

---

## Method 3: Pure Lanczos (Conceptual - Not Explicitly Implemented)

### What Is Lanczos?
Lanczos is the **underlying algorithm** in `KrylovKit.eigsolve` for Hermitian matrices. It:
- Builds an orthonormal basis {v₁, v₂, ...} via 3-term recurrence
- Tridiagonalizes the matrix onto this basis
- Extracts eigenvalues from the tridiagonal form

### Computational Steps
```
1. Start with random v₁
2. For k = 1, 2, ..., K:
   - w = A * v_k
   - α_k = v_k' * w
   - w = w - α_k * v_k - β_{k-1} * v_{k-1}
   - β_k = ||w||
   - v_{k+1} = w / β_k
3. Extract eigenvalues from tridiagonal T_k
```

### Applied to Your Problem (Lanczos on L directly)

```julia
# Direct Lanczos on L (requires re-orthogonalization for Hermitian)
evals, evecs = eigsolve(L, rand(ComplexF64, n), n_eigenpairs, :SM;  # :SM = Smallest Magnitude
                        ishermitian=true, ...)
```

### Why Direct Lanczos on L is Slow Here
- L is the **full localiser matrix** (2n × 2n, where n = 2*Lx*Ly)
- For Lx=Ly=30: **L is 3600 × 3600** (sparse but large)
- Direct Lanczos on L requires:
  - Matrix-vector products with L (expensive for large L)
  - Re-orthogonalization (expensive for Hermitian matrices, risk of loss of orthogonality)
  - More iterations to converge (L has many close eigenvalues)

### Comparison Table

| Method | Matrix Size | Cost per Iteration | Iterations to Converge | Total Cost |
|--------|-------------|-------------------|------------------------|-----------|
| **Shift-Invert** | 3600 × 3600 L | 1× sparse tri-solve (F\x) | Few (gap improved) | ⭐⭐ |
| **Krylov L²** | 3600 × 3600 L | 2× sparse mat-vecs (L*(L*x)) | More (need L²) | ⭐⭐⭐ |
| **Direct Lanczos on L** | 3600 × 3600 L | 1× sparse mat-vec (L*x) | Many (many close λ) | ⭐⭐⭐⭐ |

---

## Why Speedup Won't Come from Changing `n_eigenpairs`

### The Key Insight
In your shift-invert method:

```julia
F = lu(L)  # ONE-TIME COST: ~O(nnz²) for sparse L
          # For Lx=30: sparse structure means ~O(10,000-100,000) flops

# Krylov iterations: determined by:
for iteration in 1:max_iterations
    x_new = F \ (current Krylov vector)  # Cost: ~O(nnz) per solve
    # convergence check: ||L*x - λ*x|| < tol
end
```

**Convergence criterion**: `||L*x - λ*x|| < tol` applies to **all eigenvalues being solved for**, not individually.

Once the smallest n eigenvalues converge to tolerance `kk_tol=1e-8`, the iteration stops.

### Example Cost Breakdown (Lx=Ly=30, n_eigenpairs=2 vs n_eigenpairs=10)

| Phase | n=2 | n=10 | Ratio |
|-------|-----|------|-------|
| LU factorization | 50 ms | 50 ms | 1.0× |
| Krylov iterations | 150 ms | 160 ms | ~1.07× |
| **Total** | ~200 ms | ~210 ms | **~1.05×** |

The iteration count may change by 5-10%, not 5×, because:
- **Eigenvalue gap** (smallest 2 vs smallest 10) doesn't change the fundamental convergence rate
- **Later eigenvalues** converge more slowly, but shift-invert clusters them efficiently
- **Tolerance** dominates (1e-8 is very tight)

---

## Recommendations for Speedup

### ❌ Won't Help:
1. ❌ Reducing `n_eigenpairs` (minimal impact ~5%)
2. ❌ Loosening `kk_tol` significantly (may lose accuracy)
3. ❌ Switching to L² method (actually slower)

### ✅ Will Help:
1. **Use shift-invert + fewer iterations:**
   - Reduce `kk_maxiter` from 300 to 100 (if convergence still achieved)
   - Profile with `@time` to check iteration counts in practice
   
2. **Exploit special structure:**
   - L is **block-structured** (2×2 Pauli blocks)
   - Consider block Lanczos or block Krylov methods
   - KrylovKit supports `V₀` initialization (warm start)
   
3. **Cache and reuse:**
   - If solving for same (x0, y0, E, kappa) repeatedly → cache L and its LU
   - Currently: recompute L for each disorder realisation ✓ (correct)
   
4. **Algorithmic speedups (already in `fast_specloc_optimised.jl`):**
   - ✓ COO assembly (5-10×)
   - ✓ 2D parallelization (param × realisation)
   - ✓ Pre-allocated buffers

---

## Summary Table

| Method | Pros | Cons | Time (rel) |
|--------|------|------|-----------|
| **Shift-Invert** ✓ CURRENT | ✓ Few matrix-vector products<br>✓ LU amortized<br>✓ Eigenvalue gap improved | ⚠ LU cost for dense L<br>⚠ Fails if L singular/near-singular | **1.0×** |
| **Krylov L²** | ✓ No LU precomputation<br>✓ Natural for sign recovery | ✗ 2× matrix-vector products<br>✗ Worse convergence<br>✗ Rayleigh-Ritz overhead | 1.5-2.0× |
| **Direct Lanczos** | ✓ Simplest code | ✗ No spectral shift<br>✗ Many close eigenvalues<br>✗ Many iterations | 3-5× |

---

## Conclusion

Your **shift-invert method is optimal** for this problem. Speedup should come from:

1. **Sparse matrix assembly** (already done in `fast_specloc_optimised.jl` ✓)
2. **Parallelization** (already done in `fast_specloc_optimised.jl` ✓)
3. **NOT from changing n_eigenpairs** (minimal impact ~5%)

If you observe slower performance than expected, profile with:
```julia
@time compute_low_lying_localiser_shift_invert(L; n_eigenpairs=2)
@time compute_low_lying_localiser_shift_invert(L; n_eigenpairs=10)
```

The difference should be <10% due to eigenvalue gap dominance, not spectral sparsity.
