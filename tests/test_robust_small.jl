#!/usr/bin/env julia
"""Quick test of robust error handling with minimal parameter set"""

using LinearAlgebra
using Printf
using DataFrames
using SparseArrays
using ProgressMeter
using Statistics
using KrylovKit
using Base.Threads

println("Julia Threads: ", Threads.nthreads())
println("BLAS Threads:  ", LinearAlgebra.BLAS.get_num_threads())

# Load only the helper functions from the main script
include("notebook_calcs/fast_specloc_optimised.jl")

println("\n✓ Script loaded successfully")
println("✓ Error handling functions available")

# Quick test with TINY parameter set
println("\n" * "="^80)
println("Running MINIMAL parameter test to verify robust error handling...")
println("="^80)

Avals = [1.0]
Bvals = [1.0]
mvals = [0.0]  # Just ONE value
gammavals = [0.0, 1.0]  # Only 2 values
Ws = [0.0]
Es = [0.0]  # Just ONE energy
kappas = [0.2]
embedding_ds = [0.0]
embedding_phis = [0.0]
Lxs = [10]  # Small system
Lys = [10]
n_lowest_evals = 4
n_disorder_realisations = 1

println("\nTest parameters:")
println("  Lx, Ly: 10×10")
println("  m values: 1")
println("  gamma values: 2")
println("  E values: 1")
println("  Total grid points: $(1*1*1*2*1*1*1) = 2")

try
    results_df = fast_low_lying_localiser_spectrum_with_chern(
        Lxs, Lys;
        Avals=Avals, Bvals=Bvals, mvals=mvals, gammas=gammavals,
        perturbation_type=:asym_sin_sum, disorder_type=:none, Ws=Ws,
        x0_sym=:centre, y0_sym=:centre, Es=Es, kappas=kappas,
        phis=embedding_phis, orbital_displacements=embedding_ds,
        periodic_x=false, periodic_y=false,
        n_disorder_realisations=n_disorder_realisations,
        n_lowest_evals=n_lowest_evals,
        keep_square=true,
        scale_kappa_to_L=false,
        kappa_scales=[0.0004],
        compute_chern=true,
        chern_method=:complex
    )
    
    println("\n" * "="^80)
    println("✓ TEST PASSED")
    println("="^80)
    println("Results DataFrame:")
    println("  Rows: $(nrow(results_df))")
    println("  Columns: $(ncol(results_df))")
    println("\nColumn names:")
    for name in names(results_df)
        println("  - $name")
    end
    
    # Check for any fallback cases
    if :eval_status in names(results_df)
        status_counts = combine(groupby(results_df, :eval_status), nrow => :count)
        println("\nEigenvalue solver status breakdown:")
        println(status_counts)
    end
    
catch err
    println("\n" * "="^80)
    println("✗ TEST FAILED")
    println("="^80)
    println("Error: ", err)
    showerror(stdout, err)
    println("\nFull backtrace:")
    stacktrace(catch_backtrace())
end
