#!/usr/bin/env julia
"""Test compilation of updated fast_specloc_optimised.jl"""

using LinearAlgebra
using SparseArrays
using KrylovKit
using DataFrames

try
    # Load only the functions, not the execution section
    include("notebook_calcs/fast_specloc_optimised.jl")
    println("✓ Script compiled successfully")
    
    # Test the helper function
    test_params = (A=1.0, B=1.0, m=-1.0, E=0.5, Lx=20)
    formatted = format_params_for_log(test_params)
    println("✓ format_params_for_log works: $formatted")
    
    # Check that compute_low_lying_localiser_shift_invert has correct signature
    sig = methods(compute_low_lying_localiser_shift_invert)
    println("✓ compute_low_lying_localiser_shift_invert methods: $(length(sig)) overload(s)")
    
catch err
    println("✗ Compilation/test error:")
    showerror(stdout, err)
    stacktrace(catch_backtrace())
end
