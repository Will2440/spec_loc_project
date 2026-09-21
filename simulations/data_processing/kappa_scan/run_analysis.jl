using JLD2
using Statistics
using Printf
using Plots

"""
Quick runner script for kappa scan analysis.

Usage:
    julia run_analysis.jl <path_to_result.jld2>

This will:
  1. Load the JLD2 file
  2. Extract gap and invariant data vs κ
  3. Create plots and save to ./plots/
  4. Print summary statistics
"""

# Add project root to path so we can find load_data.jl
script_dir = dirname(abspath(@__FILE__))
push!(LOAD_PATH, script_dir)

include("load_data.jl")
include("plot_kappa_scan.jl")

function main()
    if length(ARGS) < 1
        println("""
        Spectral Localiser κ-scan Plotter
        =====================================
        
        Usage:
            julia run_analysis.jl <path_to_result.jld2>
        
        Example:
            julia run_analysis.jl ../../../data_collection/hpc/results/20260921_144508_r1/row1_chunk1_A1.0_B1.0_m-1.0_g0.0_ptsym_cos_diff_dtnone_d0.0_phi0.0_sx50_sy50.jld2
        """)
        exit(1)
    end
    
    result_file = ARGS[1]
    # If user provides a second argument, use it as custom output_dir; otherwise let plot_kappa_scan() auto-organize
    output_dir = length(ARGS) >= 2 ? ARGS[2] : ""
    
    if !isfile(result_file)
        println("❌ Error: File not found: $result_file")
        exit(1)
    end
    
    println("📂 Loading: $result_file")
    plot_kappa_scan(result_file; output_dir=output_dir)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
