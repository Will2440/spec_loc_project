"""
    plot_kappa_scan(result_file::String; output_dir::String="")

Load a spectral localiser JLD2 result and create plots of:
  1. Gap vs κ
  2. Spectral localiser invariant vs κ
  3. Signature vs κ

If output_dir is empty (default), saves to:
    ./plots/<param_summary>/

If output_dir is provided, saves to that directory directly.
"""
function plot_kappa_scan(result_file::String; output_dir::String="")
    # Load data
    result = load_specloc_result(result_file)
    
    # Extract metadata for title and directory naming
    metadata = result.metadata
    A, B, m = metadata["A"], metadata["B"], metadata["m"]
    gamma = result.axes["gammas"][1]
    W = result.axes["Ws"][1]
    pert_type = metadata["perturbation_type"]
    disorder_type = metadata["disorder_type"]
    orbital_d = metadata["orbital_displacement"]
    phi = metadata["phi"]
    
    # Determine output directory if not provided
    if isempty(output_dir)
        # Create parameter-based directory name
        param_str = @sprintf("A%.2f_B%.2f_m%.2f_g%.2f_W%.2f_pt%s_dt%s_d%.2f_phi%.2f",
                             A, B, m, gamma, W, pert_type, disorder_type, orbital_d, phi)
        output_dir = joinpath(@__DIR__, "plots", param_str)
    end
    
    # Extract data vs κ
    κs, gap_vs_κ = extract_gap_vs_kappa(result)
    _, idx_vs_κ = extract_index_vs_kappa(result)
    _, sig_vs_κ = extract_signature_vs_kappa(result)
    
    # Create output directory if needed
    mkpath(output_dir)
    
    # Generate parameter-based filenames
    param_base = @sprintf("pt%s_A%.2f_B%.2f_m%.2f_g%.2f_W%.2f_d%.2f",
                          pert_type, A, B, m, gamma, W, orbital_d)
    
    # Plot 1: Gap vs κ
    p1 = plot(
        κs, gap_vs_κ;
        xlabel="κ (log scale)", ylabel="Gap (E-averaged)",
        title="Spectral Localiser Gap vs κ\nA=$A, B=$B, m=$m, γ=$gamma, W=$W",
        marker=:circle, markersize=4, linewidth=2,
        legend=false, size=(800, 600),
        dpi=150, xscale=:log10,
        xticks=([1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1], ["10⁻⁶", "10⁻⁵", "10⁻⁴", "10⁻³", "10⁻²", "10⁻¹"])
    )
    gap_file = joinpath(output_dir, "01_gap_vs_kappa_$(param_base).png")
    savefig(p1, gap_file)
    println("✓ Saved: $(gap_file)")
    
    # Plot 2: Index vs κ
    p2 = plot(
        κs, idx_vs_κ;
        xlabel="κ (log scale)", ylabel="Spectral Localiser Index",
        title="Spectral Localiser Invariant vs κ\nA=$A, B=$B, m=$m, γ=$gamma, W=$W",
        marker=:circle, markersize=4, linewidth=2,
        legend=false, size=(800, 600),
        dpi=150, xscale=:log10,
        xticks=([1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1], ["10⁻⁶", "10⁻⁵", "10⁻⁴", "10⁻³", "10⁻²", "10⁻¹"])
    )
    idx_file = joinpath(output_dir, "02_index_vs_kappa_$(param_base).png")
    savefig(p2, idx_file)
    println("✓ Saved: $(idx_file)")
    
    # Plot 3: Signature vs κ
    p3 = plot(
        κs, sig_vs_κ;
        xlabel="κ (log scale)", ylabel="Spectral Localiser Signature",
        title="Spectral Localiser Signature vs κ\nA=$A, B=$B, m=$m, γ=$gamma, W=$W",
        marker=:circle, markersize=4, linewidth=2,
        legend=false, size=(800, 600),
        dpi=150, xscale=:log10,
        xticks=([1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1], ["10⁻⁶", "10⁻⁵", "10⁻⁴", "10⁻³", "10⁻²", "10⁻¹"])
    )
    sig_file = joinpath(output_dir, "03_signature_vs_kappa_$(param_base).png")
    savefig(p3, sig_file)
    println("✓ Saved: $(sig_file)")
    
    # Summary info
    println("\n" * "="^60)
    println("Spectral Localiser κ-scan Summary")
    println("="^60)
    println("System parameters:")
    @printf("  A=%.4f  B=%.4f  m=%.4f\n", A, B, m)
    @printf("  γ=%.4f  W=%.4f\n", gamma, W)
    println("  Perturbation type: $pert_type")
    println("  Disorder type: $disorder_type")
    println("  Orbital displacement: $(metadata["orbital_displacement"])")
    println("  Orbital angle φ: $(metadata["phi"])")
    
    println("\nData shapes:")
    @printf("  Number of κ values: %d\n", length(κs))
    @printf("  κ range: [%.6f, %.6f]\n", minimum(κs), maximum(κs))
    
    println("\nInvariant statistics:")
    @printf("  Min index:  %.6f\n", minimum(idx_vs_κ))
    @printf("  Max index:  %.6f\n", maximum(idx_vs_κ))
    @printf("  Mean index: %.6f\n", mean(idx_vs_κ))
    
    println("\nGap statistics (E-averaged):")
    @printf("  Min gap:  %.6e\n", minimum(gap_vs_κ))
    @printf("  Max gap:  %.6e\n", maximum(gap_vs_κ))
    @printf("  Mean gap: %.6e\n", mean(gap_vs_κ))
    
    println("\nPlots saved to: $output_dir")
    println("="^60)
end

# Entry point if run as script
if abspath(PROGRAM_FILE) == @__FILE__
    if length(ARGS) < 1
        println("Usage: julia plot_kappa_scan.jl <result.jld2> [output_dir]")
        exit(1)
    end
    
    result_file = ARGS[1]
    output_dir = length(ARGS) >= 2 ? ARGS[2] : "./plots"
    
    if !isfile(result_file)
        println("Error: File not found: $result_file")
        exit(1)
    end
    
    plot_kappa_scan(result_file; output_dir=output_dir)
end
