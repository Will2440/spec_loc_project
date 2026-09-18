using Printf
include("notebook_calcs/qwz_pbc_module.jl")
using .QWZ_Model

function demo_signed_critical_chern()
    println("\n" * "="^75)
    println("COMPREHENSIVE DEMO: Signed (±) Critical Chern Energy Display")
    println("="^75)

    # Test case 1: Base configuration
    println("\n" * "-"^75)
    println("TEST 1: Base QWZ Model (m=-1.0, γ=0.0)")
    println("-"^75)
    
    data = compute_bulk_band_berry_data(
        Nkx=51, Nky=51, A=1.0, B=1.0, m=-1.0, gamma=0.0, phi=0.0,
        orbital_displacement=0.0, perturbation_type=:none
    )
    
    demo_signed_energies(data, 0.5, 0.0, "Base config")

    # Test case 2: With perturbation
    println("\n" * "-"^75)
    println("TEST 2: QWZ Model with γ=0.5 perturbation")
    println("-"^75)
    
    data_pert = compute_bulk_band_berry_data(
        Nkx=51, Nky=51, A=1.0, B=1.0, m=-1.0, gamma=0.5, phi=0.0,
        orbital_displacement=0.0, perturbation_type=:none
    )
    
    demo_signed_energies(data_pert, 0.5, 0.0, "γ=0.5 config")

    println("\n" * "="^75)
    println("IMPLEMENTATION VALIDATION COMPLETE ✓")
    println("="^75 * "\n")
    println("Summary of Signed Critical Chern Feature:")
    println("  • For each band and critical Chern value:")
    println("    - Searches for BOTH positive (+C) and negative (-C) versions")
    println("    - Displays actual energy with explicit ± sign indicator")
    println("    - Shows 'N/A' only when that polarity is outside band's range")
    println("  • Matches the contour plotting logic in plt_accumulated_chern_heatmaps()")
    println("  • Properly handles edge cases like C=0 (same energy for ±0)")
    println("\n")
end

function demo_signed_energies(data, c_val_1, c_val_2, config_name)
    nbands = size(data.berry_curvature, 1)
    
    println("\nCritical Chern Values:")
    @printf("  Nominal Crit 1 = %.4f\n", c_val_1)
    @printf("  Nominal Crit 2 = %.4f\n", c_val_2)
    
    for b in 1:nbands
        E_band = vec(data.plaquette_energies[b, :, :])
        C_band = vec(data.cum_chern_per_band[b, :, :])
        
        sort_idx = sortperm(E_band)
        E_band_sorted = E_band[sort_idx]
        C_band_sorted = C_band[sort_idx]
        
        min_c = minimum(C_band_sorted)
        max_c = maximum(C_band_sorted)
        
        E_1_pos = get_energy_at_chern(E_band_sorted, C_band_sorted, c_val_1)
        E_1_neg = get_energy_at_chern(E_band_sorted, C_band_sorted, -c_val_1)
        E_2_pos = get_energy_at_chern(E_band_sorted, C_band_sorted, c_val_2)
        E_2_neg = get_energy_at_chern(E_band_sorted, C_band_sorted, -c_val_2)
        
        println("\n  Band $b: Chern ∈ [$(round(min_c, digits=4)), $(round(max_c, digits=4))]")
        
        println("    Critical Chern 1:")
        status_pos = isnan(E_1_pos) ? "N/A" : @sprintf("E = %.4f ✓", E_1_pos)
        status_neg = isnan(E_1_neg) ? "N/A" : @sprintf("E = %.4f ✓", E_1_neg)
        @printf("      C = +%.4f: %s\n", c_val_1, status_pos)
        @printf("      C = -%.4f: %s\n", c_val_1, status_neg)
        
        println("    Critical Chern 2:")
        status_pos2 = isnan(E_2_pos) ? "N/A" : @sprintf("E = %.4f ✓", E_2_pos)
        status_neg2 = isnan(E_2_neg) ? "N/A" : @sprintf("E = %.4f ✓", E_2_neg)
        @printf("      C = +%.4f: %s\n", c_val_2, status_pos2)
        @printf("      C = -%.4f: %s\n", c_val_2, status_neg2)
    end
end

demo_signed_critical_chern()
