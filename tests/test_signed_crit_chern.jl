using Printf
include("notebook_calcs/qwz_pbc_module.jl")
using .QWZ_Model

function test_signed_crit_chern()
    println("Testing Signed (±) Critical Chern Energy Display...")
    println("=" ^ 70)

    # Create test data
    data = compute_bulk_band_berry_data(
        Nkx=51,
        Nky=51,
        A=1.0,
        B=1.0,
        m=-1.0,
        gamma=0.0,
        phi=0.0,
        orbital_displacement=0.0,
        perturbation_type=:none
    )

    # Test parameters
    c_val_1 = 0.5
    c_val_2 = 0.0

    println("\nTest Configuration:")
    println("  Critical Chern 1 (nominal) = $c_val_1")
    println("  Critical Chern 2 (nominal) = $c_val_2")
    println("\n" ^ 1)
    println("=" ^ 70)

    nbands = size(data.berry_curvature, 1)
    
    println("\nBand-Resolved Signed Critical Chern Energies:")
    println("-" ^ 70)
    
    for b in 1:nbands
        println("\nBand $b:")
        
        # Get per-band data
        E_band = vec(data.plaquette_energies[b, :, :])
        C_band = vec(data.cum_chern_per_band[b, :, :])
        
        # Sort by energy
        sort_idx = sortperm(E_band)
        E_band_sorted = E_band[sort_idx]
        C_band_sorted = C_band[sort_idx]
        
        # Print range info
        @printf("  Chern range: [%.4f, %.4f]\n", minimum(C_band_sorted), maximum(C_band_sorted))
        
        # Find energies for both ± critical Chern 1
        E_crit_1_pos = get_energy_at_chern(E_band_sorted, C_band_sorted, c_val_1)
        E_crit_1_neg = get_energy_at_chern(E_band_sorted, C_band_sorted, -c_val_1)
        
        # Find energies for both ± critical Chern 2
        E_crit_2_pos = get_energy_at_chern(E_band_sorted, C_band_sorted, c_val_2)
        E_crit_2_neg = get_energy_at_chern(E_band_sorted, C_band_sorted, -c_val_2)
        
        # Display Critical Chern 1
        println("\n  Critical Chern 1:")
        @printf("    C = +%.4f: ", c_val_1)
        if isnan(E_crit_1_pos)
            println("N/A (outside range)")
        else
            @printf("E = %.4f\n", E_crit_1_pos)
        end
        
        @printf("    C = -%.4f: ", c_val_1)
        if isnan(E_crit_1_neg)
            println("N/A (outside range)")
        else
            @printf("E = %.4f\n", E_crit_1_neg)
        end
        
        # Display Critical Chern 2
        println("\n  Critical Chern 2:")
        @printf("    C = +%.4f: ", c_val_2)
        if isnan(E_crit_2_pos)
            println("N/A (outside range)")
        else
            @printf("E = %.4f\n", E_crit_2_pos)
        end
        
        @printf("    C = -%.4f: ", c_val_2)
        if isnan(E_crit_2_neg)
            println("N/A (outside range)")
        else
            @printf("E = %.4f\n", E_crit_2_neg)
        end
    end

    println("\n" ^ 1)
    println("=" ^ 70)
    println("Signed critical Chern energy computation test completed! ✓")
end

test_signed_crit_chern()
