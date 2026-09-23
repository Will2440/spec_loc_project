using Printf
include("notebook_calcs/qwz_pbc_module.jl")
using .QWZ_Model

function test_band_resolved_crit_chern()
    println("Testing Band-Resolved Energy of Critical Chern...")
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
    println("  Critical Chern 1 = $c_val_1")
    println("  Critical Chern 2 = $c_val_2")
    println("\n" ^ 1)
    println("=" ^ 70)

    nbands = size(data.berry_curvature, 1)
    
    println("\nBand-Resolved Critical Chern Energies:")
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
        @printf("  Energy range: [%.4f, %.4f]\n", minimum(E_band_sorted), maximum(E_band_sorted))
        @printf("  Chern range: [%.4f, %.4f]\n", minimum(C_band_sorted), maximum(C_band_sorted))
        
        # Find energies for both critical Chern values
        E_crit_1_band = get_energy_at_chern(E_band_sorted, C_band_sorted, c_val_1)
        E_crit_2_band = get_energy_at_chern(E_band_sorted, C_band_sorted, c_val_2)
        
        @printf("  Crit Chern 1 (C=%.4f): ", c_val_1)
        if isnan(E_crit_1_band)
            println("E = N/A (outside range)")
        else
            @printf("E = %.4f\n", E_crit_1_band)
        end
        
        @printf("  Crit Chern 2 (C=%.4f): ", c_val_2)
        if isnan(E_crit_2_band)
            println("E = N/A (outside range)")
        else
            @printf("E = %.4f\n", E_crit_2_band)
        end
    end

    println("\n" ^ 1)
    println("=" ^ 70)
    println("Band-resolved energy computation test completed! ✓")
end

test_band_resolved_crit_chern()
