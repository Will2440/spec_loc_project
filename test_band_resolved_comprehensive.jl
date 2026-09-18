using Printf
include("notebook_calcs/qwz_pbc_module.jl")
using .QWZ_Model

function test_band_resolved_with_gamma()
    println("Testing Band-Resolved Critical Chern with Gamma Perturbation...")
    println("=" ^ 70)

    # Test with and without perturbation
    test_configs = [
        (0.0, "No perturbation (base QWZ)"),
        (0.5, "γ = 0.5 (symmetric perturbation)")
    ]
    
    c_val_1 = 0.5
    c_val_2 = 0.0
    
    for (gamma_val, description) in test_configs
        println("\n" * description)
        println("-" ^ 70)
        
        data = compute_bulk_band_berry_data(
            Nkx=51,
            Nky=51,
            A=1.0,
            B=1.0,
            m=-1.0,
            gamma=gamma_val,
            phi=0.0,
            orbital_displacement=0.0,
            perturbation_type=:sym_cos_sum
        )
        
        nbands = size(data.berry_curvature, 1)
        
        for b in 1:nbands
            # Get per-band data
            E_band = vec(data.plaquette_energies[b, :, :])
            C_band = vec(data.cum_chern_per_band[b, :, :])
            
            # Sort by energy
            sort_idx = sortperm(E_band)
            E_band_sorted = E_band[sort_idx]
            C_band_sorted = C_band[sort_idx]
            
            # Find energies
            E_crit_1 = get_energy_at_chern(E_band_sorted, C_band_sorted, c_val_1)
            E_crit_2 = get_energy_at_chern(E_band_sorted, C_band_sorted, c_val_2)
            
            println("\nBand $b:")
            @printf("  Chern range: [%.4f, %.4f]\n", minimum(C_band_sorted), maximum(C_band_sorted))
            
            @printf("  Crit Chern 1 (C=%.4f): ", c_val_1)
            if isnan(E_crit_1)
                println("E = N/A (outside range)")
            else
                @printf("E = %8.4f\n", E_crit_1)
            end
            
            @printf("  Crit Chern 2 (C=%.4f): ", c_val_2)
            if isnan(E_crit_2)
                println("E = N/A (outside range)")
            else
                @printf("E = %8.4f\n", E_crit_2)
            end
        end
    end

    println("\n" * "=" ^ 70)
    println("Comprehensive band-resolved test completed! ✓")
end

test_band_resolved_with_gamma()
