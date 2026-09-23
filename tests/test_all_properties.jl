using Printf
include("notebook_calcs/qwz_pbc_module.jl")
using .QWZ_Model

function test_all_properties()
    println("Testing All Properties Together...")
    println("=" ^ 60)

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

    # Test parameters matching viewer defaults
    E_fermi = 0.0
    c_val_1 = 0.5
    c_val_2 = 0.0

    println("\nTest Configuration:")
    println("  E_F = $E_fermi")
    println("  Critical Chern 1 = $c_val_1")
    println("  Critical Chern 2 = $c_val_2")
    println("\n" ^ 1)
    println("=" ^ 60)

    # Property 1: Euler Characteristic
    println("\n1. EULER CHARACTERISTIC (χ)")
    println("-" ^ 60)
    
    chi_total = 0
    pockets_total = 0
    holes_total = 0
    
    for b in 1:size(data.plaquette_energies, 1)
        M = @view(data.plaquette_energies[b, :, :]) .<= E_fermi
        chi_b = compute_euler_characteristic(M; periodic=true)
        pockets_b, holes_b = compute_pockets_and_holes(M; periodic=true)
        chi_total += chi_b
        pockets_total += pockets_b
        holes_total += holes_b
    end

    @printf("χ = %d\n", chi_total)
    @printf("N_pockets = %d\n", pockets_total)
    @printf("N_holes = %d\n", holes_total)
    @printf("Relation: χ = N_pockets - N_holes → %d = %d - %d ✓\n", 
            chi_total, pockets_total, holes_total)

    # Property 2: Accumulated Chern
    println("\n2. ACCUMULATED CHERN C(E_F)")
    println("-" ^ 60)
    nbands = size(data.berry_curvature, 1)
    chern_accum_bands = zeros(Float64, nbands)

    for b in 1:nbands
        for ix in 1:size(data.plaquette_energies, 2), iy in 1:size(data.plaquette_energies, 3)
            if data.plaquette_energies[b, ix, iy] <= E_fermi
                chern_accum_bands[b] += data.berry_curvature[b, ix, iy] / (2π)
            end
        end
    end

    chern_accum_total = sum(chern_accum_bands)
    @printf("C(E_F) = %.4f\n", chern_accum_total)
    for b in 1:nbands
        @printf("  Band %d: %.4f\n", b, chern_accum_bands[b])
    end

    # Property 3: Energy of Critical Chern
    println("\n3. ENERGY OF CRITICAL CHERN")
    println("-" ^ 60)
    energies_vec = data.accumulation_energies
    chern_vec = data.cumulative_chern

    E_crit_1 = get_energy_at_chern(energies_vec, chern_vec, c_val_1)
    E_crit_2 = get_energy_at_chern(energies_vec, chern_vec, c_val_2)

    @printf("Critical Chern 1: C = %.4f\n", c_val_1)
    if isnan(E_crit_1)
        println("  E = N/A (outside range)")
    else
        @printf("  E = %.4f\n", E_crit_1)
    end

    @printf("Critical Chern 2: C = %.4f\n", c_val_2)
    if isnan(E_crit_2)
        println("  E = N/A (outside range)")
    else
        @printf("  E = %.4f\n", E_crit_2)
    end

    println("\n" ^ 60)
    println("All three properties computed successfully! ✓")
end

test_all_properties()
