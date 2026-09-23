using Printf
include("notebook_calcs/qwz_pbc_module.jl")
using .QWZ_Model

# Test the accumulated Chern computation
println("Testing Accumulated Chern Computation...")
println("=" ^ 60)

# Create test data with lower resolution for speed
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

# Test at various Fermi energies
test_fermi_energies = [-2.0, -1.0, -0.5, 0.0, 0.5, 1.0, 2.0]

println("\nAccumulated Chern C(E_F) at different Fermi energies:")
println("-" ^ 60)

for E_F in test_fermi_energies
    # Compute per-band accumulated Chern
    nbands = size(data.berry_curvature, 1)
    chern_accum_bands = zeros(Float64, nbands)
    
    for b in 1:nbands
        # Sum berry curvature for all k-points where E(k) <= E_F
        for ix in 1:size(data.plaquette_energies, 2), iy in 1:size(data.plaquette_energies, 3)
            if data.plaquette_energies[b, ix, iy] <= E_F
                chern_accum_bands[b] += data.berry_curvature[b, ix, iy] / (2π)
            end
        end
    end
    
    # Total accumulated Chern across all bands
    chern_accum_total = sum(chern_accum_bands)
    
    println("\nE_F = $E_F")
    @printf("  C(E_F) = %.4f\n", chern_accum_total)
    for b in 1:nbands
        @printf("  Band %d: %.4f\n", b, chern_accum_bands[b])
    end
end

println("\n" ^ 60)
println("Accumulated Chern computation test completed! ✓")
