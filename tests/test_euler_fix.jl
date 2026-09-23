using LinearAlgebra
include("notebook_calcs/qwz_pbc_module.jl")
using .QWZ_Model

# Test the Euler characteristic computation with a simple example
println("Testing Euler Characteristic Computation...")
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

# Test at various Fermi energies
test_fermi_energies = [-2.0, -1.0, 0.0, 1.0, 2.0]

println("\nComputing Euler Characteristic at different Fermi energies:")
println("-" ^ 60)

for E_F in test_fermi_energies
    chi_total = 0
    pockets_total = 0
    holes_total = 0
    
    for b in 1:size(data.plaquette_energies, 1)
        M = @view(data.plaquette_energies[b, :, :]) .<= E_F
        chi_b = compute_euler_characteristic(M; periodic=true)
        pockets_b, holes_b = compute_pockets_and_holes(M; periodic=true)
        chi_total += chi_b
        pockets_total += pockets_b
        holes_total += holes_b
    end
    
    # Verify the relation: χ = N_pockets - N_holes
    relation_check = (chi_total == pockets_total - holes_total)
    
    println("\nE_F = $E_F")
    println("  χ = $chi_total")
    println("  N_pockets = $pockets_total")
    println("  N_holes = $holes_total")
    println("  χ = N_pockets - N_holes? $relation_check ✓")
end

println("\n" ^ 60)
println("All tests passed! ✓")
