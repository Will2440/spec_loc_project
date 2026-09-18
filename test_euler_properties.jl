#!/usr/bin/env julia

include("notebook_calcs/qwz_pbc_module.jl")
using .QWZ_Model

println("Testing Euler Characteristic functions...")

# Create some test data
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

println("Computed bulk band and berry data")

# Test Euler characteristic computation
E_mesh = range(-3.0, 3.0, length=51)
euler_data = compute_euler_characteristic_vs_energy(data; E_mesh=E_mesh, periodic=true)

println("\nEuler Characteristic Test Results:")
println("Energy mesh length: $(length(euler_data.E_mesh))")
println("Euler char length: $(length(euler_data.euler_char))")
println("Pockets length: $(length(euler_data.pockets))")
println("Holes length: $(length(euler_data.holes))")

# Test a single Fermi energy
E_fermi = 0.0
euler_single = compute_euler_characteristic_vs_energy(data; E_mesh=[E_fermi], periodic=true)
chi_fermi = euler_single.euler_char[1]
pockets_fermi = euler_single.pockets[1]
holes_fermi = euler_single.holes[1]

println("\nAt Fermi Energy E_f = $E_fermi:")
println("  χ = $chi_fermi")
println("  N_pockets = $pockets_fermi")
println("  N_holes = $holes_fermi")
println("  Verification: χ = N_pockets - N_holes = $(pockets_fermi - holes_fermi)")

if chi_fermi == (pockets_fermi - holes_fermi)
    println("\n✓ Euler characteristic relation verified!")
else
    println("\n✗ Warning: χ ≠ N_pockets - N_holes")
end

println("\nTest complete!")
