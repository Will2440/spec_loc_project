# Example: Using orbital sublattice embedding (phi, orbital_displacement)
# to tune orbital position/embedding in the QWZ Hamiltonian

using LinearAlgebra

# Include the main script to get k_space_perturbed_qwz_hamiltonian
include("notebook_calcs/full_pbc_calc.jl")

# ============================================================================
# Test 1: Default (no orbital displacement)
# ============================================================================
println("Test 1: Default QWZ Hamiltonian (orbital_displacement = 0)")
H_default = k_space_perturbed_qwz_hamiltonian(
    π/4, π/4;  # kx, ky
    A=1.0, B=1.0, m=-1.0, gamma=0.0,
    orbital_displacement=0.0  # No embedding
)
println("Energy eigenvalues:", eigen(Hermitian(H_default)).values)
println()

# ============================================================================
# Test 2: Orbital displacement along x-direction (phi=0)
# ============================================================================
println("Test 2: Orbital displacement along x-direction (phi=0, magnitude=0.2)")
H_disp_x = k_space_perturbed_qwz_hamiltonian(
    π/4, π/4;
    A=1.0, B=1.0, m=-1.0, gamma=0.0,
    phi=0.0,  # Along x-axis
    orbital_displacement=0.2
)
println("Energy eigenvalues:", eigen(Hermitian(H_disp_x)).values)
println()

# ============================================================================
# Test 3: Orbital displacement along y-direction (phi=π/2)
# ============================================================================
println("Test 3: Orbital displacement along y-direction (phi=π/2, magnitude=0.2)")
H_disp_y = k_space_perturbed_qwz_hamiltonian(
    π/4, π/4;
    A=1.0, B=1.0, m=-1.0, gamma=0.0,
    phi=π/2,  # Along y-axis
    orbital_displacement=0.2
)
println("Energy eigenvalues:", eigen(Hermitian(H_disp_y)).values)
println()

# ============================================================================
# Test 4: Diagonal displacement (phi=π/4)
# ============================================================================
println("Test 4: Diagonal orbital displacement (phi=π/4, magnitude=0.2)")
H_disp_diag = k_space_perturbed_qwz_hamiltonian(
    π/4, π/4;
    A=1.0, B=1.0, m=-1.0, gamma=0.0,
    phi=π/4,  # Diagonal (45°)
    orbital_displacement=0.2
)
println("Energy eigenvalues:", eigen(Hermitian(H_disp_diag)).values)
println()

# ============================================================================
# Test 5: Using compute_bulk_band_berry_data with orbital embedding
# ============================================================================
println("Test 5: Computing band structure with orbital embedding...")
energies_embedded = compute_bulk_band_berry_data(
    Nkx=51, Nky=51,
    A=1.0, B=1.0, m=-1.0, gamma=0.0,
    phi=π/6,  # 30° rotation
    orbital_displacement=0.15  # Magnitude
)
println("Computed band structure with orbital embedding (51×51 k-points)")
println("Energy range: [$(minimum(energies_embedded)), $(maximum(energies_embedded))]")
println()

println("✓ All tests completed successfully!")
println()
println("Parameter meanings:")
println("  phi (radians): Angular position of orbital displacement vector d(phi) = |d|*(cos(phi), sin(phi))")
println("  orbital_displacement (real): Magnitude |d| of the displacement")
println()
println("Applies transformation: H' = exp(-i*theta*σ_z) * H * exp(i*theta*σ_z)")
println("                       where theta = k_x*d_x(phi) + k_y*d_y(phi)")
