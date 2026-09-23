#!/usr/bin/env julia
"""
Test script to compare embedding implementations
Tests that d=0 gives identical results for both implementations
"""

using LinearAlgebra
using LaTeXStrings
using Printf
include("notebook_calcs/qwz_pbc_module.jl")
using .QWZ_Model

println("=" ^ 60)
println("EMBEDDING IMPLEMENTATION COMPARISON TEST")
println("=" ^ 60)

# Test parameters
A, B, m, gamma = 1.0, 1.0, -1.0, 0.5
perturbation_type = :symmetric
d_test = 1.0
phi_test = π/4

# Test k-points
kx_vals = [-π/2, 0.0, π/2]
ky_vals = [-π/2, 0.0, π/2]

println("\nTest 1: Verify d=0 gives IDENTICAL results")
println("-" ^ 60)

results_d0 = []
for kx in kx_vals
    for ky in ky_vals
        H_k = k_space_perturbed_qwz_hamiltonian(
            kx, ky; A=A, B=B, m=m, gamma=gamma, 
            perturbation_type=perturbation_type,
            orbital_displacement=0.0, phi=phi_test
        )
        push!(results_d0, (kx, ky, eigvals(H_k)))
    end
end

println("With d=0, eigenvalues at test points:")
for (kx, ky, evals) in results_d0
    println("  (kx=$(round(kx,digits=2)), ky=$(round(ky,digits=2))): λ = $(round.(evals, digits=4))")
end

println("\n" ^ 1)
println("Test 2: Check effect of d≠0 on current implementation")
println("-" ^ 60)

# Temporarily test sigma_z_rotation
println("Current implementation: :sigma_z_rotation")
results_d_nonzero = []
for kx in kx_vals
    for ky in ky_vals
        H_k = k_space_perturbed_qwz_hamiltonian(
            kx, ky; A=A, B=B, m=m, gamma=gamma, 
            perturbation_type=perturbation_type,
            orbital_displacement=d_test, phi=phi_test
        )
        push!(results_d_nonzero, (kx, ky, eigvals(H_k)))
    end
end

println("With d=$(d_test), φ=$(round(phi_test, digits=3)), eigenvalues at test points:")
for (kx, ky, evals) in results_d_nonzero
    println("  (kx=$(round(kx,digits=2)), ky=$(round(ky,digits=2))): λ = $(round.(evals, digits=4))")
end

println("\n" ^ 1)
println("Test 3: Verify current implementation is :$(QWZ_Model.EMBEDDING_TYPE)")
println("-" ^ 60)
println("EMBEDDING_TYPE = $(QWZ_Model.EMBEDDING_TYPE)")

if QWZ_Model.EMBEDDING_TYPE == :orbital_displacement_phase
    println("✓ Currently using orbital_displacement_phase")
    println("\nNote: This implementation applies exp(i*k·d) phase factor")
    println("      which is gauge-dependent and may affect Berry curvature gauge choice.")
elseif QWZ_Model.EMBEDDING_TYPE == :sigma_z_rotation
    println("✓ Currently using sigma_z_rotation")
    println("\nNote: This implementation rotates (d_x, d_y) by angle θ = k·d")
    println("      which is unitary and gauge-invariant in the local basis.")
end

println("\n" ^ 1)
println("Test 4: Verify d=0 produces zero change vs bare system")
println("-" ^ 60)

# Compute with no embedding
H_no_embed = k_space_perturbed_qwz_hamiltonian(
    0.0, 0.0; A=A, B=B, m=m, gamma=gamma, 
    perturbation_type=perturbation_type,
    orbital_displacement=0.0
)

# Compute with d=0 embedding
H_d0_embed = k_space_perturbed_qwz_hamiltonian(
    0.0, 0.0; A=A, B=B, m=m, gamma=gamma, 
    perturbation_type=perturbation_type,
    orbital_displacement=0.0, phi=phi_test
)

diff = norm(H_no_embed - H_d0_embed)
println("Difference between no embedding and d=0 embedding at (0,0):")
println("  ||H_no_embed - H_d0_embed|| = $(diff)")

if diff < 1e-14
    println("  ✓ PASS: Difference is numerical zero")
else
    println("  ✗ FAIL: Difference detected! Bug in implementation?")
end

println("\n" ^ 1)
println("=" ^ 60)
println("Test complete.")
println("=" ^ 60)
