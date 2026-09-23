#!/usr/bin/env julia
"""
Diagnose Hermiticity of embedding implementations
"""

using LinearAlgebra
using LaTeXStrings
using Printf
include("notebook_calcs/qwz_pbc_module.jl")
using .QWZ_Model

println("=" ^ 70)
println("HERMITICITY DIAGNOSTIC")
println("=" ^ 70)

A, B, m, gamma = 1.0, 1.0, -1.0, 0.5
perturbation_type = :symmetric
d_test = 1.0
phi_test = π/4

# Test point
kx, ky = 0.5, 0.3

println("\nTest Point: (kx=$(kx), ky=$(ky))")
println("-" ^ 70)

# 1. No embedding - should be Hermitian
H_no_embed = k_space_perturbed_qwz_hamiltonian(
    kx, ky; A=A, B=B, m=m, gamma=gamma, 
    perturbation_type=perturbation_type,
    orbital_displacement=0.0
)

is_herm_no_embed = norm(H_no_embed - H_no_embed') < 1e-14
println("1. No embedding: Hermitian? $(is_herm_no_embed)")
println("   Eigenvalues: $(round.(eigvals(H_no_embed), digits=4))")
println("   Eigenvalues type: $(typeof.(eigvals(H_no_embed)))")

# 2. With d=0 - should still be Hermitian
H_d0 = k_space_perturbed_qwz_hamiltonian(
    kx, ky; A=A, B=B, m=m, gamma=gamma, 
    perturbation_type=perturbation_type,
    orbital_displacement=0.0, phi=phi_test
)

is_herm_d0 = norm(H_d0 - H_d0') < 1e-14
println("\n2. With d=0, φ=$(round(phi_test,digits=3)): Hermitian? $(is_herm_d0)")
println("   Eigenvalues: $(round.(eigvals(H_d0), digits=4))")

# 3. With d≠0 - CRITICAL TEST
H_d_nonzero = k_space_perturbed_qwz_hamiltonian(
    kx, ky; A=A, B=B, m=m, gamma=gamma, 
    perturbation_type=perturbation_type,
    orbital_displacement=d_test, phi=phi_test
)

is_herm_d_nonzero = norm(H_d_nonzero - H_d_nonzero') < 1e-14
println("\n3. With d=$(d_test), φ=$(round(phi_test,digits=3)): Hermitian? $(is_herm_d_nonzero)")
println("   ||H - H†|| = $(norm(H_d_nonzero - H_d_nonzero'))")
println("   Eigenvalues: $(round.(eigvals(H_d_nonzero), digits=4))")

if !is_herm_d_nonzero
    println("\n   ⚠ WARNING: Hamiltonian is NOT Hermitian!")
    println("   This indicates the embedding implementation is incorrect.")
    diff_mat = H_d_nonzero - H_d_nonzero'
    println("   H - H† = \n$(diff_mat)")
end

println("\n" ^ 1)
println("Current EMBEDDING_TYPE: $(QWZ_Model.EMBEDDING_TYPE)")
println("=" ^ 70)

if QWZ_Model.EMBEDDING_TYPE == :orbital_displacement_phase
    println("\n⚠ DIAGNOSIS: The :orbital_displacement_phase implementation breaks Hermiticity")
    println("   The phase factor e^(i*k·d) is gauge-dependent and non-Hermitian.")
    println("   Switch to :sigma_z_rotation to test the unitary transformation approach.")
end
