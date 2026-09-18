module QWZ_SpecLoc_Model

using LinearAlgebra
using Printf
using Plots
using DataFrames
using SparseArrays
using ProgressMeter
using Statistics
using LaTeXStrings
using JLD2: @save, @load
using Base.Threads

export compute_bulk_band_berry_data, plt_bandstructure_heatmap
export plt_k_resolved_F_xy_heatmaps, plt_accumulated_chern_heatmaps
export real_space_perturbed_qwz_blocks, real_space_perturbed_hamiltonian_qwz
export k_space_perturbed_qwz_hamiltonian, unpack_saved_data, debug_accumulation_extrema


## Pauli matrices
const sigma_x = [0 1; 1 0]
const sigma_y = [0 -im; im 0]
const sigma_z = [1 0; 0 -1]
const identity = [1 0; 0 1]

pi_ticks = (
                [-π, -3π/4, -π/2, -π/4, 0, π/4, π/2, 3π/4, π],
                [L"-\pi", L"-3\pi/4", L"-\pi/2", L"-\pi/4", L"0", L"\pi/4", L"\pi/2", L"3\pi/4", L"\pi"]
            )

# 1. Unified Block Generator
function real_space_perturbed_qwz_blocks(; 
    A::Real=1.0, B::Real=1.0, m::Real=0.0, gamma::Real=0.0, perturbation_type::Symbol=:none
)
    onsite = ComplexF64.((m + 2.0 * B) * sigma_z)
    tx = ComplexF64.(-0.5 * B * sigma_z - 0.5im * A * sigma_x)
    ty = ComplexF64.(-0.5 * B * sigma_z - 0.5im * A * sigma_y)

    if perturbation_type == :symmetric
        tx += ComplexF64.(0.5 * gamma * identity)
        ty += ComplexF64.(0.5 * gamma * identity)
    elseif perturbation_type == :tilt
        tx += ComplexF64.(0.5im * gamma * identity)
    elseif perturbation_type != :none
        error("Unknown perturbation type: $perturbation_type")
    end

    return onsite, tx, ty
end

site_index_qwz(x::Int, y::Int, orb::Int, Lx::Int, Ly::Int) = 2 * ((y - 1) * Lx + (x - 1)) + orb

# 2. Updated Real-Space Hamiltonian Builder
function real_space_perturbed_hamiltonian_qwz(
    Lx::Int, Ly::Int; A::Real=1.0, B::Real=1.0, m::Real=0.0, gamma::Real=0.0,
    perturbation_type::Symbol=:none, periodic_x::Bool=false, periodic_y::Bool=false, sparse_output::Bool=true
)
    Lx > 0 || error("Lx must be positive.")
    Ly > 0 || error("Ly must be positive.")

    onsite, tx, ty = real_space_perturbed_qwz_blocks(; A=A, B=B, m=m, gamma=gamma, perturbation_type=perturbation_type)
    nsites = Lx * Ly
    dim = 2 * nsites
    H = sparse_output ? spzeros(ComplexF64, dim, dim) : zeros(ComplexF64, dim, dim)

    function add_block!(mat, row_site::Tuple{Int, Int}, col_site::Tuple{Int, Int}, block::AbstractMatrix{<:Number})
        (xr, yr) = row_site
        (xc, yc) = col_site
        row_base = site_index_qwz(xr, yr, 1, Lx, Ly)
        col_base = site_index_qwz(xc, yc, 1, Lx, Ly)
        @inbounds for a in 0:1, b in 0:1
            mat[row_base + a, col_base + b] += ComplexF64(block[a + 1, b + 1])
        end
        return nothing
    end

    for y in 1:Ly, x in 1:Lx
        add_block!(H, (x, y), (x, y), onsite)

        if x < Lx
            add_block!(H, (x + 1, y), (x, y), tx)
            add_block!(H, (x, y), (x + 1, y), tx')
        elseif periodic_x
            add_block!(H, (1, y), (x, y), tx)
            add_block!(H, (x, y), (1, y), tx')
        end

        if y < Ly
            add_block!(H, (x, y + 1), (x, y), ty)
            add_block!(H, (x, y), (x, y + 1), ty')
        elseif periodic_y
            add_block!(H, (x, 1), (x, y), ty)
            add_block!(H, (x, y), (x, 1), ty')
        end
    end

    return H
end

function k_space_perturbed_qwz_hamiltonian(
    kx::Real, ky::Real; A::Real=1.0, B::Real=1.0, m::Real=0.0, gamma::Real=0.0,
    perturbation_type::Symbol=:none, winding_number::Int=0
)::Matrix{ComplexF64}

    winding_number >= 0 || error("winding_number must be non-negative.")

    d_x = A * sin(kx)
    d_y = A * sin(ky)
    d_z = m + 2B - B * cos(kx) - B * cos(ky)
    scalar_term = 0.0

    if perturbation_type == :symmetric
        scalar_term += gamma * (cos(kx) + cos(ky))
    elseif perturbation_type == :tilt
        scalar_term += gamma * sin(kx)
    elseif perturbation_type != :none
        error("Unknown perturbation type: $perturbation_type")
    end

    if winding_number > 0
        n = winding_number
        cos_nkx = cos(n * kx)
        sin_nkx = sin(n * kx)
        d_x, d_y = (cos_nkx * d_x - sin_nkx * d_y, cos_nkx * d_y + sin_nkx * d_x)
    end

    H_k = scalar_term * identity + d_x * sigma_x + d_y * sigma_y + d_z * sigma_z
    return ComplexF64.(H_k)
end

function compute_bulk_band_berry_data(; Nkx::Int=101, Nky::Int=101, kwargs...)
    kx_vals = collect(range(-π, π; length=Nkx + 1))[1:end-1]
    ky_vals = collect(range(-π, π; length=Nky + 1))[1:end-1]
    nbands = 2

    energies = zeros(Float64, nbands, Nkx, Nky)
    eigenvectors = Array{ComplexF64}(undef, 2, nbands, Nkx, Nky)

    k_space_kwargs = (; [k => v for (k, v) in pairs(kwargs) if k ∈ (:A, :B, :m, :gamma, :perturbation_type, :winding_number)]...)

    for (ix, kx) in enumerate(kx_vals), (iy, ky) in enumerate(ky_vals)
        spectrum = eigen(Hermitian(k_space_perturbed_qwz_hamiltonian(kx, ky; k_space_kwargs...)))
        energies[:, ix, iy] .= spectrum.values
        eigenvectors[:, :, ix, iy] .= spectrum.vectors
    end

    berry_curvature = zeros(Float64, nbands, Nkx, Nky)
    plaquette_energies = zeros(Float64, nbands, Nkx, Nky)

    for ix in 1:Nkx, iy in 1:Nky
        ix_next = ix == Nkx ? 1 : (ix + 1)
        iy_next = iy == Nky ? 1 : (iy + 1)

        for band in 1:nbands
            u00 = @view eigenvectors[:, band, ix, iy]
            u10 = @view eigenvectors[:, band, ix_next, iy]
            u11 = @view eigenvectors[:, band, ix_next, iy_next]
            u01 = @view eigenvectors[:, band, ix, iy_next]

            link_x_00 = dot(u00, u10)
            link_y_10 = dot(u10, u11)
            link_x_01 = dot(u01, u11)
            link_y_00 = dot(u00, u01)

            berry_curvature[band, ix, iy] = angle(link_x_00 * link_y_10 / (link_x_01 * link_y_00))
            plaquette_energies[band, ix, iy] = 0.25 * (
                energies[band, ix, iy] + energies[band, ix_next, iy] +
                energies[band, ix_next, iy_next] + energies[band, ix, iy_next]
            )
        end
    end

    band_chern_numbers = vec(sum(berry_curvature, dims=(2, 3)) ./ (2π))

    cum_chern_per_band = zeros(Float64, nbands, Nkx, Nky)
    for band in 1:nbands
        flat_E = vec(plaquette_energies[band, :, :])
        flat_F = vec(berry_curvature[band, :, :]) ./ (2π)
        
        order = sortperm(flat_E)
        cum_sorted = cumsum(flat_F[order])
        
        cum_band_2d = zeros(Float64, Nkx, Nky)
        cum_band_2d[order] .= cum_sorted
        cum_chern_per_band[band, :, :] .= cum_band_2d
    end

    flat_energies_all = vec(plaquette_energies)
    flat_weights_all = vec(berry_curvature) ./ (2π)
    global_order = sortperm(flat_energies_all)

    accumulation_energies = flat_energies_all[global_order]
    accumulation_weights = flat_weights_all[global_order]
    cumulative_chern = cumsum(flat_weights_all[global_order])

    return (
        kx_vals=kx_vals, ky_vals=ky_vals, energies=energies, berry_curvature=berry_curvature,
        plaquette_energies=plaquette_energies, band_chern_numbers=band_chern_numbers,
        cum_chern_per_band=cum_chern_per_band, accumulation_energies=accumulation_energies,
        accumulation_weights=accumulation_weights, cumulative_chern=cumulative_chern,
    )
end




end # end module QWZ_SpecLoc_Model