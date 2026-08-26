module QWZ_Model

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

function plt_bandstructure_heatmap(kx_vals::Vector{Float64}, ky_vals::Vector{Float64}, energies::Array{Float64, 3}; title::LaTeXString=L"", xlabel::LaTeXString=L"k_x", ylabel::LaTeXString=L"k_y", colour=:curl, share_colour_scale::Bool=false)
    c_limits = share_colour_scale ? extrema(energies[1:2, :, :]) : :auto
    plt_band1 = heatmap(kx_vals, ky_vals, energies[1, :, :], xlabel=xlabel, ylabel=ylabel, title="Bandstrucutre (Band 1 Lower)", xlims=(-pi, pi), ylims=(-pi, pi), xticks=pi_ticks, yticks=pi_ticks, color=colour, clims=c_limits, colorbar_title=L"E(k_x, k_y)", aspect_ratio=:equal)
    plt_band2 = heatmap(kx_vals, ky_vals, energies[2, :, :], xlabel=xlabel, ylabel=ylabel, title="Bandstrucutre (Band 2 Upper)", xlims=(-pi, pi), ylims=(-pi, pi), xticks=pi_ticks, yticks=pi_ticks, color=colour, clims=c_limits, colorbar_title=L"E(k_x, k_y)", aspect_ratio=:equal)
    return plot(plt_band1, plt_band2, layout=(1, 2), size=(1600, 600))
end

function plt_k_resolved_F_xy_heatmaps(kx_vals::Vector{Float64}, ky_vals::Vector{Float64}, berry_curvature::Array{Float64, 3}; title::LaTeXString=L"", xlabel::LaTeXString=L"k_x", ylabel::LaTeXString=L"k_y", colour=:RdBu, shift_to_centers::Bool=false, share_colour_scale::Bool=false)
    chern_density = berry_curvature ./ (2π)
    dkx = kx_vals[2] - kx_vals[1]
    dky = ky_vals[2] - ky_vals[1]
    x_coords = shift_to_centers ? (kx_vals .+ dkx / 2) : kx_vals
    y_coords = shift_to_centers ? (ky_vals .+ dky / 2) : ky_vals
    global_max = share_colour_scale ? maximum(abs, chern_density[1:2, :, :]) : 0.0

    plots = map(1:2) do band_index
        density_b = chern_density[band_index, :, :]
        c_limits = share_colour_scale ? (-global_max, global_max) : let max_val = maximum(abs, density_b); (-max_val, max_val) end

        c_min, c_max = c_limits[1], c_limits[2]
        c_positions = range(c_min, c_max, length=5)

        # Format each position into a scientific notation string
        c_labels = [@sprintf("%.1e", val) for val in c_positions]
        c_ticks = (c_positions, c_labels)

        total_chern = round(sum(density_b), digits=4)
        sub_title = "Flux k-resolved: " * L"F_{xy} / 2\pi \ C = %$(total_chern)" #L"\\text{Flux k-resolved: } F_{xy} / 2\pi C = %$(total_chern))"
        heatmap(
            x_coords, 
            y_coords, 
            density_b', 
            xlabel=xlabel, 
            ylabel=ylabel, 
            title=sub_title, 
            color=colour,
            clims=c_limits, 
            colorbar_ticks=c_ticks,
            xlims=(-pi, pi), 
            ylims=(-pi, pi), 
            xticks=pi_ticks, 
            yticks=pi_ticks, 
            aspect_ratio=:equal, 
            colorbar_title=L"F_{xy}(k_x, k_y) / 2\pi"
            )
    end
    return plot(plots[1], plots[2], layout=(1, 2), size=(1600, 600))
end

function plt_accumulated_chern_heatmaps(
    kx_vals::Vector{Float64}, 
    ky_vals::Vector{Float64}, 
    cum_chern_per_band::Array{Float64, 3}; 
    title::LaTeXString=L"", 
    xlabel::LaTeXString=L"k_x",
    ylabel::LaTeXString=L"k_y", 
    colour=:viridis, 
    shift_to_centers::Bool=false,
    add_contour::Bool=true,
    contour_value::Float64=0.5
)

    dkx = kx_vals[2] - kx_vals[1]
    dky = ky_vals[2] - ky_vals[1]
    x_coords = shift_to_centers ? (kx_vals .+ dkx / 2) : kx_vals
    y_coords = shift_to_centers ? (ky_vals .+ dky / 2) : ky_vals

    plots = map(1:2) do band_index
        cum_transposed = cum_chern_per_band[band_index, :, :]'
        min_val = minimum(cum_transposed)
        max_val = maximum(cum_transposed)
        total_chern = abs(max_val) > abs(min_val) ? round(max_val, digits=4) : round(min_val, digits=4)
        sub_title = "Chern number accumulated: " * L"C(E, k_x, k_y)" #L"\\text{Chern number accumulated: } C(E, k_x, k_y)"
        p = heatmap(
            x_coords, 
            y_coords, 
            cum_transposed, 
            xlabel=xlabel, 
            ylabel=ylabel, 
            title=sub_title, 
            color=colour, 
            clims=(min_val, max_val), 
            xlims=(-pi, pi), 
            ylims=(-pi, pi), 
            xticks=pi_ticks, 
            yticks=pi_ticks, 
            aspect_ratio=:equal, 
            colorbar_title=L"C(E, k_x, k_y)"
            )

        # Overlay ± Critical Contour Lines
        if add_contour
            target_val = abs(contour_value)
            
            # Collect both +val and -val if they exist in the band's value range
            levels_to_draw = Float64[]
            if min_val <= target_val <= max_val
                push!(levels_to_draw, target_val)
            end
            if min_val <= -target_val <= max_val
                push!(levels_to_draw, -target_val)
            end

            # Plot contours if any valid targets were found
            if !isempty(levels_to_draw)
                contour!(
                    p,
                    x_coords,
                    y_coords,
                    cum_transposed,
                    levels=levels_to_draw,        # Evaluates +c and/or -c isolines
                    color=:white,                 # High-contrast isoline color
                    linewidth=1.5,
                    linestyle=:solid,
                    contour_labels=true,          # Annotates each line with +0.5 or -0.5
                    colorbar_entry=false
                )
            end
        end

        return p

    end
    return plot(plots[1], plots[2], layout=(1, 2), size=(1600, 600))
end

end # end module