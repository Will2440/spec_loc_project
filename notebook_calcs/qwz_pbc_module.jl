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

export compute_bulk_band_berry_data, plt_bandstructure_heatmap, plt_bandstructure_fermi_surface_heatmap
export plt_k_resolved_F_xy_heatmaps, plt_accumulated_chern_heatmaps
export real_space_perturbed_qwz_blocks, real_space_perturbed_hamiltonian_qwz
export k_space_perturbed_qwz_hamiltonian, unpack_saved_data, debug_accumulation_extrema
export compute_euler_characteristic, compute_pockets_and_holes, get_energy_at_chern, get_energy_at_chern_2d


## Pauli matrices
const sigma_x = [0 1; 1 0]
const sigma_y = [0 -im; im 0]
const sigma_z = [1 0; 0 -1]
const identity = [1 0; 0 1]

pi_ticks = (
    [-π, -3π/4, -π/2, -π/4, 0, π/4, π/2, 3π/4, π],
    [L"-\pi", L"-3\pi/4", L"-\pi/2", L"-\pi/4", L"0", L"\pi/4", L"\pi/2", L"3\pi/4", L"\pi"]
)

# EMBEDDING IMPLEMENTATION TOGGLE (hardcoded for testing)
# Options: :sigma_z_rotation (current), :orbital_displacement_phase (new)
const EMBEDDING_TYPE = :sigma_z_rotation

# # 1. Unified Block Generator
# function real_space_perturbed_qwz_blocks(; 
#     A::Real=1.0, B::Real=1.0, m::Real=0.0, gamma::Real=0.0, perturbation_type::Symbol=:none
# )
#     onsite = ComplexF64.((m + 2.0 * B) * sigma_z)
#     tx = ComplexF64.(-0.5 * B * sigma_z - 0.5im * A * sigma_x)
#     ty = ComplexF64.(-0.5 * B * sigma_z - 0.5im * A * sigma_y)

#     # Apply generic scalar perturbation to real-space blocks
#     # (Real-space representation is approximate for k-space-dependent perturbations)
#     if perturbation_type in (:sym_cos_sum, :sym_cos_diff, :sym_cos_add, :sym_cos_sub)
#         tx += ComplexF64.(0.5 * gamma * identity)
#         ty += ComplexF64.(0.5 * gamma * identity)
#     elseif perturbation_type in (:asym_sin_sum, :asym_sin_diff, :asym_sin_add, :asym_sin_sub)
#         tx += ComplexF64.(0.5im * gamma * identity)
#         ty += ComplexF64.(0.5im * gamma * identity)
#     elseif perturbation_type != :none
#         error("Unknown perturbation type: $perturbation_type")
#     end

#     return onsite, tx, ty
# end

# site_index_qwz(x::Int, y::Int, orb::Int, Lx::Int, Ly::Int) = 2 * ((y - 1) * Lx + (x - 1)) + orb

# # 2. Updated Real-Space Hamiltonian Builder
# function real_space_perturbed_hamiltonian_qwz(
#     Lx::Int, Ly::Int; A::Real=1.0, B::Real=1.0, m::Real=0.0, gamma::Real=0.0,
#     perturbation_type::Symbol=:none, periodic_x::Bool=false, periodic_y::Bool=false, sparse_output::Bool=true
# )
#     Lx > 0 || error("Lx must be positive.")
#     Ly > 0 || error("Ly must be positive.")

#     onsite, tx, ty = real_space_perturbed_qwz_blocks(; A=A, B=B, m=m, gamma=gamma, perturbation_type=perturbation_type)
#     nsites = Lx * Ly
#     dim = 2 * nsites
#     H = sparse_output ? spzeros(ComplexF64, dim, dim) : zeros(ComplexF64, dim, dim)

#     function add_block!(mat, row_site::Tuple{Int, Int}, col_site::Tuple{Int, Int}, block::AbstractMatrix{<:Number})
#         (xr, yr) = row_site
#         (xc, yc) = col_site
#         row_base = site_index_qwz(xr, yr, 1, Lx, Ly)
#         col_base = site_index_qwz(xc, yc, 1, Lx, Ly)
#         @inbounds for a in 0:1, b in 0:1
#             mat[row_base + a, col_base + b] += ComplexF64(block[a + 1, b + 1])
#         end
#         return nothing
#     end

#     for y in 1:Ly, x in 1:Lx
#         add_block!(H, (x, y), (x, y), onsite)

#         if x < Lx
#             add_block!(H, (x + 1, y), (x, y), tx)
#             add_block!(H, (x, y), (x + 1, y), tx')
#         elseif periodic_x
#             add_block!(H, (1, y), (x, y), tx)
#             add_block!(H, (x, y), (1, y), tx')
#         end

#         if y < Ly
#             add_block!(H, (x, y + 1), (x, y), ty)
#             add_block!(H, (x, y), (x, y + 1), ty')
#         elseif periodic_y
#             add_block!(H, (x, 1), (x, y), ty)
#             add_block!(H, (x, y), (x, 1), ty')
#         end
#     end

#     return H
# end

function k_space_perturbed_qwz_hamiltonian(
    kx::Real, ky::Real; A::Real=1.0, B::Real=1.0, m::Real=0.0, gamma::Real=0.0,
    perturbation_type::Symbol=:none, winding_number::Int=0, phi::Real=0.0, orbital_displacement::Real=0.0
)::Matrix{ComplexF64}

    winding_number >= 0 || error("winding_number must be non-negative.")

    d_x = A * sin(kx)
    d_y = A * sin(ky)
    d_z = m + 2B - B * cos(kx) - B * cos(ky)
    scalar_term = 0.0

    # Symmetric perturbations (cos variants)
    if perturbation_type == :sym_cos_sum
        scalar_term += gamma * (cos(kx) + cos(ky))
    elseif perturbation_type == :sym_cos_diff
        scalar_term += gamma * (cos(kx) - cos(ky))
    elseif perturbation_type == :sym_cos_add
        scalar_term += gamma * cos(kx + ky)
    elseif perturbation_type == :sym_cos_sub
        scalar_term += gamma * cos(kx - ky)
    # Antisymmetric perturbations (sin variants)
    elseif perturbation_type == :asym_sin_sum
        scalar_term += gamma * (sin(kx) + sin(ky))
    elseif perturbation_type == :asym_sin_diff
        scalar_term += gamma * (sin(kx) - sin(ky))
    elseif perturbation_type == :asym_sin_add
        scalar_term += gamma * sin(kx + ky)
    elseif perturbation_type == :asym_sin_sub
        scalar_term += gamma * sin(kx - ky)
    elseif perturbation_type != :none
        error("Unknown perturbation type: $perturbation_type")
    end

    if winding_number > 0
        n = winding_number
        cos_nkx = cos(n * kx)
        sin_nkx = sin(n * kx)
        d_x, d_y = (cos_nkx * d_x - sin_nkx * d_y, cos_nkx * d_y + sin_nkx * d_x)
    end

    # Build Hamiltonian before embedding
    H_k = scalar_term * identity + d_x * sigma_x + d_y * sigma_y + d_z * sigma_z

    # Apply orbital embedding based on selected implementation
    if abs(orbital_displacement) > 1e-14
        dx_phi = orbital_displacement * cos(phi)
        dy_phi = orbital_displacement * sin(phi)
        
        if EMBEDDING_TYPE == :sigma_z_rotation
            # ============================================
            # IMPLEMENTATION 1: σ_z basis rotation
            # Applies: U(k) = exp(-i*θ*σ_z) where θ = k·d
            # Rotates d-vector in (d_x, d_y) plane by angle θ
            # ============================================
            theta = kx * dx_phi + ky * dy_phi
            cos_theta = cos(theta)
            sin_theta = sin(theta)
            
            # Extract current d-vector components
            d_x_current = d_x
            d_y_current = d_y
            
            # Apply rotation: (d_x, d_y) -> (d_x*cos(θ) - d_y*sin(θ), d_x*sin(θ) + d_y*cos(θ))
            d_x_rot = cos_theta * d_x_current - sin_theta * d_y_current
            d_y_rot = sin_theta * d_x_current + cos_theta * d_y_current
            
            # Rebuild with rotated components
            H_k = scalar_term * identity + d_x_rot * sigma_x + d_y_rot * sigma_y + d_z * sigma_z
            
        elseif EMBEDDING_TYPE == :orbital_displacement_phase
            # ============================================
            # IMPLEMENTATION 2: Orbital displacement phase
            # Applies phase factor exp(i*k·d) to encode 
            # physical orbital position difference
            # ============================================
            phase = kx * dx_phi + ky * dy_phi
            phase_factor = exp(im * phase)
            
            # Apply phase to entire Hamiltonian
            # This encodes the orbital displacement through k-dependent gauge
            H_k = phase_factor * H_k
            
        else
            error("Unknown EMBEDDING_TYPE: $EMBEDDING_TYPE")
        end
    end

    return ComplexF64.(H_k)
end

function compute_bulk_band_berry_data(; Nkx::Int=101, Nky::Int=101, kwargs...)
    kx_vals = collect(range(-π, π; length=Nkx + 1))[1:end-1]
    ky_vals = collect(range(-π, π; length=Nky + 1))[1:end-1]
    nbands = 2

    energies = zeros(Float64, nbands, Nkx, Nky)
    eigenvectors = Array{ComplexF64}(undef, 2, nbands, Nkx, Nky)

    k_space_kwargs = (; [k => v for (k, v) in pairs(kwargs) if k ∈ (:A, :B, :m, :gamma, :perturbation_type, :winding_number, :phi, :orbital_displacement)]...)

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
    # FIXED: Add transpose to match contour and other heatmap functions (Plots.jl expects z[i,j] at (x[j], y[i]))
    plt_band1 = heatmap(kx_vals, ky_vals, energies[1, :, :]', xlabel=xlabel, ylabel=ylabel, title="Bandstrucutre (Band 1 Lower)", xlims=(-pi, pi), ylims=(-pi, pi), xticks=pi_ticks, yticks=pi_ticks, color=colour, clims=c_limits, colorbar_title=L"E(k_x, k_y)", aspect_ratio=:equal)
    plt_band2 = heatmap(kx_vals, ky_vals, energies[2, :, :]', xlabel=xlabel, ylabel=ylabel, title="Bandstrucutre (Band 2 Upper)", xlims=(-pi, pi), ylims=(-pi, pi), xticks=pi_ticks, yticks=pi_ticks, color=colour, clims=c_limits, colorbar_title=L"E(k_x, k_y)", aspect_ratio=:equal)
    return plot(plt_band1, plt_band2, layout=(1, 2), size=(1600, 600))
end

function plt_bandstructure_fermi_surface_heatmap(
    kx_vals::Vector{Float64}, 
    ky_vals::Vector{Float64}, 
    energies::Array{Float64, 3}; 
    fermi_energy::Real=0.0, 
    title::LaTeXString=L"", 
    xlabel::LaTeXString=L"k_x", 
    ylabel::LaTeXString=L"k_y", 
    colour=:curl, 
    share_colour_scale::Bool=false,
    show_fermi_contour::Bool=false,
    mask_fermi_surface::Bool=false
)

    c_limits = share_colour_scale ? extrema(energies[1:2, :, :]) : :auto

    # Create masked energy arrays: values > fermi_energy become NaN
    b1_data = mask_fermi_surface ? ifelse.(energies[1, :, :] .<= fermi_energy, energies[1, :, :], NaN) : energies[1, :, :]
    b1_data = b1_data'  # FIXED: Transpose to match Plots.jl convention and contour overlay
    b2_data = mask_fermi_surface ? ifelse.(energies[2, :, :] .<= fermi_energy, energies[2, :, :], NaN) : energies[2, :, :]
    b2_data = b2_data'  # FIXED: Transpose to match Plots.jl convention and contour overlay

    # Band 1 Plot
    plt_band1 = plot(background_color_subplot=:gray80, aspect_ratio=:equal)
    heatmap!(plt_band1, kx_vals, ky_vals, b1_data, xlabel=xlabel, ylabel=ylabel, 
             title="Bandstructure (Band 1 Lower)", xlims=(-pi, pi), ylims=(-pi, pi), 
             xticks=pi_ticks, yticks=pi_ticks, color=colour, clims=c_limits, colorbar_title=L"E(k_x, k_y)")
    
    if show_fermi_contour
        contour!(plt_band1, kx_vals, ky_vals, energies[1, :, :]', levels=[fermi_energy], color=:steelblue, linewidth=2.0)
    end

    # Band 2 Plot
    plt_band2 = plot(background_color_subplot=:gray80, aspect_ratio=:equal)
    heatmap!(plt_band2, kx_vals, ky_vals, b2_data, xlabel=xlabel, ylabel=ylabel, 
             title="Bandstructure (Band 2 Upper)", xlims=(-pi, pi), ylims=(-pi, pi), 
             xticks=pi_ticks, yticks=pi_ticks, color=colour, clims=c_limits, colorbar_title=L"E(k_x, k_y)")
    
    if show_fermi_contour
        contour!(plt_band2, kx_vals, ky_vals, energies[2, :, :]', levels=[fermi_energy], color=:steelblue, linewidth=2.0)
    end

    return plot(plt_band1, plt_band2, layout=(1, 2), size=(1600, 600))
end

function plt_k_resolved_F_xy_heatmaps(kx_vals::Vector{Float64}, ky_vals::Vector{Float64}, berry_curvature::Array{Float64, 3}; 
    title::LaTeXString=L"", 
    xlabel::LaTeXString=L"k_x", 
    ylabel::LaTeXString=L"k_y", 
    colour=:RdBu, 
    shift_to_centers::Bool=false, 
    share_colour_scale::Bool=false,
    energies::Union{Array{Float64, 3}, Nothing}=nothing,
    fermi_energy::Real=0.0,
    show_fermi_contour::Bool=false,
    mask_fermi_surface::Bool=false
)
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
        
        # Apply fermi surface masking if requested
        display_data = density_b
        if mask_fermi_surface && !isnothing(energies)
            display_data = ifelse.(energies[band_index, :, :] .<= fermi_energy, density_b, NaN)
        end
        
        ##################################################################
        # Calculate direct sum of flux in masked region (for diagnostic purposes)
        masked_flux_sum = 0.0
        for i in eachindex(display_data)
            if !isnan(display_data[i])
                masked_flux_sum += display_data[i]
            end
        end
        masked_flux_sum = round(masked_flux_sum, digits=4)
        ##################################################################

        p = plot(background_color_subplot=:gray80, aspect_ratio=:equal)
        heatmap!(p,
            x_coords, 
            y_coords, 
            display_data', 
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
        
        ##################################################################
        # Annotate direct flux sum in center of plot
        if mask_fermi_surface
            annotate!(p, (0.5, 0.5), text(L"∑F = %$(masked_flux_sum)", 12, :center, :black))
        end
        ##################################################################

        # Add fermi contour if requested
        if show_fermi_contour && !isnothing(energies)
            contour!(p, x_coords, y_coords, energies[band_index, :, :]', levels=[fermi_energy], color=:steelblue, linewidth=2.0)
        end
        
        return p
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
    contour_value_1::Float64=0.5,
    contour_value_2::Float64=0.0,
    tol::Float64=0.01
)

    dkx = kx_vals[2] - kx_vals[1]
    dky = ky_vals[2] - ky_vals[1]
    x_coords = shift_to_centers ? (kx_vals .+ dkx / 2) : kx_vals
    y_coords = shift_to_centers ? (ky_vals .+ dky / 2) : ky_vals

    pi_ticks = (
        [-π, -π/2, 0, π/2, π],
        [L"-\pi", L"-\pi/2", L"0", L"\pi/2", L"\pi"]
    )

    plots = map(1:2) do band_index
        cum_transposed = cum_chern_per_band[band_index, :, :]'
        min_val = minimum(cum_transposed)
        max_val = maximum(cum_transposed)
        total_chern = abs(max_val) > abs(min_val) ? round(max_val, digits=4) : round(min_val, digits=4)
        sub_title = "Chern number accumulated: " * L"C(E, k_x, k_y)"
        
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
            get_levels = val -> begin
                # When val is 0.0, use tol as the offset to probe into + or - band regions
                eff_val = abs(val) < 1e-8 ? tol : abs(val)
                
                levels = Float64[]
                if min_val <= eff_val <= max_val
                    push!(levels, eff_val)
                end
                if min_val <= -eff_val <= max_val
                    push!(levels, -eff_val)
                end
                return unique(levels)
            end

            # Contour 1: Solid White Lines
            levels_1 = get_levels(contour_value_1)
            if !isempty(levels_1)
                contour!(
                    p,
                    x_coords,
                    y_coords,
                    cum_transposed,
                    levels=levels_1,
                    color=:white,
                    linewidth=1.5,
                    linestyle=:solid,
                    contour_labels=true,
                    colorbar_entry=false
                )
            end

            # Contour 2: Dashed Red Lines
            levels_2 = get_levels(contour_value_2)
            if !isempty(levels_2)
                contour!(
                    p,
                    x_coords,
                    y_coords,
                    cum_transposed,
                    levels=levels_2,
                    color=:red,
                    linewidth=1.5,
                    linestyle=:dash,
                    contour_labels=true,
                    colorbar_entry=false
                )
            end
        end

        return p
    end

    return plot(plots[1], plots[2], layout=(1, 2), size=(1600, 600))
end

# Euler Characteristic Computation Functions
function compute_euler_characteristic(
    M::AbstractMatrix{Bool}; 
    periodic::Bool=true
)
    """
    Computes the Euler characteristic χ = V - E + F on a 2D lattice mask M 
    for the occupied Fermi sea (E(k) <= E_F).
    Uses the standard topological formula: χ = V (vertices) - E (edges) + F (faces)
    """
    Nx, Ny = size(M)
    V = count(M)
    Eh, Ev, F = 0, 0, 0

    for j in 1:Ny, i in 1:Nx
        if M[i, j]
            i_next = (periodic && i == Nx) ? 1 : i + 1
            j_next = (periodic && j == Ny) ? 1 : j + 1

            # Count horizontal edges (occupied-to-unoccupied transitions in i direction)
            if j_next <= Ny && M[i, j_next]
                Eh += 1
            end
            # Count vertical edges (occupied-to-unoccupied transitions in j direction)
            if i_next <= Nx && M[i_next, j]
                Ev += 1
            end
            # Count 2x2 faces/plaquettes (fully occupied 2x2 blocks)
            if i_next <= Nx && j_next <= Ny && M[i_next, j] && M[i, j_next] && M[i_next, j_next]
                F += 1
            end
        end
    end

    return V - (Eh + Ev) + F
end

function compute_pockets_and_holes(M::AbstractMatrix{Bool}; periodic::Bool=true)
    Nx, Ny = size(M)
    chi = compute_euler_characteristic(M; periodic=periodic)
    
    visited = fill(false, Nx, Ny)
    n_pockets = 0

    # 4-neighbor direction offsets
    dirs = ((1, 0), (-1, 0), (0, 1), (0, -1))

    for j in 1:Ny, i in 1:Nx
        if M[i, j] && !visited[i, j]
            # BFS to find connected components (pockets)
            queue = [(i, j)]
            visited[i, j] = true

            while !isempty(queue)
                (ci, cj) = popfirst!(queue)
                for (di, dj) in dirs
                    ni = ci + di
                    nj = cj + dj
                    
                    if periodic
                        # Use mod1 for 1-based periodic indexing
                        ni = mod1(ni, Nx)
                        nj = mod1(nj, Ny)
                    else
                        if ni < 1 || ni > Nx || nj < 1 || nj > Ny
                            continue
                        end
                    end

                    if M[ni, nj] && !visited[ni, nj]
                        visited[ni, nj] = true
                        push!(queue, (ni, nj))
                    end
                end
            end
            
            n_pockets += 1
        end
    end

    n_holes = n_pockets - chi
    return n_pockets, n_holes
end

"""
    get_energy_at_chern(energies::Vector, chern_values::Vector, target_chern::Float64) -> Float64

    Find the energy corresponding to a given accumulated Chern value using linear interpolation.
    Returns NaN if target_chern is outside the range of chern_values.

    Arguments:
    - energies: Sorted energy values
    - chern_values: Cumulative Chern values corresponding to energies
    - target_chern: Target accumulated Chern value to find

    Returns:
    - Interpolated energy at the target Chern value, or NaN if outside range
"""
function get_energy_at_chern(energies::Vector, chern_values::Vector, target_chern::Float64)
    if isempty(energies) || isempty(chern_values)
        return NaN
    end
    
    min_chern = minimum(chern_values)
    max_chern = maximum(chern_values)
    
    # Check if target is outside range
    if target_chern < min_chern || target_chern > max_chern
        return NaN
    end
    
    # Find the two indices where chern_values brackets target_chern
    idx = searchsortedlast(chern_values, target_chern)
    
    if idx == 0
        return energies[1]
    elseif idx == length(chern_values)
        return energies[end]
    else
        # Linear interpolation between idx and idx+1
        c1, c2 = chern_values[idx], chern_values[idx + 1]
        e1, e2 = energies[idx], energies[idx + 1]
        
        # Avoid division by zero
        if abs(c2 - c1) < 1e-14
            return e1
        end
        
        # Linear interpolation: E = E1 + (target - C1) * (E2 - E1) / (C2 - C1)
        return e1 + (target_chern - c1) * (e2 - e1) / (c2 - c1)
    end
end

# Extract energy from 2D contour (robust to non-monotonic C(E) curves)
function get_energy_at_chern_2d(energies_2d::AbstractMatrix, chern_2d::AbstractMatrix, target_chern::Float64; tol::Float64=0.02)
    """
    Extract energy from 2D C(E,k) heatmap by finding contour points where C ≈ target_chern.
    
    Args:
        energies_2d: 2D matrix of energies at each k-point [Nkx, Nky]
        chern_2d: 2D matrix of accumulated Chern at each k-point [Nkx, Nky]
        target_chern: Target Chern value to find
        tol: Tolerance for contour matching (default 0.02)
    
    Returns:
        Median energy of all k-points where |C - target| < tol, or NaN if insufficient points
    """
    if isempty(energies_2d) || isempty(chern_2d)
        return NaN
    end
    
    # Find all k-points within tolerance of target Chern value
    contour_energies = Float64[]
    
    for i in eachindex(chern_2d)
        if abs(chern_2d[i] - target_chern) < tol
            push!(contour_energies, energies_2d[i])
        end
    end
    
    # Need at least 3 points to be confident about the contour
    if length(contour_energies) < 3
        return NaN
    end
    
    # Return median energy (robust to outliers)
    return median(contour_energies)
end


end # end module