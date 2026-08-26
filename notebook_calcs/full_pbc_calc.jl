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


## Pauli matrices
sigma_x = [0 1; 1 0]
sigma_y = [0 -im; im 0]
sigma_z = [1 0; 0 -1]
identity = [1 0; 0 1]

# 1. Unified Block Generator
function real_space_perturbed_qwz_blocks(; 
    A::Real=1.0, 
    B::Real=1.0, 
    m::Real=0.0, 
    gamma::Real=0.0, 
    perturbation_type::Symbol=:none
)
    # Base QWZ blocks
    onsite = ComplexF64.((m + 2.0 * B) * sigma_z)
    tx = ComplexF64.(-0.5 * B * sigma_z - 0.5im * A * sigma_x)
    ty = ComplexF64.(-0.5 * B * sigma_z - 0.5im * A * sigma_y)

    # Apply perturbations
    if perturbation_type == :symmetric
        # gamma * (cos(kx) + cos(ky)) * I
        tx += ComplexF64.(0.5 * gamma * identity)
        ty += ComplexF64.(0.5 * gamma * identity)
    elseif perturbation_type == :tilt
        # gamma * sin(kx) * I
        tx += ComplexF64.(0.5im * gamma * identity)
        # ty is unperturbed by the 1D tilt
    elseif perturbation_type == :none
        # Base QWZ, do nothing
    else
        error("Unknown perturbation type: $perturbation_type")
    end

    return onsite, tx, ty
end

# Map a site and orbital to a matrix index in the 2LxLy basis.
site_index_qwz(x::Int, y::Int, orb::Int, Lx::Int, Ly::Int) = 2 * ((y - 1) * Lx + (x - 1)) + orb


# 2. Updated Real-Space Hamiltonian Builder
function real_space_perturbed_hamiltonian_qwz(
    Lx::Int,
    Ly::Int;
    A::Real=1.0,
    B::Real=1.0,
    m::Real=0.0,
    gamma::Real=0.0,
    perturbation_type::Symbol=:none,
    periodic_x::Bool=false,
    periodic_y::Bool=false,
    sparse_output::Bool=true
)
    Lx > 0 || error("Lx must be positive.")
    Ly > 0 || error("Ly must be positive.")

    # Get the perturbed blocks
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
    kx::Real,
    ky::Real;
    A::Real=1.0,
    B::Real=1.0,
    m::Real=0.0,
    gamma::Real=0.0,
    perturbation_type::Symbol=:none,
    winding_number::Int=0
)::Matrix{ComplexF64}

    winding_number >= 0 || error("winding_number must be non-negative.")

    # Base QWZ Hamiltonian in k-space
    d_x = A * sin(kx)
    d_y = A * sin(ky)
    d_z = m + 2B - B * cos(kx) - B * cos(ky)
    scalar_term = 0.0

    # Apply scalar perturbations independently of the winding distortion.
    if perturbation_type == :symmetric
        scalar_term += gamma * (cos(kx) + cos(ky))
    elseif perturbation_type == :tilt
        scalar_term += gamma * sin(kx)
    elseif perturbation_type == :none
        # Do nothing
    else
        error("Unknown perturbation type: $perturbation_type")
    end

    if winding_number > 0
        n = winding_number
        cos_nkx = cos(n * kx)
        sin_nkx = sin(n * kx)

        d_x, d_y = (
            cos_nkx * d_x - sin_nkx * d_y,
            cos_nkx * d_y + sin_nkx * d_x,
        )
    end

    H_k = scalar_term * identity + d_x * sigma_x + d_y * sigma_y + d_z * sigma_z

    return ComplexF64.(H_k)
end

function compute_bulk_band_berry_data(
    ;
    Nkx::Int=101,
    Nky::Int=101,
    kwargs...
)
    kx_vals = collect(range(-π, π; length=Nkx + 1))[1:end-1]
    ky_vals = collect(range(-π, π; length=Nky + 1))[1:end-1]
    nbands = 2

    energies = zeros(Float64, nbands, Nkx, Nky)
    eigenvectors = Array{ComplexF64}(undef, 2, nbands, Nkx, Nky)

    # Filter kwargs to only include parameters accepted by k_space_perturbed_qwz_hamiltonian
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
                energies[band, ix, iy] +
                energies[band, ix_next, iy] +
                energies[band, ix_next, iy_next] +
                energies[band, ix, iy_next]
            )
        end
    end

    band_chern_numbers = vec(sum(berry_curvature, dims=(2, 3)) ./ (2π))

    # --- Per-Band Energy Accumulation (3D Spatial Field) ---
    cum_chern_per_band = zeros(Float64, nbands, Nkx, Nky)
    for band in 1:nbands
        flat_E = vec(plaquette_energies[band, :, :])
        flat_F = vec(berry_curvature[band, :, :]) ./ (2π)
        
        order = sortperm(flat_E)
        cum_sorted = cumsum(flat_F[order])
        
        # Map sorted 1D sum back onto 2D spatial grid for this band
        cum_band_2d = zeros(Float64, Nkx, Nky)
        cum_band_2d[order] .= cum_sorted
        cum_chern_per_band[band, :, :] .= cum_band_2d
    end

    # --- Global System Accumulation (1D Vectors across all bands) ---
    flat_energies_all = vec(plaquette_energies)
    flat_weights_all = vec(berry_curvature) ./ (2π)
    global_order = sortperm(flat_energies_all)

    accumulation_energies = flat_energies_all[global_order]
    accumulation_weights = flat_weights_all[global_order]
    cumulative_chern = cumsum(flat_weights_all[global_order])

    return (
        kx_vals=kx_vals,
        ky_vals=ky_vals,
        energies=energies,
        berry_curvature=berry_curvature,
        plaquette_energies=plaquette_energies,
        band_chern_numbers=band_chern_numbers,
        cum_chern_per_band=cum_chern_per_band,
        accumulation_energies=accumulation_energies,
        accumulation_weights=accumulation_weights,
        cumulative_chern=cumulative_chern,
    )
end

function edge_state_prediction(

)
    """
    predicts 
    (1) at what (kx, ky) and E the edge states cross within the bulk gap
        given by E = ± sqrt(A^2 * (sin(kx)^2 + sin(ky)^2) + (m + 2B - B * cos(kx) - B * cos(ky))^2)
    (2) the extent of the edge states in k
        given by cos(kx) + cos(ky) = (m + 2B) / B
    """

end

############################################################################

function plt_bandstructure_heatmap(
    kx_vals::Vector{Float64},
    ky_vals::Vector{Float64},
    energies::Array{Float64, 3};
    title::LaTeXString="",
    xlabel::LaTeXString=L"k_x",
    ylabel::LaTeXString=L"k_y",
    colour=:curl,
    share_colour_scale::Bool=true,
)
    """
    Plots two heatmaps of the bandstructure side-by-side for bands 1 and 2,
    with an optional shared colour scale.
    """
    c_limits = share_colour_scale ? extrema(energies[1:2, :, :]) : :auto

    plt_band1 = heatmap(
        kx_vals,
        ky_vals,
        energies[1, :, :],
        xlabel=xlabel,
        ylabel=ylabel,
        title=isempty(title) ? "Band 1" : "$title (Band 1)",
        color=colour,
        clims=c_limits,
        aspect_ratio=:equal
    )

    plt_band2 = heatmap(
        kx_vals,
        ky_vals,
        energies[2, :, :],
        xlabel=xlabel,
        ylabel=ylabel,
        title=isempty(title) ? "Band 2" : "$title (Band 2)",
        color=colour,
        clims=c_limits,
        aspect_ratio=:equal
    )

    return plot(plt_band1, plt_band2, layout=(1, 2), size=(800, 400))
end

function plt_k_resolved_F_xy_heatmaps(
    kx_vals::Vector{Float64},
    ky_vals::Vector{Float64},
    berry_curvature::Array{Float64, 3};
    title::LaTeXString=L"",
    xlabel::LaTeXString=L"k_x",
    ylabel::LaTeXString=L"k_y",
    colour=:RdBu,
    shift_to_centers::Bool=false,
    share_colour_scale::Bool=true,
)
    """
    Plots the local k-resolved contribution to the Chern number (F_{xy} / 2π)
    for both bands side-by-side in a single figure.
    """
    # Convert phase angle into Chern number density per plaquette
    chern_density = berry_curvature ./ (2π)

    # Grid coordinate setup
    dkx = kx_vals[2] - kx_vals[1]
    dky = ky_vals[2] - ky_vals[1]
    x_coords = shift_to_centers ? (kx_vals .+ dkx / 2) : kx_vals
    y_coords = shift_to_centers ? (ky_vals .+ dky / 2) : ky_vals

    # Unified color limit calculation if shared
    global_max = share_colour_scale ? maximum(abs, chern_density[1:2, :, :]) : 0.0

    plots = map(1:2) do band_index
        density_b = chern_density[band_index, :, :]
        c_limits = share_colour_scale ? 
            (-global_max, global_max) : 
            let max_val = maximum(abs, density_b)
                (-max_val, max_val)
            end

        total_chern = round(sum(density_b), digits=4)
        sub_title = isempty(title) ? 
            L"F_{xy} / 2\pi (Band %$(band_index), C = %$(total_chern))" : 
            L"%$(title) (Band %$(band_index), C = %$(total_chern))"

        heatmap(
            x_coords,
            y_coords,
            density_b',
            xlabel=xlabel,
            ylabel=ylabel,
            title=sub_title,
            color=colour,
            clims=c_limits,
            aspect_ratio=:equal
        )
    end

    return plot(plots[1], plots[2], layout=(1, 2), size=(900, 400))
end

# function plt_berry_curvature_heatmaps(
#     kx_vals::Vector{Float64},
#     ky_vals::Vector{Float64},
#     berry_curvature::Array{Float64, 3},
#     plaquette_energies::Array{Float64, 3};
#     title::LaTeXString=L"",
#     xlabel::LaTeXString=L"k_x",
#     ylabel::LaTeXString=L"k_y",
#     colour=:RdBu,
#     shift_to_centers::Bool=false,
#     share_colour_scale::Bool=true,
#     add_contours::Bool=true,
# )
#     """
#     Plots the local k-resolved contribution to the Chern number (F_{xy} / 2π)
#     for both bands side-by-side, with optional contours showing energy-cumulative Chern values.
#     """
#     chern_density = berry_curvature ./ (2π)

#     dkx = kx_vals[2] - kx_vals[1]
#     dky = ky_vals[2] - ky_vals[1]
#     x_coords = shift_to_centers ? (kx_vals .+ dkx / 2) : kx_vals
#     y_coords = shift_to_centers ? (ky_vals .+ dky / 2) : ky_vals

#     global_max = share_colour_scale ? maximum(abs, chern_density[1:2, :, :]) : 0.0

#     plots = map(1:2) do band_index
#         density_b = chern_density[band_index, :, :]
#         energies_b = plaquette_energies[band_index, :, :]

#         # 1. Compute Cumulative Chern Field as a function of local energy E(kx, ky)
#         flat_E = vec(energies_b)
#         flat_F = vec(density_b)
#         order = sortperm(flat_E)
        
#         sorted_E = flat_E[order]
#         cum_chern_sorted = cumsum(flat_F[order])
        
#         # Map sorted cumulative values back to original 2D (kx, ky) grid layout
#         cum_chern_2d = zeros(Float64, size(energies_b))
#         cum_chern_2d[order] .= cum_chern_sorted

#         # 2. Setup Subplot Properties
#         c_limits = share_colour_scale ? 
#             (-global_max, global_max) : 
#             let max_val = maximum(abs, density_b)
#                 (-max_val, max_val)
#             end

#         total_chern = round(sum(density_b), digits=4)
#         sub_title = isempty(title) ? 
#             L"F_{xy} / 2\pi \text{ (Band %$(band_index), C = %$(total_chern))}" : 
#             L"%$(title) (Band %$(band_index), C = %$(total_chern))"

#         # 3. Base Heatmap Plot
#         p = heatmap(
#             x_coords,
#             y_coords,
#             density_b',
#             xlabel=xlabel,
#             ylabel=ylabel,
#             title=sub_title,
#             color=colour,
#             clims=c_limits,
#             aspect_ratio=:equal
#         )

#         # 4. Overlay Energy-Accumulated Chern Contours
#         if add_contours
#             contour!(
#                 p,
#                 x_coords,
#                 y_coords,
#                 cum_chern_2d',
#                 levels=-1.0:0.1:1.0,
#                 color=:black,
#                 linestyle=:dash,
#                 linewidth=0.8,
#                 contour_labels=true,
#                 colorbar_entry=false
#             )
#         end

#         return p
#     end

#     return plot(plots[1], plots[2], layout=(1, 2), size=(900, 400))
# end

# function plt_berry_curvature_heatmaps(
#     kx_vals::Vector{Float64},
#     ky_vals::Vector{Float64},
#     berry_curvature::Array{Float64, 3},
#     plaquette_energies::Array{Float64, 3};
#     title::LaTeXString=L"",
#     xlabel::LaTeXString=L"k_x",
#     ylabel::LaTeXString=L"k_y",
#     colour=:RdBu,
#     shift_to_centers::Bool=false,
#     share_colour_scale::Bool=true,
#     add_contours::Bool=true,
#     n_contour_levels::Int=6
# )
#     """
#     Plots the local k-resolved contribution to the Chern number (F_{xy} / 2π)
#     for both bands side-by-side, with band-specific contours showing energy-cumulative Chern values.
#     """
#     chern_density = berry_curvature ./ (2π)

#     dkx = kx_vals[2] - kx_vals[1]
#     dky = ky_vals[2] - ky_vals[1]
#     x_coords = shift_to_centers ? (kx_vals .+ dkx / 2) : kx_vals
#     y_coords = shift_to_centers ? (ky_vals .+ dky / 2) : ky_vals

#     global_max = share_colour_scale ? maximum(abs, chern_density[1:2, :, :]) : 0.0

#     plots = map(1:2) do band_index
#         density_b = chern_density[band_index, :, :]
#         energies_b = plaquette_energies[band_index, :, :]

#         # 1. Prepare Transposed Arrays for Plotting Alignment
#         # Plots.jl expects matrix z such that z[iy, ix] maps to (x[ix], y[iy])
#         density_transposed = density_b'
#         energies_transposed = energies_b'

#         # 2. Compute Cumulative Chern Field on Transposed Grid
#         flat_E = vec(energies_transposed)
#         flat_F = vec(density_transposed)
#         order = sortperm(flat_E)
        
#         cum_chern_sorted = cumsum(flat_F[order])
        
#         # Map sorted cumulative values directly into the transposed 2D grid layout
#         cum_chern_2d_transposed = zeros(Float64, size(energies_transposed))
#         cum_chern_2d_transposed[order] .= cum_chern_sorted

#         # 3. Setup Subplot Properties
#         c_limits = share_colour_scale ? 
#             (-global_max, global_max) : 
#             let max_val = maximum(abs, density_b)
#                 (-max_val, max_val)
#             end

#         total_chern = round(sum(density_b), digits=4)
#         sub_title = isempty(title) ? 
#             L"F_{xy} / 2\pi \text{ (Band %$(band_index), C = %$(total_chern))}" : 
#             L"%$(title) \text{ (Band %$(band_index), C = %$(total_chern))}"

#         # 4. Base Heatmap Plot
#         p = heatmap(
#             x_coords,
#             y_coords,
#             density_transposed,
#             xlabel=xlabel,
#             ylabel=ylabel,
#             title=sub_title,
#             color=colour,
#             clims=c_limits,
#             aspect_ratio=:equal
#         )

#         # 5. Overlay Energy-Accumulated Chern Contours
#         if add_contours
#             # Ensure levels are sorted in ascending order for Plots.jl / GR
#             raw_levels = range(0.0, total_chern, length=n_contour_levels)
#             contour_levels = sort(collect(raw_levels))

#             contour!(
#                 p,
#                 x_coords,
#                 y_coords,
#                 cum_chern_2d_transposed, # Already transposed to match density_transposed!
#                 levels=contour_levels,
#                 color=:black,
#                 linestyle=:dash,
#                 linewidth=0.8,
#                 contour_labels=true,
#                 colorbar_entry=false
#             )
#         end

#         return p
#     end

#     return plot(plots[1], plots[2], layout=(1, 2), size=(900, 400))
# end


function plt_berry_curvature_heatmaps(
    kx_vals::Vector{Float64},
    ky_vals::Vector{Float64},
    berry_curvature::Array{Float64, 3},
    cum_chern_per_band::Array{Float64, 3}; # Directly accept pre-computed 3D field
    title::LaTeXString=L"",
    xlabel::LaTeXString=L"k_x",
    ylabel::LaTeXString=L"k_y",
    colour=:RdBu,
    shift_to_centers::Bool=false,
    share_colour_scale::Bool=true,
    add_contours::Bool=true,
    n_contour_levels::Int=6
)
    chern_density = berry_curvature ./ (2π)

    dkx = kx_vals[2] - kx_vals[1]
    dky = ky_vals[2] - ky_vals[1]
    x_coords = shift_to_centers ? (kx_vals .+ dkx / 2) : kx_vals
    y_coords = shift_to_centers ? (ky_vals .+ dky / 2) : ky_vals

    global_max = share_colour_scale ? maximum(abs, chern_density[1:2, :, :]) : 0.0

    plots = map(1:2) do band_index
        density_transposed = chern_density[band_index, :, :]'
        cum_chern_transposed = cum_chern_per_band[band_index, :, :]' # Transposed to match

        c_limits = share_colour_scale ? 
            (-global_max, global_max) : 
            let max_val = maximum(abs, density_transposed)
                (-max_val, max_val)
            end

        total_chern = round(sum(density_transposed), digits=4)
        sub_title = isempty(title) ? 
            L"F_{xy} / 2\pi (Band %$(band_index), C = %$(total_chern))" : 
            L"%$(title) (Band %$(band_index), C = %$(total_chern))"

        p = heatmap(
            x_coords, y_coords, density_transposed,
            xlabel=xlabel, ylabel=ylabel, title=sub_title,
            color=colour, clims=c_limits, aspect_ratio=:equal
        )

        if add_contours
            raw_levels = range(0.0, total_chern, length=n_contour_levels)
            contour_levels = sort(collect(raw_levels))

            contour!(
                p, x_coords, y_coords, cum_chern_transposed,
                levels=contour_levels, color=:black, linestyle=:dash,
                linewidth=0.8, contour_labels=true, colorbar_entry=false
            )
        end

        return p
    end

    return plot(plots[1], plots[2], layout=(1, 2), size=(900, 400))
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
    """
    Plots the 2D energy-accumulated Chern number C(kx, ky) directly as a 
    heatmap in k-space for both bands side-by-side.
    """
    dkx = kx_vals[2] - kx_vals[1]
    dky = ky_vals[2] - ky_vals[1]
    x_coords = shift_to_centers ? (kx_vals .+ dkx / 2) : kx_vals
    y_coords = shift_to_centers ? (ky_vals .+ dky / 2) : ky_vals

    plots = map(1:2) do band_index
        # Transpose matrix to match Plots.jl (y=row, x=col) alignment: z[iy, ix]
        cum_transposed = cum_chern_per_band[band_index, :, :]'

        min_val = minimum(cum_transposed)
        max_val = maximum(cum_transposed)
        
        # Determine total integrated Chern value at peak energy
        total_chern = abs(max_val) > abs(min_val) ? round(max_val, digits=4) : round(min_val, digits=4)

        sub_title = isempty(title) ? 
            L"\mathcal{C}(k_x, k_y) (Band %$(band_index), C = %$(total_chern))" : 
            L"%$(title) (Band %$(band_index), C = %$(total_chern))"

        p = heatmap(
            x_coords,
            y_coords,
            cum_transposed,
            xlabel=xlabel,
            ylabel=ylabel,
            title=sub_title,
            color=colour,
            clims=(min_val, max_val),
            aspect_ratio=:equal,
            colorbar_title=L"\mathcal{C}(k_x, k_y)"
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

    return plot(plots[1], plots[2], layout=(1, 2), size=(950, 400))
end


############################################################################

 function unpack_saved_data(
    folder_path::String
)
    """
    Unpacks all saved data files in the specified folder into a dictionary.
    """
    data_dict = Dict{String, Any}()

    for file in readdir(folder_path)
        if endswith(file, ".jld2")
            file_path = joinpath(folder_path, file)
            @load file_path data
            data_dict[file] = data
        end
    end

    return data_dict
end


function debug_accumulation_extrema(data)
    """
    Prints debugging info identifying where in k-space and at what energies
    the accumulated Chern numbers reach their extrema for each band.
    """
    nbands = size(data.cum_chern_per_band, 1)

    println("="^60)
    println("ACCUMULATED CHERN EXTREMA DEBUG INFO")
    println("="^60)

    for band in 1:nbands
        cum_field = data.cum_chern_per_band[band, :, :]
        energies_field = data.plaquette_energies[band, :, :]
        total_chern = data.band_chern_numbers[band]

        # Determine target extremum (max for positive Chern, min for negative)
        is_pos = total_chern >= 0
        target_val = is_pos ? maximum(cum_field) : minimum(cum_field)
        target_idx = is_pos ? argmax(cum_field) : argmin(cum_field)

        ix, iy = target_idx[1], target_idx[2]
        kx = data.kx_vals[ix]
        ky = data.ky_vals[iy]
        E_at_target = energies_field[ix, iy]

        # Minimum energy state in band (where accumulation starts)
        min_E_idx = argmin(energies_field)
        min_E_kx = data.kx_vals[min_E_idx[1]]
        min_E_ky = data.ky_vals[min_E_idx[2]]
        E_start = energies_field[min_E_idx]

        # Maximum energy state in band
        max_E_idx = argmax(energies_field)
        max_E_kx = data.kx_vals[max_E_idx[1]]
        max_E_ky = data.ky_vals[max_E_idx[2]]
        E_end = energies_field[max_E_idx]

        println("BAND $band:")
        println("  Total Band Chern Number C_$band : $total_chern")
        println("  Energy Range [E_min, E_max]   : [$(round(E_start, digits=4)), $(round(E_end, digits=4))]")
        println("    -> E_min at (kx, ky)        : ($(round(min_E_kx, digits=4)), $(round(min_E_ky, digits=4)))")
        println("    -> E_max at (kx, ky)        : ($(round(max_E_kx, digits=4)), $(round(max_E_ky, digits=4)))")
        println("  Accumulated Value Peak        : $(round(target_val, digits=5))")
        println("    -> Grid Index (ix, iy)     : ($ix, $iy)")
        println("    -> Coordinates (kx, ky)     : ($(round(kx, digits=4)), $(round(ky, digits=4)))")
        println("    -> Local Plaquette Energy   : $(round(E_at_target, digits=4))")
        println("-"^60)
    end
end


############################################################################
"""
Computing the PBS bandstructure and its berry curvature for distorted QWZ model
"""

As = [1.0]
Bs = [1.0]
ms = [-1.0]
gammas = collect(0.0:0.1:3.0)
distortion_type = :symmetric
Nkx = 101
Nky = 101

data_folder = joinpath("data", "bulk_band_berry_data")
plot_folder = joinpath("plots", "bulk_band_berry_data")
isdir(data_folder) || mkpath(data_folder)
isdir(plot_folder) || mkpath(plot_folder)

generate_new_data = true
generate_new_plots = true

if generate_new_data
    println("Generating new data...")
    @showprogress for (A, B, m, gamma) in Iterators.product(As, Bs, ms, gammas)
        data = compute_bulk_band_berry_data(
            Nkx=101,
            Nky=101,
            A=A,
            B=B,
            m=m,
            gamma=gamma,
            perturbation_type=distortion_type
        )

        save_path = joinpath(data_folder, "bandstrucutre_data")
        isdir(save_path) || mkpath(save_path)
        @save joinpath(save_path, "bulk_band_berry_data_A$(A)_B$(B)_m$(m)_gamma$(gamma)_perturbation_$(distortion_type).jld2") data

        # debug_accumulation_extrema(data)
    end
else
    println("Using existing data...")
end



if generate_new_plots
    data_to_unpack = joinpath(data_folder, "bandstrucutre_data")
    println("Plotting data from $(data_to_unpack)...")

    @showprogress for (A, B, m, gamma) in Iterators.product(As, Bs, ms, gammas)
        # 1. Build the target filename directly
        filename = "bulk_band_berry_data_A$(A)_B$(B)_m$(m)_gamma$(gamma)_perturbation_$(distortion_type).jld2"
        filepath = joinpath(data_to_unpack, filename)
        
        # 2. Skip if file is missing
        if !isfile(filepath)
            @warn "File not found: $filepath, skipping..."
            continue
        end

        # 3. Load the NamedTuple directly
        @load filepath data

        # 4. Bare bandstructure heatmaps for both bands
        bandstruct_heatmaps = plt_bandstructure_heatmap(
            data.kx_vals,
            data.ky_vals,
            data.energies;
            title = L"\gamma = %$(gamma)",
            colour = :curl,
            share_colour_scale = false
        )
        
        savefig(
            bandstruct_heatmaps, 
            joinpath(plot_folder, "bulk_bandstruct_A$(A)_B$(B)_m$(m)_gamma$(gamma)_perturbation_$(distortion_type).png")
        )

        # 5. Plot the Berry curvature heatmaps. F_xy flux per plaquette resolved for (kx, ky) for both bands.
        F_xy_heatmaps = plt_k_resolved_F_xy_heatmaps(
            data.kx_vals,
            data.ky_vals,
            data.berry_curvature;
            title = L"\gamma = %$(gamma), F_{xy} / 2\pi",
            # colour = :RdBu,
            # share_colour_scale = false
        )

        savefig(
            F_xy_heatmaps, 
            joinpath(plot_folder, "bulk_band_dEdC_A$(A)_B$(B)_m$(m)_gamma$(gamma)_perturbation_$(distortion_type).png")
        )

        # # 6. Plot the same as 6. but including contours of the energy-accumulated Chern number 
        # # BROKEN CONTOURS
        # chern_heatmaps = plt_berry_curvature_heatmaps(
        #     data.kx_vals,
        #     data.ky_vals,
        #     data.berry_curvature,
        #     data.plaquette_energies;
        #     title = L"\gamma = %$(gamma), F_{xy} / 2\pi",
        #     add_contours = true
        # )

        # savefig(
        #     chern_heatmaps, 
        #     joinpath(plot_folder, "bulk_band_dEdC_contours_A$(A)_B$(B)_m$(m)_gamma$(gamma)_perturbation_$(distortion_type).png")
        # )

        # 7. Plot the accumulated chern number as a heatmap on (kx, ky) grid
        cum_chern_heatmaps = plt_accumulated_chern_heatmaps(
            data.kx_vals,
            data.ky_vals,
            data.cum_chern_per_band;
            title = L"\gamma = %$(gamma), \mathcal{C}(k_x, k_y)",
            colour = :viridis
        )

        savefig(
            cum_chern_heatmaps, 
            joinpath(plot_folder, "bulk_band_accumulated_chern_A$(A)_B$(B)_m$(m)_gamma$(gamma)_perturbation_$(distortion_type).png")
        )
    end
end