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
using Peaks


## Pauli matrices
sigma_x = [0 1; 1 0]
sigma_y = [0 -im; im 0]
sigma_z = [1 0; 0 -1]
identity = [1 0; 0 1]

pi_ticks = (
                [-π, -3π/4, -π/2, -π/4, 0, π/4, π/2, 3π/4, π],
                [L"-\pi", L"-3\pi/4", L"-\pi/2", L"-\pi/4", L"0", L"\pi/4", L"\pi/2", L"3\pi/4", L"\pi"]
            )

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

function edge_state_prediction(;
    Nkx::Int=101,
    Nky::Int=101,
    A::Real=1.0,
    B::Real=1.0,
    m::Real=-1.0,
    gamma::Real=0.0,
    perturbation_type::Symbol=:none,
    kwargs...
)
    """
    predicts 
    (1) at what (kx, ky) and E the edge states cross within the bulk gap
        given by E = ± sqrt(A^2 * (sin(kx)^2 + sin(ky)^2) + (m + 2B - B * cos(kx) - B * cos(ky))^2)
    (2) the extent of the edge states in k
        given by cos(kx) + cos(ky) = (m + 2B) / B
    (3) the bandgap in (kx, ky)
    """

    kx_vals = collect(range(-π, π; length=Nkx + 1))[1:end-1]
    ky_vals = collect(range(-π, π; length=Nky + 1))[1:end-1]

    # 1. Bulk energy spectrum surface
    energies = zeros(Float64, 2, Nkx, Nky)
    dz_locus = zeros(Float64, Nkx, Nky)
    edge_crossings = Vector{NamedTuple{(:kx, :ky, :E), Tuple{Float64, Float64, Float64}}}()
    bulk_gap = zeros(Float64, Nkx, Nky)

    for (ix, kx) in enumerate(kx_vals), (iy, ky) in enumerate(ky_vals)
        # Scalar perturbation term
        scalar_term = 0.0
        if perturbation_type == :symmetric
            scalar_term = gamma * (cos(kx) + cos(ky))
        elseif perturbation_type == :tilt
            scalar_term = gamma * sin(kx)
        end

        dx = A * sin(kx)
        dy = A * sin(ky)
        dz = m + 2B - B * cos(kx) - B * cos(ky)

        # Bulk dispersion: E_± = scalar ± sqrt(dx^2 + dy^2 + dz^2)
        gap_magnitude = sqrt(dx^2 + dy^2 + dz^2)
        energies[1, ix, iy] = scalar_term - gap_magnitude # Lower band
        energies[2, ix, iy] = scalar_term + gap_magnitude # Upper band

        # Extent surface criterion dz(kx, ky) = 0 => cos(kx) + cos(ky) = (m + 2B) / B
        dz_locus[ix, iy] = dz

        # High-symmetry crossing detection (where gap closes, dx = dy = dz = 0)
        if isapprox(dx, 0.0; atol=1e-5) && isapprox(dy, 0.0; atol=1e-5) && isapprox(dz, 0.0; atol=1e-5)
            push!(edge_crossings, (kx=kx, ky=ky, E=scalar_term))
        end

        # Bulk bandgap calculation
        bulk_gap[ix, iy] = energies[2, ix, iy] - energies[1, ix, iy]
    end

    return (
        kx_vals = kx_vals,
        ky_vals = ky_vals,
        energies = energies,
        bulk_gap = bulk_gap,
        dz_locus = dz_locus,
        edge_crossings = edge_crossings,
        is_topological = (-4B < m < 0)
    )
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

function plt_edge_state_prediction(
    data; 
    title::LaTeXString=L"", 
    colour=:thermal
)
    gap_transposed = data.bulk_gap'

    # Base heatmap of the bulk gap
    p = heatmap(
        data.kx_vals,
        data.ky_vals,
        gap_transposed,
        xlabel=L"k_x",
        ylabel=L"k_y",
        title=isempty(title) ? L"\text{Bulk Gap } \Delta E(k_x, k_y)" : title,
        color=colour,
        aspect_ratio=:equal,
        colorbar_title=L"\Delta E"
    )

    # Contour dz = 0 (Extents of the edge state boundary locus)
    contour!(
        p,
        data.kx_vals,
        data.ky_vals,
        data.dz_locus',
        levels=[0.0],
        color=:white,
        linewidth=2,
        linestyle=:dash,
        label=L"d_z(k) = 0 \text{ (Edge Extent)}"
    )

    # Overlay band gap crossing points (where edge states close/cross)
    if !isempty(data.edge_crossings)
        x_pts = [pt.kx for pt in data.edge_crossings]
        y_pts = [pt.ky for pt in data.edge_crossings]

        scatter!(
            p,
            x_pts,
            y_pts,
            markersize=6,
            markercolor=:red,
            # markershape=:xcross,
            markerstrokewidth=2,
            label=L"\text{Dirac/Crossing Points}"
        )
    end

    return p
end

function plt_accumulated_chern_energy_vs_gamma(
    gamma_data_pairs::Vector{<:Tuple{Real, Any}}; 
    target_p::Float64=0.5, 
    target_q::Float64=0.0,
    tol::Float64=0.01
)
    # Buffers for Band 1 (±p, ±q)
    g_b1_pos_p, e_b1_pos_p = Float64[], Float64[]
    g_b1_neg_p, e_b1_neg_p = Float64[], Float64[]
    g_b1_pos_q, e_b1_pos_q = Float64[], Float64[]
    g_b1_neg_q, e_b1_neg_q = Float64[], Float64[]

    # Buffers for Band 2 (±p, ±q)
    g_b2_pos_p, e_b2_pos_p = Float64[], Float64[]
    g_b2_neg_p, e_b2_neg_p = Float64[], Float64[]
    g_b2_pos_q, e_b2_pos_q = Float64[], Float64[]
    g_b2_neg_q, e_b2_neg_q = Float64[], Float64[]

    q_is_zero = abs(target_q) < 1e-8

    for (gamma, band_data) in gamma_data_pairs
        g = Float64(gamma)

        # 2D field extractions per band
        cum_1 = band_data.cum_chern_per_band[1, :, :]
        E_1   = band_data.plaquette_energies[1, :, :]

        cum_2 = band_data.cum_chern_per_band[2, :, :]
        E_2   = band_data.plaquette_energies[2, :, :]

        # --- Band 1 (±p and ±q) ---
        for i in findall(v -> abs(v - target_p) <= tol, cum_1)
            push!(g_b1_pos_p, g); push!(e_b1_pos_p, E_1[i])
        end
        for i in findall(v -> abs(v + target_p) <= tol, cum_1)
            push!(g_b1_neg_p, g); push!(e_b1_neg_p, E_1[i])
        end
        for i in findall(v -> abs(v - target_q) <= tol, cum_1)
            push!(g_b1_pos_q, g); push!(e_b1_pos_q, E_1[i])
        end
        if !q_is_zero
            for i in findall(v -> abs(v + target_q) <= tol, cum_1)
                push!(g_b1_neg_q, g); push!(e_b1_neg_q, E_1[i])
            end
        end

        # --- Band 2 (±p and ±q) ---
        for i in findall(v -> abs(v - target_p) <= tol, cum_2)
            push!(g_b2_pos_p, g); push!(e_b2_pos_p, E_2[i])
        end
        for i in findall(v -> abs(v + target_p) <= tol, cum_2)
            push!(g_b2_neg_p, g); push!(e_b2_neg_p, E_2[i])
        end
        for i in findall(v -> abs(v - target_q) <= tol, cum_2)
            push!(g_b2_pos_q, g); push!(e_b2_pos_q, E_2[i])
        end
        if !q_is_zero
            for i in findall(v -> abs(v + target_q) <= tol, cum_2)
                push!(g_b2_neg_q, g); push!(e_b2_neg_q, E_2[i])
            end
        end
    end

    # Dynamic label formatting for q
    lbl_q_b1_pos = q_is_zero ? "Band 1, " * L" %$(target_q)" : "Band 1, " * L" +%$(target_q)"
    lbl_q_b2_pos = q_is_zero ? "Band 2, " * L" %$(target_q)" : "Band 2, " * L" +%$(target_q)"

    # Series metadata: (Buffers, Label, Color, Shape)
    series_configs = [
        # Band 1
        (g_b1_pos_p, e_b1_pos_p, "Band 1, " * L" +%$(target_p)", :steelblue, :circle),
        (g_b1_neg_p, e_b1_neg_p, "Band 1, " * L" -%$(target_p)", :steelblue, :diamond),
        (g_b1_pos_q, e_b1_pos_q, lbl_q_b1_pos,                  :blue,      :utriangle),
        (g_b1_neg_q, e_b1_neg_q, "Band 1, " * L" -%$(target_q)", :blue,      :dtriangle),
        
        # Band 2
        (g_b2_pos_p, e_b2_pos_p, "Band 2, " * L" +%$(target_p)", :red,    :circle),
        (g_b2_neg_p, e_b2_neg_p, "Band 2, " * L" -%$(target_p)", :red,    :diamond),
        (g_b2_pos_q, e_b2_pos_q, lbl_q_b2_pos,                  :orange, :utriangle),
        (g_b2_neg_q, e_b2_neg_q, "Band 2, " * L" -%$(target_q)", :orange, :dtriangle)
    ]

    # Title logic
    title_str = if target_p == target_q
        "Energy at Contour " * L"\mathcal{C}(E, k_x, k_y) = \pm %$(target_p)"
    elseif q_is_zero
        "Energy at Contours " * L"\mathcal{C}(E, k_x, k_y) = \pm %$(target_p), %$(target_q)"
    else
        "Energy at Contours " * L"\mathcal{C}(E, k_x, k_y) = \pm %$(target_p), \pm %$(target_q)"
    end

    # Aggregate global energy bounds for y-ticks
    all_energies = Float64[]
    for (g_pts, e_pts, _, _, _) in series_configs
        if !isempty(g_pts)
            append!(all_energies, e_pts)
        end
    end

    global_yticks = if !isempty(all_energies)
        floor(Int, minimum(all_energies)):1:ceil(Int, maximum(all_energies))
    else
        :auto
    end

    p = nothing

    # Render non-empty series
    for (g_pts, e_pts, lbl, clr, shp) in series_configs
        if !isempty(g_pts)
            if isnothing(p)
                p = scatter(
                    g_pts, e_pts,
                    xlabel = L"\gamma",
                    ylabel = L"E",
                    title = title_str,
                    label = lbl,
                    yticks = global_yticks,
                    color = clr,
                    markershape = shp,
                    markersize = 3,
                    markerstrokewidth = 0,
                    alpha = 0.6
                )
            else
                scatter!(
                    p,
                    g_pts, e_pts,
                    label = lbl,
                    color = clr,
                    markershape = shp,
                    markersize = 3,
                    markerstrokewidth = 0,
                    alpha = 0.6
                )
            end
        end
    end

    # Empty canvas fallback
    if isnothing(p)
        fallback_title = q_is_zero ? 
            "No contours found for p = ±$(target_p), q = $(target_q)" :
            "No contours found for p = ±$(target_p), q = ±$(target_q)"
        p = plot(
            xlabel = L"\gamma",
            ylabel = L"E",
            title = fallback_title
        )
    end

    return p
end



function plt_joint_band_accumulated_chern_heatmap(
    kx_vals::Vector{Float64},
    ky_vals::Vector{Float64},
    berry_curvature::Array{Float64, 3},
    plaquette_energies::Array{Float64, 3},
    E_cut::Float64,
    gamma::Float64;
    global_clims::Union{Nothing, Tuple{Float64, Float64}}=nothing,
    title::LaTeXString=L"",
    xlabel::LaTeXString=L"k_x",
    ylabel::LaTeXString=L"k_y",
    colour=:RdBu,
    shift_to_centers::Bool=false
)
    n_bands, nkx, nky = size(plaquette_energies)
    cum_2d = zeros(Float64, nkx, nky)
    
    # Sum over all bands at each k-point up to cutoff energy E_cut
    for n in 1:n_bands, i in 1:nkx, j in 1:nky
        if plaquette_energies[n, i, j] <= E_cut
            cum_2d[i, j] += berry_curvature[n, i, j] / (2π)
        end
    end

    dkx = kx_vals[2] - kx_vals[1]
    dky = ky_vals[2] - ky_vals[1]
    x_coords = shift_to_centers ? (kx_vals .+ dkx / 2) : kx_vals
    y_coords = shift_to_centers ? (ky_vals .+ dky / 2) : ky_vals

    cum_transposed = cum_2d'
    total_chern = round(sum(cum_2d), digits=4)

    E_str = @sprintf("%.3f", E_cut)
    t_str = L"C(k_x, k_y; E < %$(E_str)), \gamma=%$(gamma), C_{tot} = %$(total_chern)"

    # 1. Base scaling exponent on global_clims if passed, else fallback to current frame
    ref_max = isnothing(global_clims) ? maximum(abs, cum_transposed) : maximum(abs, global_clims)
    exponent = ref_max == 0 ? 0 : floor(Int, log10(ref_max))
    scale_factor = 10.0^(-exponent)

    # 2. Scale both the heatmap matrix and colorbar limits consistently
    scaled_data = cum_transposed .* scale_factor
    scaled_clims = isnothing(global_clims) ? 
        (-maximum(abs, scaled_data), maximum(abs, scaled_data)) : 
        (global_clims[1] * scale_factor, global_clims[2] * scale_factor)

    cbar_label = exponent == 0 ? 
        L"\mathcal{C}(k_x, k_y)" : 
        L"\mathcal{C}(k_x, k_y) \times 10^{%$(exponent)}"

    return heatmap(
        x_coords, y_coords, scaled_data, #cum_transposed,
        xlabel=xlabel, ylabel=ylabel, title=t_str,
        color=colour, aspect_ratio=:equal, colorbar_title=cbar_label, #colorbar_title=L"\mathcal{C}(k_x, k_y)",
        clims=scaled_clims,
        xlims=(minimum(x_coords), maximum(x_coords)),
        ylims=(minimum(y_coords), maximum(y_coords)),
        xticks=pi_ticks, yticks=pi_ticks,
    )
end


function plt_joint_band_accumulated_chern_energy_vs_gamma(
    gamma_data_pairs::Vector{<:Tuple{Real, Any}}; 
    target_p::Float64=0.5, 
    target_q::Float64=0.0,
    tol::Float64=0.01,
    ylims::Union{Nothing, Tuple{Float64, Float64}}=nothing
)
    # Global buffers unified across all bands (grouped by contour target)
    g_pos_p, e_pos_p = Float64[], Float64[]
    g_neg_p, e_neg_p = Float64[], Float64[]
    g_pos_q, e_pos_q = Float64[], Float64[]
    g_neg_q, e_neg_q = Float64[], Float64[]

    q_is_zero = abs(target_q) < 1e-8

    for (gamma, band_data) in gamma_data_pairs
        g = Float64(gamma)

        # 1. Flatten energy and Berry curvature across all bands
        E_flat = vec(band_data.plaquette_energies)

        # 2. Compute global energy-accumulated Chern number across all overlapping bands
        cum_flat = if hasproperty(band_data, :cum_chern_global)
            vec(band_data.cum_chern_global)
        elseif hasproperty(band_data, :berry_curvature)
            # Sort all plaquettes globally by energy
            b_flat = vec(band_data.berry_curvature) ./ (2π)
            sort_idx = sortperm(E_flat)
            
            # Accumulate along ascending energy order
            cum_sorted = cumsum(b_flat[sort_idx])
            
            # Restore to original index mapping
            cum_resorted = similar(cum_sorted)
            cum_resorted[sort_idx] = cum_sorted
            cum_resorted
        else
            # Fallback: flatten existing accumulated array
            vec(band_data.cum_chern_per_band)
        end

        # 3. Filter points matching target contours (independent of band index)
        for i in findall(v -> abs(v - target_p) <= tol, cum_flat)
            push!(g_pos_p, g); push!(e_pos_p, E_flat[i])
        end
        for i in findall(v -> abs(v + target_p) <= tol, cum_flat)
            push!(g_neg_p, g); push!(e_neg_p, E_flat[i])
        end
        for i in findall(v -> abs(v - target_q) <= tol, cum_flat)
            push!(g_pos_q, g); push!(e_pos_q, E_flat[i])
        end
        if !q_is_zero
            for i in findall(v -> abs(v + target_q) <= tol, cum_flat)
                push!(g_neg_q, g); push!(e_neg_q, E_flat[i])
            end
        end
    end

    # Dynamic label formatting
    lbl_q_pos = q_is_zero ? L"C = %$(target_q)" : L"C = +%$(target_q)"

    # Series configurations grouped by target Chern value (Color/Shape = Target Value)
    series_configs = [
        (g_pos_p, e_pos_p, L"C = +%$(target_p)", :steelblue,  :circle),
        (g_neg_p, e_neg_p, L"C = -%$(target_p)", :crimson,    :diamond),
        (g_pos_q, e_pos_q, lbl_q_pos,                      :forestgreen,:utriangle),
        (g_neg_q, e_neg_q, L"C = -%$(target_q)", :darkorange, :dtriangle)
    ]

    # Title generation
    title_str = if target_p == target_q
        "Global Energy Contours at " * L"C(E) = \pm %$(target_p)"
    elseif q_is_zero
        "Global Energy Contours at " * L"C(E) = \pm %$(target_p), %$(target_q)"
    else
        "Global Energy Contours at " * L"C(E) = \pm %$(target_p), \pm %$(target_q)"
    end

    # Axis tick bounds
    all_energies = Float64[]
    for (g_pts, e_pts, _, _, _) in series_configs
        if !isempty(g_pts)
            append!(all_energies, e_pts)
        end
    end

    global_yticks = if !isempty(all_energies)
        floor(Int, minimum(all_energies)):1:ceil(Int, maximum(all_energies))
    else
        :auto
    end

    p = nothing

    # Render unified series
    for (g_pts, e_pts, lbl, clr, shp) in series_configs
        if !isempty(g_pts)
            if isnothing(p)
                p = scatter(
                    g_pts, e_pts,
                    xlabel = L"\gamma",
                    ylabel = L"E",
                    title = title_str,
                    label = lbl,
                    yticks = global_yticks,
                    color = clr,
                    markershape = shp,
                    markersize = 3,
                    markerstrokewidth = 0,
                    alpha = 0.6
                )
            else
                scatter!(
                    p,
                    g_pts, e_pts,
                    label = lbl,
                    color = clr,
                    markershape = shp,
                    markersize = 3,
                    markerstrokewidth = 0,
                    alpha = 0.6
                )
            end
        end
    end

    # Fallback canvas if empty
    if isnothing(p)
        fallback_title = q_is_zero ? 
            "No global contours found for p = ±$(target_p), q = $(target_q)" :
            "No global contours found for p = ±$(target_p), q = ±$(target_q)"
        p = plot(
            xlabel = L"\gamma",
            ylabel = L"E",
            title = fallback_title
        )
    end

    return p
end


function plt_joint_band_accumulated_chern_extrema_vs_gamma(
    gamma_data_pairs::Vector{<:Tuple{Real, Any}}; 
    tol::Float64=0.01,
    ylims::Union{Nothing, Tuple{Float64, Float64}}=nothing
)
    # Buffers for accepted max and min Chern states
    g_max, e_max = Float64[], Float64[]
    g_min, e_min = Float64[], Float64[]

    for (gamma, band_data) in gamma_data_pairs
        g = Float64(gamma)

        # 1. Flatten energy across all bands
        E_flat = vec(band_data.plaquette_energies)

        # 2. Compute global energy-accumulated Chern numbers across overlapping bands
        cum_flat = if hasproperty(band_data, :cum_chern_global)
            vec(band_data.cum_chern_global)
        elseif hasproperty(band_data, :berry_curvature)
            b_flat = vec(band_data.berry_curvature) ./ (2π)
            sort_idx = sortperm(E_flat)
            cum_sorted = cumsum(b_flat[sort_idx])
            cum_resorted = similar(cum_sorted)
            cum_resorted[sort_idx] = cum_sorted
            cum_resorted
        else
            vec(band_data.cum_chern_per_band)
        end

        # 3. Dynamic extrema detection for this specific gamma
        c_max_val = maximum(cum_flat)
        c_min_val = minimum(cum_flat)

        # 4. Filter indices within tolerance of max and min
        idx_max = findall(v -> abs(v - c_max_val) <= tol, cum_flat)
        idx_min = findall(v -> abs(v - c_min_val) <= tol, cum_flat)

        for i in idx_max
            push!(g_max, g)
            push!(e_max, E_flat[i])
        end

        # Prevent duplicate series overlay if max and min coincide (e.g., trivial zero bands)
        if abs(c_max_val - c_min_val) > tol
            for i in idx_min
                push!(g_min, g)
                push!(e_min, E_flat[i])
            end
        end
    end

    # Series configurations
    series_configs = [
        (g_max, e_max, L"\mathcal{C}_{\text{max}}", :crimson, :circle),
        (g_min, e_min, L"\mathcal{C}_{\text{min}}", :steelblue, :diamond)
    ]

    p = plot(
        xlabel = L"\gamma",
        ylabel = L"E",
        title = L"\text{Energies at } \mathcal{C}_{\text{max}} \text{ and } \mathcal{C}_{\text{min}} \text{ vs } \gamma",
        legend = :best
    )

    if !isnothing(ylims)
        plot!(p, ylims = ylims)
    end

    # Overlay non-empty extrema series
    for (g_pts, e_pts, lbl, clr, shp) in series_configs
        if !isempty(g_pts)
            scatter!(
                p,
                g_pts, e_pts,
                label = lbl,
                color = clr,
                markershape = shp,
                markersize = 3,
                markerstrokewidth = 0,
                alpha = 0.6
            )
        end
    end

    return p
end


function plt_joint_band_cum_chern_vs_energy(
    gamma_data_pair::Tuple{Real, Any}; 
    ylims::Union{Nothing, Tuple{Float64, Float64}}=nothing,
    peak_tol::Float64=0.01,
    e_cluster_tol::Float64=0.1,
    flat_span_tol::Float64=0.02
)
    gamma, band_data = gamma_data_pair
    g = Float64(gamma)

    p = plot(
        xlabel = L"E",
        ylabel = L"\mathcal{C}(E)",
        title = "Accumulated Chern Number vs Energy " * L"(\gamma = %$(g))",
        legend = :best,
        ylims = (-0.2, 1.01)
    )

    # 1. Flatten energy and accumulated Chern
    E_flat = vec(band_data.plaquette_energies)
    cum_flat = if hasproperty(band_data, :cum_chern_global)
        vec(band_data.cum_chern_global)
    elseif hasproperty(band_data, :berry_curvature)
        b_flat = vec(band_data.berry_curvature) ./ (2π)
        sort_idx = sortperm(E_flat)
        cum_sorted = cumsum(b_flat[sort_idx])
        cum_resorted = similar(cum_sorted)
        cum_resorted[sort_idx] = cum_sorted
        cum_resorted
    else
        vec(band_data.cum_chern_per_band)
    end

    # 2. Sort by energy
    sort_idx = sortperm(E_flat)
    E_sorted = E_flat[sort_idx]
    cum_sorted = cum_flat[sort_idx]

    # 3. Plot curve
    plot!(p, E_sorted, cum_sorted, label = L"\mathcal{C}(E, \gamma=%$(g))", alpha = 0.6)

    # --- 4. Peak Detection via Peaks.jl for candidate local peaks ---
    pks = findmaxima(cum_sorted)
    valid_indices = Int[]
    if !isempty(pks.indices)
        prom_res = peakproms(pks.indices, cum_sorted)
        proms = prom_res isa Tuple ? prom_res[2] : (hasproperty(prom_res, :proms) ? prom_res.proms : prom_res)

        range_c = maximum(cum_sorted) - minimum(cum_sorted)
        valid_indices = pks.indices[proms .>= peak_tol * range_c]

        sort!(valid_indices, by = idx -> cum_sorted[idx], rev = true)
    end

    # --- Global Maximum (Captures flat top left & right edges) ---
    c_max_val = maximum(cum_sorted)
    idx_max_all = findall(v -> abs(v - c_max_val) <= 1e-3, cum_sorted)

    if !isempty(idx_max_all)
        # Group contiguous index runs (plateaus)
        max_clusters = Vector{Int}[]
        curr = [idx_max_all[1]]
        for idx in idx_max_all[2:end]
            if idx == curr[end] + 1
                push!(curr, idx)
            else
                push!(max_clusters, curr)
                curr = [idx]
            end
        end
        push!(max_clusters, curr)

        # Extract left and right boundary indices for flat tops
        max_boundary_indices = Int[]
        for cl in max_clusters
            push!(max_boundary_indices, cl[1]) # Left edge
            if length(cl) > 1 && abs(E_sorted[cl[end]] - E_sorted[cl[1]]) > 1e-4
                push!(max_boundary_indices, cl[end]) # Right edge
            end
        end
        unique!(max_boundary_indices)

        for (i, idx) in enumerate(max_boundary_indices)
            E_max = E_sorted[idx]
            lbl = if length(max_boundary_indices) == 1
                L"E(\mathcal{C}_{\text{max}}) = %$(round(E_max, digits=3))"
            elseif length(max_boundary_indices) == 2 && length(max_clusters) == 1
                i == 1 ? L"E(\mathcal{C}_{\text{max, left}}) = %$(round(E_max, digits=3))" :
                         L"E(\mathcal{C}_{\text{max, right}}) = %$(round(E_max, digits=3))"
            else
                L"E(\mathcal{C}_{\text{max, %$i}}) = %$(round(E_max, digits=3))"
            end

            vline!(
                p, [E_max], linestyle = :dash, color = :red, alpha = 0.7,
                label = lbl
            )
        end
    end

    # --- Second Local Maximum ---
    # Filter out any candidates that belong to the global maximum plateau
    second_max_candidates = filter(idx -> abs(cum_sorted[idx] - c_max_val) > 1e-4, valid_indices)
    if !isempty(second_max_candidates)
        idx_2nd = second_max_candidates[1]
        E_2nd = E_sorted[idx_2nd]
        vline!(
            p, [E_2nd], linestyle = :dash, color = :orange, alpha = 0.7,
            label = L"E(\mathcal{C}_{\text{2nd max}}) = %$(round(E_2nd, digits=3))"
        )
    end

    # # --- Global Minimum (Captures flat bottoms & multiple distinct minima) ---
    # c_min_val = minimum(cum_sorted)
    # idx_min_all = findall(v -> abs(v - c_min_val) <= 1e-5, cum_sorted)

    # if !isempty(idx_min_all)
    #     min_clusters = Vector{Int}[]
    #     curr = [idx_min_all[1]]
    #     for idx in idx_min_all[2:end]
    #         if idx == curr[end] + 1
    #             push!(curr, idx)
    #         else
    #             push!(min_clusters, curr)
    #             curr = [idx]
    #         end
    #     end
    #     push!(min_clusters, curr)

    #     min_boundary_indices = Int[]
    #     for cl in min_clusters
    #         push!(min_boundary_indices, cl[1]) # Left edge / 1st minimum
    #         if length(cl) > 1 && abs(E_sorted[cl[end]] - E_sorted[cl[1]]) > 1e-4
    #             push!(min_boundary_indices, cl[end]) # Right edge of flat bottom
    #         end
    #     end
    #     unique!(min_boundary_indices)

    #     for (i, idx) in enumerate(min_boundary_indices)
    #         E_min = E_sorted[idx]
    #         lbl = if length(min_boundary_indices) == 1
    #             L"E(\mathcal{C}_{\text{min}}) = %$(round(E_min, digits=3))"
    #         elseif length(min_boundary_indices) == 2 && length(min_clusters) == 1
    #             i == 1 ? L"E(\mathcal{C}_{\text{min, left}}) = %$(round(E_min, digits=3))" :
    #                      L"E(\mathcal{C}_{\text{min, right}}) = %$(round(E_min, digits=3))"
    #         else
    #             L"E(\mathcal{C}_{\text{min, %$i}}) = %$(round(E_min, digits=3))"
    #         end

    #         vline!(
    #             p, [E_min], linestyle = :dash, color = :blue, alpha = 0.7,
    #             label = lbl
    #         )
    #     end
    # end


    # --- Global Minimum (Grouped by Energy Distance) ---
    c_min_val = minimum(cum_sorted)
    idx_min_all = findall(v -> abs(v - c_min_val) <= 1e-5, cum_sorted)

    if !isempty(idx_min_all)
        # Cluster candidate points by energy gap rather than index gap
        min_clusters = Vector{Int}[]
        curr = [idx_min_all[1]]
        
        for idx in idx_min_all[2:end]
            if E_sorted[idx] - E_sorted[curr[end]] <= e_cluster_tol
                push!(curr, idx)
            else
                push!(min_clusters, curr)
                curr = [idx]
            end
        end
        push!(min_clusters, curr)

        # Extract representative boundaries or center points
        min_boundary_indices = Int[]
        for cl in min_clusters
            e_span = E_sorted[cl[end]] - E_sorted[cl[1]]
            if e_span > flat_span_tol
                # Genuine flat bottom plateau: keep both left and right edges
                push!(min_boundary_indices, cl[1])
                push!(min_boundary_indices, cl[end])
            else
                # Single localized minimum region: pick the middle index
                push!(min_boundary_indices, cl[div(length(cl) + 1, 2)])
            end
        end
        unique!(min_boundary_indices)

        for (i, idx) in enumerate(min_boundary_indices)
            E_min = E_sorted[idx]
            lbl = if length(min_boundary_indices) == 1
                L"E(\mathcal{C}_{\text{min}}) = %$(round(E_min, digits=3))"
            elseif length(min_clusters) == 1 && length(min_boundary_indices) == 2
                i == 1 ? L"E(\mathcal{C}_{\text{min, left}}) = %$(round(E_min, digits=3))" :
                         L"E(\mathcal{C}_{\text{min, right}}) = %$(round(E_min, digits=3))"
            else
                L"E(\mathcal{C}_{\text{min, %$i}}) = %$(round(E_min, digits=3))"
            end

            vline!(
                p, [E_min], linestyle = :dash, color = :blue, alpha = 0.7,
                label = lbl
            )
        end
    end

    # --- 6. Mark C=0 crossing points and boundary zero values ---
    zero_tol = 1e-3
    is_zero = abs.(cum_sorted) .<= zero_tol
    sign_crossings = findall(i -> cum_sorted[i] * cum_sorted[i+1] < 0, 1:(length(cum_sorted)-1))
    zero_candidates = sort(unique(vcat(findall(is_zero), sign_crossings)))

    zero_indices = Int[]
    if !isempty(zero_candidates)
        block_start = zero_candidates[1]
        prev_idx = zero_candidates[1]
        
        for idx in zero_candidates[2:end]
            if idx > prev_idx + 1
                push!(zero_indices, block_start)
                block_start = idx
            end
            prev_idx = idx
        end
        
        if prev_idx == length(cum_sorted) && abs(cum_sorted[end]) <= zero_tol
            push!(zero_indices, length(cum_sorted))
        else
            push!(zero_indices, block_start)
        end
        unique!(zero_indices)
    end

    for (i, idx) in enumerate(zero_indices)
        E_zero = E_sorted[idx]
        vline!(
            p, [E_zero],
            linestyle = :dot, color = :black, alpha = 0.5,
            label = L"E(\mathcal{C}=0) = %$(round(E_zero, digits=3))"
        )
    end

    if !isnothing(ylims)
        plot!(p, ylims = ylims)
    end

    return p
end

# Vector overload: returns an array of plots for multiple pairs
function plt_joint_band_cum_chern_vs_energy(
    gamma_data_pairs::Vector{<:Tuple{Real, Any}}; 
    ylims::Union{Nothing, Tuple{Float64, Float64}}=nothing,
    peak_tol::Float64=0.01
)
    return [plt_joint_band_cum_chern_vs_energy(pair; ylims=ylims, peak_tol=peak_tol) for pair in gamma_data_pairs]
end


function plt_joint_band_smart_extrema_accumulated_chern_energy_vs_gamma(
    gamma_data_pairs::Vector{<:Tuple{Real, Any}}; 
    peak_tol::Float64=0.01,
    e_cluster_tol::Float64=0.1,
    flat_span_tol::Float64=0.02,
    zero_tol::Float64=1e-4,
    ylims::Union{Nothing, Tuple{Float64, Float64}}=nothing
)
    # Global buffers to accumulate feature points across all gamma values
    g_max, e_max = Float64[], Float64[]
    g_2nd, e_2nd = Float64[], Float64[]
    g_min, e_min = Float64[], Float64[]
    g_zero, e_zero = Float64[], Float64[]

    for (gamma, band_data) in gamma_data_pairs
        g = Float64(gamma)

        # 1. Flatten energy and accumulated Chern
        E_flat = vec(band_data.plaquette_energies)
        cum_flat = if hasproperty(band_data, :cum_chern_global)
            vec(band_data.cum_chern_global)
        elseif hasproperty(band_data, :berry_curvature)
            b_flat = vec(band_data.berry_curvature) ./ (2π)
            sort_idx = sortperm(E_flat)
            cum_sorted = cumsum(b_flat[sort_idx])
            cum_resorted = similar(cum_sorted)
            cum_resorted[sort_idx] = cum_sorted
            cum_resorted
        else
            vec(band_data.cum_chern_per_band)
        end

        # 2. Sort by energy
        sort_idx = sortperm(E_flat)
        E_sorted = E_flat[sort_idx]
        cum_sorted = cum_flat[sort_idx]

        # --- Peak Detection via Peaks.jl ---
        pks = findmaxima(cum_sorted)
        valid_indices = Int[]
        if !isempty(pks.indices)
            prom_res = peakproms(pks.indices, cum_sorted)
            proms = prom_res isa Tuple ? prom_res[2] : (hasproperty(prom_res, :proms) ? prom_res.proms : prom_res)

            range_c = maximum(cum_sorted) - minimum(cum_sorted)
            valid_indices = pks.indices[proms .>= peak_tol * range_c]
            sort!(valid_indices, by = idx -> cum_sorted[idx], rev = true)
        end

        # --- Global Maximum (Captures flat top left & right edges) ---
        c_max_val = maximum(cum_sorted)
        idx_max_all = findall(v -> abs(v - c_max_val) <= 1e-3, cum_sorted)

        if !isempty(idx_max_all)
            max_clusters = Vector{Int}[]
            curr = [idx_max_all[1]]
            for idx in idx_max_all[2:end]
                if idx == curr[end] + 1
                    push!(curr, idx)
                else
                    push!(max_clusters, curr)
                    curr = [idx]
                end
            end
            push!(max_clusters, curr)

            max_boundary_indices = Int[]
            for cl in max_clusters
                push!(max_boundary_indices, cl[1]) # Left edge
                if length(cl) > 1 && abs(E_sorted[cl[end]] - E_sorted[cl[1]]) > 1e-4
                    push!(max_boundary_indices, cl[end]) # Right edge
                end
            end
            unique!(max_boundary_indices)

            for idx in max_boundary_indices
                push!(g_max, g)
                push!(e_max, E_sorted[idx])
            end
        end

        # --- Second Local Maximum ---
        second_max_candidates = filter(idx -> abs(cum_sorted[idx] - c_max_val) > 1e-4, valid_indices)
        if !isempty(second_max_candidates)
            idx_2nd = second_max_candidates[1]
            push!(g_2nd, g)
            push!(e_2nd, E_sorted[idx_2nd])
        end

        # --- Global Minimum (Grouped by Energy Distance) ---
        c_min_val = minimum(cum_sorted)
        idx_min_all = findall(v -> abs(v - c_min_val) <= 1e-5, cum_sorted)

        if !isempty(idx_min_all)
            min_clusters = Vector{Int}[]
            curr = [idx_min_all[1]]
            for idx in idx_min_all[2:end]
                if E_sorted[idx] - E_sorted[curr[end]] <= e_cluster_tol
                    push!(curr, idx)
                else
                    push!(min_clusters, curr)
                    curr = [idx]
                end
            end
            push!(min_clusters, curr)

            min_boundary_indices = Int[]
            for cl in min_clusters
                e_span = E_sorted[cl[end]] - E_sorted[cl[1]]
                if e_span > flat_span_tol
                    # Genuine flat bottom: include left and right boundary edges
                    push!(min_boundary_indices, cl[1])
                    push!(min_boundary_indices, cl[end])
                else
                    # Localized minimum region: pick central point
                    push!(min_boundary_indices, cl[div(length(cl) + 1, 2)])
                end
            end
            unique!(min_boundary_indices)

            for idx in min_boundary_indices
                push!(g_min, g)
                push!(e_min, E_sorted[idx])
            end
        end

        # --- Mark C=0 crossing points and boundary zero values ---
        is_zero = abs.(cum_sorted) .<= zero_tol
        sign_crossings = findall(i -> cum_sorted[i] * cum_sorted[i+1] < 0, 1:(length(cum_sorted)-1))
        zero_candidates = sort(unique(vcat(findall(is_zero), sign_crossings)))

        zero_indices = Int[]
        if !isempty(zero_candidates)
            block_start = zero_candidates[1]
            prev_idx = zero_candidates[1]
            
            for idx in zero_candidates[2:end]
                if idx > prev_idx + 1
                    push!(zero_indices, block_start)
                    block_start = idx
                end
                prev_idx = idx
            end
            
            if prev_idx == length(cum_sorted) && abs(cum_sorted[end]) <= zero_tol
                push!(zero_indices, length(cum_sorted))
            else
                push!(zero_indices, block_start)
            end
            unique!(zero_indices)
        end

        for idx in zero_indices
            push!(g_zero, g)
            push!(e_zero, E_sorted[idx])
        end
    end

    # Series configurations (Buffers, Label, Color, Shape)
    series_configs = [
        (g_max, e_max, L"E(\mathcal{C}_{\text{max}})", :crimson, :star5),
        (g_2nd, e_2nd, L"E(\mathcal{C}_{\text{2nd max}})", :darkorange, :diamond),
        (g_min, e_min, L"E(\mathcal{C}_{\text{min}})", :royalblue, :rect),
        (g_zero, e_zero, L"E(\mathcal{C}=0)", :black, :circle)
    ]

    p = scatter(
        xlabel = L"\gamma",
        ylabel = L"E",
        title = "Characteristic Chern Energies vs " * L"\gamma",
        legend = :best
    )

    for (g_pts, e_pts, lbl, clr, shp) in series_configs
        if !isempty(g_pts)
            scatter!(
                p, g_pts, e_pts,
                label = lbl,
                color = clr,
                markershape = shp,
                markersize = 4,
                markerstrokewidth = 0,
                alpha = 0.7
            )
        end
    end

    if !isnothing(ylims)
        plot!(p, ylims = ylims)
    end

    return p
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
ms = [-1.9]
gammas = collect(-3.0:0.1:3.0)
distortion_type = :symmetric
Nkx = 101
Nky = 101

data_folder = joinpath("data", "bulk_band_berry_data")
plot_folder = joinpath("plots", "bulk_band_berry_data")
isdir(data_folder) || mkpath(data_folder)
isdir(plot_folder) || mkpath(plot_folder)

generate_new_data = true
generate_new_plots = false
generate_gamma_scatter = false
generate_cum_chern_extrema = false

if generate_new_data
    println("Generating new data...")
    @showprogress for (A, B, m, gamma) in Iterators.product(As, Bs, ms, gammas)
        band_data = compute_bulk_band_berry_data(
            Nkx=101,
            Nky=101,
            A=A,
            B=B,
            m=m,
            gamma=gamma,
            perturbation_type=distortion_type
        )

        edge_state_data = edge_state_prediction(
            Nkx=101,
            Nky=101,
            A=A,
            B=B,
            m=m,
            gamma=gamma,
            perturbation_type=distortion_type
        )

        data = (
            band = band_data,
            edge = edge_state_data
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

        # # 4. Bare bandstructure heatmaps for both bands
        # bandstruct_heatmaps = plt_bandstructure_heatmap(
        #     data.band.kx_vals,
        #     data.band.ky_vals,
        #     data.band.energies;
        #     title = L"\gamma = %$(gamma)",
        #     colour = :curl,
        #     share_colour_scale = false
        # )
        
        # savefig(
        #     bandstruct_heatmaps, 
        #     joinpath(plot_folder, "bulk_bandstruct_A$(A)_B$(B)_m$(m)_gamma$(gamma)_perturbation_$(distortion_type).png")
        # )

        # # 5. Plot the Berry curvature heatmaps. F_xy flux per plaquette resolved for (kx, ky) for both bands.
        # F_xy_heatmaps = plt_k_resolved_F_xy_heatmaps(
        #     data.band.kx_vals,
        #     data.band.ky_vals,
        #     data.band.berry_curvature;
        #     title = L"\gamma = %$(gamma), F_{xy} / 2\pi",
        #     # colour = :RdBu,
        #     # share_colour_scale = false
        # )

        # savefig(
        #     F_xy_heatmaps, 
        #     joinpath(plot_folder, "bulk_band_dEdC_A$(A)_B$(B)_m$(m)_gamma$(gamma)_perturbation_$(distortion_type).png")
        # )

        # # 6. Plot the accumulated chern number as a heatmap on (kx, ky) grid
        # cum_chern_heatmaps = plt_accumulated_chern_heatmaps(
        #     data.band.kx_vals,
        #     data.band.ky_vals,
        #     data.band.cum_chern_per_band;
        #     title = L"\gamma = %$(gamma), \mathcal{C}(k_x, k_y)",
        #     colour = :viridis
        # )

        # savefig(
        #     cum_chern_heatmaps, 
        #     joinpath(plot_folder, "bulk_band_accumulated_chern_A$(A)_B$(B)_m$(m)_gamma$(gamma)_perturbation_$(distortion_type).png")
        # )

        # # 7. Plot the bandgap heatmap and edge state predictions
        # edge_state_heatmap = plt_edge_state_prediction(
        #     data.edge;
        #     title = L"\gamma = %$(gamma), \Delta E(k_x, k_y) \text{ and Edge State Predictions}",
        #     colour = :thermal
        # )

        # savefig(
        #     edge_state_heatmap, 
        #     joinpath(plot_folder, "bulk_band_edge_state_prediction_A$(A)_B$(B)_m$(m)_gamma$(gamma)_perturbation_$(distortion_type).png")
        # )

        # 8. Plot the heatmap of accumulated berry curvature given overlapping bands at fermi energy E_cut
        energies = data.band.plaquette_energies
        E_cuts = collect(range(minimum(energies), stop=maximum(energies), length=50))

        # Pre-compute global clims across ALL E_cuts for this specific gamma
        n_bands, nkx, nky = size(energies)
        max_abs_gamma = 0.0

        for E_cut in E_cuts
            cum_2d = zeros(Float64, nkx, nky)
            for n in 1:n_bands, i in 1:nkx, j in 1:nky
                if energies[n, i, j] <= E_cut
                    cum_2d[i, j] += data.band.berry_curvature[n, i, j] / (2π)
                end
            end
            max_abs_gamma = max(max_abs_gamma, maximum(abs, cum_2d))
        end

        # Set symmetric bounds around zero
        global_clims = (-max_abs_gamma, max_abs_gamma)

        for (idx, E_cut) in enumerate(E_cuts)

            cum_chern_E_cut_heatmap = plt_joint_band_accumulated_chern_heatmap(
                data.band.kx_vals,
                data.band.ky_vals,
                data.band.berry_curvature,
                data.band.plaquette_energies,
                E_cut,
                gamma;
                global_clims=global_clims,
                colour = :RdBu
            )

            savefig(
                cum_chern_E_cut_heatmap, 
                joinpath(plot_folder, "$(idx)_bulk_band_accumulated_chern_Ecut_$(round(E_cut, digits=2))_A$(A)_B$(B)_m$(m)_gamma$(gamma)_perturbation_$(distortion_type).png")
            )
        end
    end
    
end



if generate_gamma_scatter
    println("Generating gamma scatter plot...")
    gamma_data_list = Tuple{Float64, Any}[]
    data_to_unpack = joinpath(data_folder, "bandstrucutre_data")

    for (A, B, m, gamma) in Iterators.product(As, Bs, ms, gammas)
        filename = "bulk_band_berry_data_A$(A)_B$(B)_m$(m)_gamma$(gamma)_perturbation_$(distortion_type).jld2"
        filepath = joinpath(data_to_unpack, filename)

        if isfile(filepath)
            @load filepath data

            # Robust unpacking for NamedTuple or Dict[cite: 2]
            band_data = if hasproperty(data, :band)
                data.band
            elseif data isa AbstractDict && haskey(data, :band)
                data[:band]
            else
                data
            end

            push!(gamma_data_list, (gamma, band_data))
        end
    end

    if !isempty(gamma_data_list)
        band_separated_plt = plt_accumulated_chern_energy_vs_gamma(
            gamma_data_list; 
            target_p = 0.5, 
            target_q = 0.05,
            tol = 0.005
        )

        m_val = ms[1]
        A_val = As[1]
        B_val = Bs[1]

        # gamma_range_str = "$(minimum(gammas))_to_$(maximum(gammas))"
        savefig(band_separated_plt, joinpath(plot_folder, "energy_at_accumulated_chern_0.5_vs_gamma$(minimum(gammas))-$(maximum(gammas))_A$(A_val)_B$(B_val)_m$(m_val).png"))


        band_joint_plot = plt_joint_band_accumulated_chern_energy_vs_gamma(
            gamma_data_list; 
            target_p = 0.5, 
            target_q = 0.0,
            tol = 0.005
        )

        savefig(band_joint_plot, joinpath(plot_folder, "energy_at_joint_band_accumulated_chern_0.5_vs_gamma$(minimum(gammas))-$(maximum(gammas))_A$(A_val)_B$(B_val)_m$(m_val).png"))

        band_extrema_plot = plt_joint_band_accumulated_chern_extrema_vs_gamma(
            gamma_data_list; 
            tol = 0.005
        )

        savefig(band_extrema_plot, joinpath(plot_folder, "energy_at_joint_band_accumulated_chern_extrema_vs_gamma$(minimum(gammas))-$(maximum(gammas))_A$(A_val)_B$(B_val)_m$(m_val).png"))

        joint_tracking_plot = plt_joint_band_smart_extrema_accumulated_chern_energy_vs_gamma(
            gamma_data_list; 
            peak_tol = 0.01,
            e_cluster_tol = 0.1,
            flat_span_tol = 0.02,
            zero_tol = 1e-3
        )

        savefig(joint_tracking_plot, joinpath(plot_folder, "energy_at_joint_band_accumulated_chern_extrema_and_zero_vs_gamma$(minimum(gammas))-$(maximum(gammas))_A$(A_val)_B$(B_val)_m$(m_val).png"))
    end
end


if generate_cum_chern_extrema
    println("Generating cumulative chern line plots with extrema marked...")
    subfolder = joinpath(plot_folder, "cumulative_chern_extrema")
    isdir(subfolder) || mkpath(subfolder)

    data_to_unpack = joinpath(data_folder, "bandstrucutre_data")

    for (A, B, m, gamma) in Iterators.product(As, Bs, ms, gammas)
        filename = "bulk_band_berry_data_A$(A)_B$(B)_m$(m)_gamma$(gamma)_perturbation_$(distortion_type).jld2"
        filepath = joinpath(data_to_unpack, filename)

        if isfile(filepath)
            @load filepath data

            band_data = if hasproperty(data, :band)
                data.band
            elseif data isa AbstractDict && haskey(data, :band)
                data[:band]
            else
                data
            end

            # Generate individual plot for this gamma
            cum_chern_extrema_plt = plt_joint_band_cum_chern_vs_energy((gamma, band_data); peak_tol=0.01)

            out_filename = "cumulative_chern_vs_energy_with_extrema_marked_A$(A)_B$(B)_m$(m)_gamma$(gamma).png"
            savefig(cum_chern_extrema_plt, joinpath(subfolder, out_filename))
        end
    end
end