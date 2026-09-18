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

function compute_euler_characteristic(
    M::AbstractMatrix{Bool}; 
    periodic::Bool=true
)
    """
        compute_euler_characteristic(M::AbstractMatrix{Bool}; periodic::Bool=true)

    Computes the Euler characteristic χ = V - E + F on a 2D lattice mask M 
    for the occupied Fermi sea (E(k) <= E_F).
    """

    Nx, Ny = size(M)
    V = count(M)
    Eh, Ev, F = 0, 0, 0

    for j in 1:Ny, i in 1:Nx
        if M[i, j]
            i_next = (periodic && i == Nx) ? 1 : i + 1
            j_next = (periodic && j == Ny) ? 1 : j + 1

            # Horizontal & vertical edge connectivity
            if j_next <= Ny && M[i, j_next]
                Eh += 1
            end
            if i_next <= Nx && M[i_next, j]
                Ev += 1
            end
            # 2x2 Plaquette / Face
            if i_next <= Nx && j_next <= Ny && M[i_next, j] && M[i, j_next] && M[i_next, j_next]
                F += 1
            end
        end
    end

    return V - (Eh + Ev) + F
end

function plt_euler_characteristic_heatmap_vs_gamma(
    gamma_data_pairs::Vector{<:Tuple{Real, Any}};
    n_energies::Int = 300,
    energy_lims::Union{Nothing, Tuple{Float64, Float64}} = nothing,
    periodic::Bool = true,
    add_contours::Bool = true,
)
    gammas = Float64[g for (g, _) in gamma_data_pairs]

    # 1. Determine global energy range across all plaquette data
    min_E, max_E = if isnothing(energy_lims)
        g_min, g_max = Inf, -Inf
        for (_, band_data) in gamma_data_pairs
            E_flat = vec(band_data.plaquette_energies)
            g_min = min(g_min, minimum(E_flat))
            g_max = max(g_max, maximum(E_flat))
        end
        (g_min, g_max)
    else
        energy_lims
    end

    energy_grid = range(min_E, max_E, length=n_energies)
    chi_matrix = zeros(Int, n_energies, length(gamma_data_pairs))

    # 2. Compute χ(E_F) across all bands on the k-space grid
    for (g_idx, (_, band_data)) in enumerate(gamma_data_pairs)
        E_data = band_data.plaquette_energies
        n_dims = ndims(E_data)
        n_bands = n_dims == 3 ? size(E_data, 1) : 1

        for (e_idx, E_F) in enumerate(energy_grid)
            chi_total = 0
            if n_dims == 3
                for b in 1:n_bands
                    M = @view(E_data[b, :, :]) .<= E_F
                    chi_total += compute_euler_characteristic(M; periodic=periodic)
                end
            else
                M = E_data .<= E_F
                chi_total += compute_euler_characteristic(M; periodic=periodic)
            end
            chi_matrix[e_idx, g_idx] = chi_total
        end
    end

    # 3. Render 2D heatmap plot
    p = heatmap(
        gammas,
        energy_grid,
        chi_matrix,
        xlabel = L"\gamma",
        ylabel = L"E_F",
        title = "Fermi Sea Euler Characteristic " * L" \chi(E_F, \gamma)",
        color = :viridis,
        colorbar_title = L"\chi",
        clims = (minimum(chi_matrix), maximum(chi_matrix))
    )

    if add_contours
        contour!(
            p,
            gammas,
            energy_grid,
            chi_matrix,
            levels=unique(chi_matrix),
            color=:red,
            linewidth=1.0,
            linestyle=:dash,
            contour_labels=true,
            colorbar_entry=false
        )
    end

    return p
end

############################################################################
"""
Computing the PBS bandstructure and its berry curvature for distorted QWZ model
"""

As = [1.0]
Bs = [1.0]
ms = [-5.0]
gammas = collect(-0.5:0.001:1.5)
distortion_type = :symmetric
Nkx = 501
Nky = 501

data_folder = joinpath("data", "bulk_band_berry_data")
plot_folder = joinpath("plots", "bulk_band_berry_data")
isdir(data_folder) || mkpath(data_folder)
isdir(plot_folder) || mkpath(plot_folder)

generate_new_data = true
generate_new_plots = false
generate_gamma_scatter = true
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
