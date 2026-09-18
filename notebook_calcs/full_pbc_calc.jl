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
    elseif perturbation_type == :antisymmetric
        # gamma * sin(kx+ky) * I
        tx += ComplexF64.(0.5im * gamma * identity)
        ty += ComplexF64.(0.5im * gamma * identity)
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
    winding_number::Int=0,
    phi::Real=0.0,
    orbital_displacement::Real=0.0
)::Matrix{ComplexF64}
    # Orbital sublattice embedding parameters:
    #   phi::Real - Angular position of orbital displacement in 2D k-space (radians)
    #   orbital_displacement::Real - Magnitude of displacement vector d(phi) = |d|*(cos(phi), sin(phi))
    # Applies unitary transformation: exp(-i*theta*sigma_z) * H * exp(i*theta*sigma_z)
    # where theta = k_x*d_x(phi) + k_y*d_y(phi)

    winding_number >= 0 || error("winding_number must be non-negative.")

    # Base QWZ Hamiltonian in k-space
    d_x = A * sin(kx)
    d_y = A * sin(ky)
    d_z = m + 2B - B * cos(kx) - B * cos(ky)
    scalar_term = 0.0

    # Apply scalar perturbations independently of the winding distortion.
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
    # Legacy simplified names (for backward compatibility)
    elseif perturbation_type == :symmetric
        scalar_term += gamma * (cos(kx) + cos(ky))
    elseif perturbation_type == :antisymmetric
        scalar_term += gamma * sin(kx + ky)
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

    # Build Hamiltonian before embedding
    H_k = scalar_term * identity + d_x * sigma_x + d_y * sigma_y + d_z * sigma_z

    # Apply orbital embedding based on selected implementation
    if abs(orbital_displacement) > 1e-14
        dx_phi = orbital_displacement * cos(phi)
        dy_phi = orbital_displacement * sin(phi)
        
        # ============================================
        # IMPLEMENTATION: σ_z basis rotation
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
            
    end

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
        elseif perturbation_type == :antisymmetric
            scalar_term = gamma * sin(kx + ky)
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

# #only the euler char
# function compute_euler_characteristic_vs_energy(
#     bulk_data; 
#     E_mesh = range(-3.0, 3.0, length=101)
# )
#     # Get all bands' energies and sum Euler characteristic across all bands (matches original plotting function)
#     E_data = bulk_data.plaquette_energies
#     n_dims = ndims(E_data)
#     n_bands = n_dims == 3 ? size(E_data, 1) : 1
    
#     euler_vs_E = Float64[]
    
#     for EF in E_mesh
#         chi_total = 0
#         if n_dims == 3
#             for b in 1:n_bands
#                 M = @view(E_data[b, :, :]) .<= EF
#                 chi_total += compute_euler_characteristic(M; periodic=true)
#             end
#         else
#             M = E_data .<= EF
#             chi_total = compute_euler_characteristic(M; periodic=true)
#         end
#         push!(euler_vs_E, chi_total)
#     end
    
#     return (E_mesh = collect(E_mesh), euler_char = euler_vs_E)
# end

## including the N_pockets and N_hoels decomp
function compute_euler_characteristic_vs_energy(
    bulk_data; 
    E_mesh = range(-3.0, 3.0, length=101),
    periodic::Bool = true
)
    # Extract energy data across band dimensions
    E_data = hasproperty(bulk_data, :plaquette_energies) ? bulk_data.plaquette_energies : bulk_data[:plaquette_energies]
    n_dims = ndims(E_data)
    n_bands = n_dims == 3 ? size(E_data, 1) : 1
    
    euler_vs_E = Int[]
    pockets_vs_E = Int[]
    holes_vs_E = Int[]
    
    for EF in E_mesh
        chi_total = 0
        pockets_total = 0
        holes_total = 0

        for b in 1:n_bands
            M = n_dims == 3 ? @view(E_data[b, :, :]) .<= EF : E_data .<= EF
            
            # Compute independent topological counts
            n_pockets, n_holes = compute_pockets_and_holes(M; periodic=periodic)
            chi = compute_euler_characteristic(M; periodic=periodic)

            pockets_total += n_pockets
            holes_total += n_holes
            chi_total += chi
        end
        
        push!(euler_vs_E, chi_total)
        push!(pockets_vs_E, pockets_total)
        push!(holes_vs_E, holes_total)
    end
    
    return (
        E_mesh = collect(E_mesh), 
        euler_char = euler_vs_E,
        n_pockets = pockets_vs_E,
        n_holes = holes_vs_E
    )
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
            n_pockets += 1
            visited[i, j] = true
            queue = [(i, j)]

            # BFS to traverse the connected Fermi pocket
            while !isempty(queue)
                curr_i, curr_j = popfirst!(queue)

                for (di, dj) in dirs
                    ni = curr_i + di
                    nj = curr_j + dj

                    # Apply periodic boundary conditions or grid bounds
                    if periodic
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
        end
    end

    n_holes = n_pockets - chi
    return n_pockets, n_holes
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

function plt_accumulated_chern_2d_inflections(
    kx_vals::Vector{Float64}, 
    ky_vals::Vector{Float64}, 
    cum_chern_per_band::Array{Float64, 3}; 
    title::LaTeXString=L"", 
    xlabel::LaTeXString=L"k_x",
    ylabel::LaTeXString=L"k_y", 
    colour=:viridis, 
    shift_to_centers::Bool=false,
    add_contour::Bool=true,
    smooth_window::Int=15,       # 2D moving-average window size (odd integer)
    min_slope_tol::Float64=1e-3  # Gradient magnitude tolerance to mask out flat plateau noise
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
        # Transpose matrix to match Plots.jl (y=row, x=col) alignment
        cum_transposed = cum_chern_per_band[band_index, :, :]'
        
        min_val = minimum(cum_transposed)
        max_val = maximum(cum_transposed)
        total_chern = abs(max_val) > abs(min_val) ? round(max_val, digits=4) : round(min_val, digits=4)
        
        sub_title = "Chern accumulated: " * L"\mathcal{C}(k_x, k_y)" * " (Band $band_index, C = $total_chern)"
        
        p = heatmap(
            x_coords, 
            y_coords, 
            cum_transposed, 
            xlabel=xlabel, 
            ylabel=ylabel, 
            title=sub_title, 
            color=colour, 
            clims=(min_val, max_val), 
            xlims=(-π, π), 
            ylims=(-π, π), 
            xticks=pi_ticks, 
            yticks=pi_ticks, 
            aspect_ratio=:equal, 
            colorbar_title=L"\mathcal{C}(k_x, k_y)"
        )

        if add_contour
            # --- 1. Apply 2D Moving Average Box Filter ---
            Z_smooth = copy(cum_transposed)
            w = max(1, smooth_window)
            r = div(w, 2)
            ny, nx = size(cum_transposed)

            if r > 0
                for iy in 1:ny, ix in 1:nx
                    # Periodic index wrapping across Brillouin Zone boundaries
                    y_idxs = mod1.((iy - r):(iy + r), ny)
                    x_idxs = mod1.((ix - r):(ix + r), nx)
                    Z_smooth[iy, ix] = mean(cum_transposed[y_idxs, x_idxs])
                end
            end

            # --- 2. Periodic Finite Differences on Smoothed Array ---
            Z_left  = circshift(Z_smooth, (0, 1))
            Z_right = circshift(Z_smooth, (0, -1))
            Z_up    = circshift(Z_smooth, (1, 0))
            Z_down  = circshift(Z_smooth, (-1, 0))

            # First Derivatives (Central Difference for Gradient Magnitude)
            C_x = (Z_right .- Z_left) ./ (2 * dkx)
            C_y = (Z_down .- Z_up) ./ (2 * dky)
            grad_mag = sqrt.(C_x.^2 .+ C_y.^2)

            # Second Derivatives
            C_xx = (Z_right .- 2 .* Z_smooth .+ Z_left) ./ (dkx^2)
            C_yy = (Z_down  .- 2 .* Z_smooth .+ Z_up)   ./ (dky^2)

            # 2D Laplacian
            laplacian = C_xx .+ C_yy

            # --- 3. Mask Laplacian with Gradient Tolerance ---
            # Set laplacian to NaN where gradient is below tolerance to suppress noise in flat regions
            laplacian_masked = copy(laplacian)
            laplacian_masked[grad_mag .< min_slope_tol] .= NaN

            # --- 4. Overlay Inflection Contour (∇²C = 0) ---
            contour!(
                p,
                x_coords,
                y_coords,
                laplacian_masked,
                levels=[0.0],
                color=:white,
                linewidth=1.5,
                linestyle=:dash,
                contour_labels=false,
                colorbar_entry=false
            )
        end

        return p
    end

    return plot(plots[1], plots[2], layout=(1, 2), size=(1000, 450))
end

function plt_accumulated_chern_1d_inflections(
    kx_vals::Vector{Float64}, 
    ky_vals::Vector{Float64}, 
    cum_chern_per_band::Array{Float64, 3}; 
    plaquette_energies::Union{Array{Float64, 3}, Nothing}=nothing,
    title::LaTeXString=L"", 
    xlabel::LaTeXString=L"k_x",
    ylabel::LaTeXString=L"k_y", 
    colour=:viridis, 
    shift_to_centers::Bool=false,
    add_contour::Bool=true,
    smooth_window::Int=25,        # Window size to smooth 1D accumulation curve
    min_slope_tol::Float64=1e-4   # Gradient tolerance to ignore flat plateaus
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
        sub_title = "Chern accumulated: " * L"C(E, k_x, k_y)" * " (Band $band_index)"

        p = heatmap(
            x_coords, 
            y_coords, 
            cum_transposed, 
            xlabel=xlabel, 
            ylabel=ylabel, 
            title=sub_title, 
            color=colour, 
            clims=(min_val, max_val), 
            xlims=(-π, π), 
            ylims=(-π, π), 
            xticks=pi_ticks, 
            yticks=pi_ticks, 
            aspect_ratio=:equal, 
            colorbar_title=L"C(E, k_x, k_y)"
        )

        if add_contour
            c_flat = vec(cum_transposed)

            # 1. Order Chern values by energy (if provided) or monotonic sorting
            if !isnothing(plaquette_energies)
                e_flat = vec(plaquette_energies[band_index, :, :]')
                c_sorted = c_flat[sortperm(e_flat)]
            else
                c_sorted = sort(c_flat)
            end

            # 2. Moving average smoothing on 1D accumulated curve
            w = max(3, smooth_window)
            half_w = div(w, 2)
            N = length(c_sorted)
            c_smoothed = copy(c_sorted)
            for i in (half_w + 1):(N - half_w)
                c_smoothed[i] = mean(c_sorted[(i - half_w):(i + half_w)])
            end

            # 3. Compute 1D Finite Differences
            d1 = diff(c_smoothed)
            d2 = diff(d1)

            # 4. Find 1D Inflection Points (d²C/dE² = 0 with sign change)
            inflection_levels = Float64[]
            eps_shift = 1e-4

            for i in 1:(length(d2) - 1)
                if d2[i] * d2[i+1] <= 0 && abs(d1[i]) > min_slope_tol
                    lvl = c_smoothed[i + 1]
                    # Clamp levels slightly inside [min_val, max_val] domain boundary[cite: 1]
                    clamped_lvl = clamp(lvl, min_val + eps_shift, max_val - eps_shift)
                    if min_val < clamped_lvl < max_val
                        push!(inflection_levels, clamped_lvl)
                    end
                end
            end

            # Deduplicate near-identical inflection values
            inflection_levels = unique(round.(inflection_levels, digits=3))

            # 5. Overlay 1D inflection contours onto 2D heatmap
            if !isempty(inflection_levels)
                contour!(
                    p,
                    x_coords,
                    y_coords,
                    cum_transposed,
                    levels=inflection_levels,
                    color=:white,
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

function plt_euler_characteristic_heatmap_vs_gamma_precalc(
    gamma_data_pairs::Vector{<:Tuple{Real, Any}};
    add_contours::Bool = true,
)
    """
    Creates a heatmap of Euler characteristic vs gamma using pre-calculated data.
    Each element in gamma_data_pairs should be (gamma, full_data) where full_data 
    contains the :euler_char field with pre-computed E_mesh and euler_char.
    """
    gammas = Float64[g for (g, _) in gamma_data_pairs]
    
    # Extract pre-calculated data
    E_mesh = nothing
    chi_matrix = []
    
    for (g_idx, (gamma, full_data)) in enumerate(gamma_data_pairs)
        # Extract euler_char data from full_data structure
        euler_char_data = if hasproperty(full_data, :euler_char)
            full_data.euler_char
        elseif full_data isa AbstractDict && haskey(full_data, :euler_char)
            full_data[:euler_char]
        else
            error("Data does not contain :euler_char field for gamma=$gamma")
        end
        
        # Initialize chi_matrix on first iteration
        if E_mesh === nothing
            E_mesh = euler_char_data.E_mesh
            chi_matrix = zeros(Int, length(E_mesh), length(gamma_data_pairs))
        end
        
        # Copy pre-calculated euler characteristic values
        chi_matrix[:, g_idx] = euler_char_data.euler_char
    end
    
    # Render 2D heatmap plot
    p = heatmap(
        gammas,
        E_mesh,
        chi_matrix,
        xlabel = L"\gamma",
        ylabel = L"E_F",
        title = "Fermi Sea Euler Characteristic (Pre-calculated) " * L" \chi(E_F, \gamma)",
        color = :viridis,
        colorbar_title = L"\chi",
        clims = (minimum(chi_matrix), maximum(chi_matrix))
    )

    if add_contours
        contour!(
            p,
            gammas,
            E_mesh,
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

function plt_cum_chern_vs_energy(
    accumulation_energies::Vector{Vector{Float64}},
    cumulative_chern::Vector{Vector{Float64}};
    title::LaTeXString=L"",
    ylabel::LaTeXString=L"E",
    xlabel::LaTeXString=L"\mathcal{C}(E)",
    colors::Vector{Symbol}=[:royalblue, :crimson],
    joint_color::Symbol=:darkorange,
    linewidth::Float64=2.0,
    band_labels::Vector{String}=["Band 1 (lower)", "Band 2 (upper)"]
)
    """
    Plots energy-accumulated Chern number C(E) vs energy for both bands individually 
    and their combined joint sum.
    """
    p = plot(
        xlabel=xlabel,
        ylabel=ylabel,
        title=title,
        legend=:best,
        grid=true
    )

    nbands = min(length(accumulation_energies), length(cumulative_chern))
    
    # 1. Plot individual band contributions
    for band_idx in 1:nbands
        lbl = band_idx <= length(band_labels) ? band_labels[band_idx] : "Band $(band_idx)"
        clr = colors[mod1(band_idx, length(colors))]

        plot!(
            p,
            cumulative_chern[band_idx],
            accumulation_energies[band_idx],
            label=lbl,
            color=clr,
            linewidth=linewidth
        )
    end

    # 2. Compute and plot global joint contribution
    all_energies = reduce(vcat, accumulation_energies[1:nbands])
    
    # Calculate differential Chern steps for each band
    diff_chern = map(1:nbands) do b
        c = cumulative_chern[b]
        isempty(c) ? Float64[] : vcat(c[1], diff(c))
    end
    all_diff_chern = reduce(vcat, diff_chern)

    # Sort unified energy spectrum to accumulate joint Chern number
    sort_idx = sortperm(all_energies)
    joint_energies = all_energies[sort_idx]
    joint_cum_chern = cumsum(all_diff_chern[sort_idx])

    plot!(
        p,
        joint_cum_chern,
        joint_energies,
        label="Joint (total)",
        color=joint_color,
        linewidth=linewidth,
        linestyle=:dash
    )

    return p
end

# Overload for 3D Arrays [nbands × kx × ky]
function plt_cum_chern_vs_energy(
    accumulation_energies::Array{Float64, 3},
    cumulative_chern::Array{Float64, 3};
    kwargs...
)
    nbands = size(accumulation_energies, 1)
    energies_vec = Vector{Vector{Float64}}(undef, nbands)
    chern_vec = Vector{Vector{Float64}}(undef, nbands)

    for band in 1:nbands
        e_flat = vec(accumulation_energies[band, :, :])
        c_flat = vec(cumulative_chern[band, :, :])

        sort_order = sortperm(e_flat)
        energies_vec[band] = e_flat[sort_order]
        chern_vec[band] = c_flat[sort_order]
    end

    return plt_cum_chern_vs_energy(energies_vec, chern_vec; kwargs...)
end

# Overload for 2D Matrix [nbands × N_steps]
function plt_cum_chern_vs_energy(
    accumulation_energies::Matrix{Float64},
    cumulative_chern::Matrix{Float64};
    kwargs...
)
    energies_vec = [accumulation_energies[i, :] for i in 1:size(accumulation_energies, 1)]
    chern_vec = [cumulative_chern[i, :] for i in 1:size(cumulative_chern, 1)]
    return plt_cum_chern_vs_energy(energies_vec, chern_vec; kwargs...)
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
Helper function to extract accumulated Chern number at a specific Fermi energy
"""
function get_accumulated_chern_at_energy(band_data::NamedTuple, E_F::Float64)
    """
    Extracts C(E_F) for each band at the specified Fermi energy.
    
    Returns: (c_band1, c_band2, c_total)
        - c_band1: accumulated Chern at E_F for lower band
        - c_band2: accumulated Chern at E_F for upper band  
        - c_total: total accumulated Chern at E_F
    """
    nbands = size(band_data.cum_chern_per_band, 1)
    
    chern_at_EF = Float64[]
    
    for band in 1:nbands
        # Flatten the 3D arrays
        E_flat = vec(band_data.plaquette_energies[band, :, :])
        C_flat = vec(band_data.cum_chern_per_band[band, :, :])
        
        # Sort by energy
        sort_idx = sortperm(E_flat)
        E_sorted = E_flat[sort_idx]
        C_sorted = C_flat[sort_idx]
        
        # Interpolate C at E_F using linear interpolation
        if E_F < E_sorted[1] || E_F > E_sorted[end]
            @warn "Energy E_F=$E_F is outside the band energy range [$(E_sorted[1]), $(E_sorted[end])] for band $band"
            # Extrapolate with nearest value
            C_at_EF = if E_F < E_sorted[1]
                C_sorted[1]
            else
                C_sorted[end]
            end
        else
            # Find the surrounding indices
            idx_upper = searchsortedfirst(E_sorted, E_F)
            idx_lower = idx_upper - 1
            
            if idx_lower < 1
                C_at_EF = C_sorted[1]
            elseif idx_upper > length(E_sorted)
                C_at_EF = C_sorted[end]
            else
                # Linear interpolation
                E_lo = E_sorted[idx_lower]
                E_hi = E_sorted[idx_upper]
                C_lo = C_sorted[idx_lower]
                C_hi = C_sorted[idx_upper]
                
                if abs(E_hi - E_lo) < 1e-14
                    C_at_EF = C_lo
                else
                    C_at_EF = C_lo + (C_hi - C_lo) * (E_F - E_lo) / (E_hi - E_lo)
                end
            end
        end
        
        push!(chern_at_EF, C_at_EF)
    end
    
    # Total accumulated Chern (sum of both bands)
    c_total = sum(chern_at_EF)
    
    return chern_at_EF[1], chern_at_EF[2], c_total
end


"""
Plotting function for accumulated Chern number vs orbital displacement d at fixed Fermi energy
"""
function plt_accumulated_chern_vs_d_at_EF(
    data_folder::String,
    ds::Vector{Float64};
    A::Float64=1.0,
    B::Float64=1.0,
    m::Float64=-1.0,
    gamma::Float64=0.5,
    distortion_type::Symbol=:symmetric,
    phi::Float64=0.0,
    E_F::Float64=1.1,
    ylims::Union{Nothing, Tuple{Float64, Float64}}=nothing,
    title::Union{String, LaTeXString}="Accumulated Chern Number at \$E_F = $E_F\$ vs Orbital Displacement \$d\$"
)
    """
    Plots C(E_F) vs d for both bands and their total.
    
    Parameters:
        - data_folder: path to folder containing bandstrucutre_data subdirectory
        - ds: vector of orbital displacement values
        - A, B, m, gamma, distortion_type, phi: parameters matching data generation
        - E_F: Fermi energy at which to evaluate accumulated Chern
        - title: plot title
    
    Returns: Plot object
    """
    
    data_path = joinpath(data_folder, "bandstrucutre_data")
    
    # Buffers for results
    C_band1_vals = Float64[]
    C_band2_vals = Float64[]
    C_total_vals = Float64[]
    valid_ds = Float64[]
    
    # Load data for each d value
    for d in ds
        # Build filename matching the data generation pattern
        filename = "bulk_band_berry_data_A$(A)_B$(B)_m$(m)_gamma$(gamma)_perturbation_$(distortion_type)_phi$(phi)_d$(d).jld2"
        filepath = joinpath(data_path, filename)
        
        # Check if file exists
        if !isfile(filepath)
            @warn "File not found: $filepath, skipping d=$d"
            continue
        end
        
        # Load data
        @load filepath data
        band_data = data.band
        
        # Extract C(E_F) for this d value
        try
            c1, c2, c_tot = get_accumulated_chern_at_energy(band_data, E_F)
            push!(valid_ds, d)
            push!(C_band1_vals, c1)
            push!(C_band2_vals, c2)
            push!(C_total_vals, c_tot)
        catch e
            @warn "Error processing d=$d: $e"
            continue
        end
    end
    
    if ylims !== nothing
        y_min, y_max = ylims
    else
        y_min = minimum(vcat(C_band1_vals, C_band2_vals, C_total_vals)) - 0.1
        y_max = maximum(vcat(C_band1_vals, C_band2_vals, C_total_vals)) + 0.1
    end

    # Create the plot
    p = plot(
        xlabel=L"d \mathrm{ (orbital displacement)}",
        ylabel=L"\mathcal{C}(E_F)",
        title=title,
        legend=:right,
        grid=true,
        linewidth=2.5,
        marker=:circle,
        markersize=4,
        ylims=(y_min, y_max),
    )
    
    # Plot the three lines
    plot!(p, valid_ds, C_band1_vals, 
        label=L"\mathrm{Lower band } \ \mathcal{C}_1(E_F)",
        color=:royalblue,
        linewidth=2.5,
        # marker=:circle
    )
    
    plot!(p, valid_ds, C_band2_vals,
        label=L"\mathrm{Upper band } \ \mathcal{C}_2(E_F)",
        color=:crimson,
        linewidth=2.5,
        # marker=:square
    )
    
    plot!(p, valid_ds, C_total_vals,
        label=L"\mathrm{Total } \mathcal{C}_{\mathrm{tot}}(E_F)",
        color=:darkorange,
        linewidth=2.5,
        # marker=:diamond,
        linestyle=:dash
    )
    
    return p
end

############################################################################
"""
Helper function to extract energy from 2D accumulated Chern landscape
"""
function get_energy_at_chern_2d(energies_2d::AbstractMatrix, chern_2d::AbstractMatrix, target_chern::Float64; tol::Float64=0.02)
    """
    Extract energy from 2D C(E,k) heatmap by finding contour points where C ≈ target_chern.
    This is the correct approach since accumulated Chern is a 2D landscape in k-space,
    not a simple function of energy.
    
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

############################################################################
"""
Plotting function for energy at critical Chern values vs orbital displacement d
"""
function plt_energy_at_critical_chern_vs_d(
    data_folder::String,
    ds::Vector{Float64};
    A::Float64=1.0,
    B::Float64=1.0,
    m::Float64=-1.0,
    gamma::Float64=0.5,
    distortion_type::Symbol=:symmetric,
    phi::Float64=0.0,
    ylims::Union{Nothing, Tuple{Float64, Float64}}=nothing,
    title::Union{String, LaTeXString}="Energy at Critical Chern Values vs Orbital Displacement \$d\$",
    plt_vs_inverse_d::Bool=false
)
    """
    Plots E(C=±0.5) vs d for both bands.
    
    Parameters:
        - data_folder: path to folder containing bandstrucutre_data subdirectory
        - ds: vector of orbital displacement values
        - A, B, m, gamma, distortion_type, phi: parameters matching data generation
        - title: plot title
    
    Returns: Plot object
    """
    
    data_path = joinpath(data_folder, "bandstrucutre_data")
    
    # Buffers for results
    E_band1_pos = Float64[]   # Energy at C = +0.5 for band 1
    E_band1_neg = Float64[]   # Energy at C = -0.5 for band 1
    E_band2_pos = Float64[]   # Energy at C = +0.5 for band 2
    E_band2_neg = Float64[]   # Energy at C = -0.5 for band 2
    valid_ds = Float64[]
    
    # Load data for each d value
    for d in ds
        # Build filename matching the data generation pattern
        filename = "bulk_band_berry_data_A$(A)_B$(B)_m$(m)_gamma$(gamma)_perturbation_$(distortion_type)_phi$(phi)_d$(d).jld2"
        filepath = joinpath(data_path, filename)
        
        # Check if file exists
        if !isfile(filepath)
            continue
        end
        
        # Load data
        @load filepath data
        band_data = data.band
        
        # Extract critical energies for this d value
        try
            # Process each band
            nbands = size(band_data.cum_chern_per_band, 1)
            
            for band in 1:nbands
                # Use 2D contour approach: find all k-points where C ≈ target
                # (critical Chern is a 2D landscape, not a 1D function of E)
                energies_2d = band_data.plaquette_energies[band, :, :]
                chern_2d = band_data.cum_chern_per_band[band, :, :]
                
                # Get energy at C = +0.5 (find contour points and take median)
                E_at_pos_half = get_energy_at_chern_2d(energies_2d, chern_2d, 0.5; tol=0.02)
                
                # Get energy at C = -0.5
                E_at_neg_half = get_energy_at_chern_2d(energies_2d, chern_2d, -0.5; tol=0.02)
                
                if band == 1
                    push!(E_band1_pos, E_at_pos_half)
                    push!(E_band1_neg, E_at_neg_half)
                elseif band == 2
                    push!(E_band2_pos, E_at_pos_half)
                    push!(E_band2_neg, E_at_neg_half)
                end
            end
            
            push!(valid_ds, d)
            
        catch e
            @warn "Error processing d=$d: $e"
            continue
        end
    end
    
    if ylims !== nothing
        y_min, y_max = ylims
    else
        all_energies = vcat(
            filter(!isnan, E_band1_pos), 
            filter(!isnan, E_band1_neg),
            filter(!isnan, E_band2_pos), 
            filter(!isnan, E_band2_neg)
        )
        if !isempty(all_energies)
            y_min = minimum(all_energies) - 0.2
            y_max = maximum(all_energies) + 0.2
        else
            y_min, y_max = -2.0, 2.0
        end
    end

    if plt_vs_inverse_d
        xlabel=L"1/d \mathrm{ (inverse orbital displacement)}"
    else
        xlabel=L"d \mathrm{ (orbital displacement)}"
    end

    # Create the plot
    p = plot(
        xlabel=xlabel,
        ylabel=L"E \mathrm{ (critical Chern)}",
        title=title,
        legend=:best,
        grid=true,
        linewidth=2.5,
        marker=:circle,
        markersize=4,
        ylims=(y_min, y_max),
    )
    
    if plt_vs_inverse_d
        plot!(p, 1.0 ./ valid_ds, E_band1_pos,
            label=L"\mathrm{Band\ 1:\ } \mathcal{C} = +0.5",
            color=:royalblue,
            linewidth=2.5,
        )
        
        plot!(p, 1.0 ./ valid_ds, E_band1_neg,
            label=L"\mathrm{Band\ 1:\ } \mathcal{C} = -0.5",
            color=:royalblue,
            linewidth=2.5,
            linestyle=:dash,
        )
        
        plot!(p, 1.0 ./ valid_ds, E_band2_pos,
            label=L"\mathrm{Band\ 2:\ } \mathcal{C} = +0.5",
            color=:crimson,
            linewidth=2.5,
        )
        
        plot!(p, 1.0 ./ valid_ds, E_band2_neg,
            label=L"\mathrm{Band\ 2:\ } \mathcal{C} = -0.5",
            color=:crimson,
            linewidth=2.5,
            linestyle=:dash,
        )
    else
    
        # Plot the four lines
        plot!(p, valid_ds, E_band1_pos,
            label=L"\mathrm{Band\ 1:\ } \mathcal{C} = +0.5",
            color=:royalblue,
            linewidth=2.5,
        )
        
        plot!(p, valid_ds, E_band1_neg,
            label=L"\mathrm{Band\ 1:\ } \mathcal{C} = -0.5",
            color=:royalblue,
            linewidth=2.5,
            linestyle=:dash,
        )
        
        plot!(p, valid_ds, E_band2_pos,
            label=L"\mathrm{Band\ 2:\ } \mathcal{C} = +0.5",
            color=:crimson,
            linewidth=2.5,
        )
        
        plot!(p, valid_ds, E_band2_neg,
            label=L"\mathrm{Band\ 2:\ } \mathcal{C} = -0.5",
            color=:crimson,
            linewidth=2.5,
            linestyle=:dash,
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
ms = [-1.0]
gammas = [2.0] #collect(0.0:0.01:2.0)
distortion_type = :sym_cos_sum
phis = [pi/4] #collect(0.0:0.25π:2π)
ds = collect(0.0:1.0:20.0)
Nkx = 201
Nky = 201

data_folder = joinpath("data", "bulk_band_berry_data", "embedding")
plot_folder = joinpath("plots", "bulk_band_berry_data", "embedding")
isdir(data_folder) || mkpath(data_folder)
isdir(plot_folder) || mkpath(plot_folder)
println("Plot folder: $(plot_folder)")

generate_new_data = true 
generate_new_plots = false
generate_gamma_scatter = false
generate_cum_chern_extrema = false
generate_cum_chern_vs_d_at_EF = false
generate_energy_at_critical_chern_vs_d = true

if generate_new_data
    println("Generating new data...")
    @showprogress for (A, B, m, gamma, phi, d) in Iterators.product(As, Bs, ms, gammas, phis, ds)
        band_data = compute_bulk_band_berry_data(
            Nkx=Nkx,
            Nky=Nky,
            A=A,
            B=B,
            m=m,
            gamma=gamma,
            phi=phi,
            orbital_displacement=d,
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

        
        euler_char_data = compute_euler_characteristic_vs_energy(
            band_data;
            E_mesh=range(-0.0, 3.0, length=501)
        )

        hasproperty(euler_char_data, :E_mesh) || error("Euler payload missing E_mesh")
        hasproperty(euler_char_data, :euler_char) || error("Euler payload missing euler_char")
        length(euler_char_data.E_mesh) == length(euler_char_data.euler_char) ||
            error("Euler payload length mismatch: E_mesh and euler_char")

        data = (
            band = band_data,
            edge = edge_state_data,
            euler_char = euler_char_data,
            euler_char_schema_version = 2
        ) 

        save_path = joinpath(data_folder, "bandstrucutre_data")
        isdir(save_path) || mkpath(save_path)
        @save joinpath(save_path, "bulk_band_berry_data_A$(A)_B$(B)_m$(m)_gamma$(gamma)_perturbation_$(distortion_type)_phi$(phi)_d$(d).jld2") data

        # debug_accumulation_extrema(data)
    end
else
    println("Using existing data...")
end



if generate_new_plots
    data_to_unpack = joinpath(data_folder, "bandstrucutre_data")
    println("Plotting data from $(data_to_unpack)...")

    @showprogress for (A, B, m, gamma, phi, d) in Iterators.product(As, Bs, ms, gammas, phis, ds)
        # 1. Build the target filename directly
        filename = "bulk_band_berry_data_A$(A)_B$(B)_m$(m)_gamma$(gamma)_perturbation_$(distortion_type)_phi$(phi)_d$(d).jld2"
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

        # 5. Plot the Berry curvature heatmaps. F_xy flux per plaquette resolved for (kx, ky) for both bands.
        F_xy_heatmaps = plt_k_resolved_F_xy_heatmaps(
            data.band.kx_vals,
            data.band.ky_vals,
            data.band.berry_curvature;
            title = L"\gamma = %$(gamma), F_{xy} / 2\pi, \phi = %$(phi), d = %$(d)",
            # colour = :RdBu,
            # share_colour_scale = false
        )

        savefig(
            F_xy_heatmaps, 
            joinpath(plot_folder, "bulk_band_dEdC_A$(A)_B$(B)_m$(m)_gamma$(gamma)_perturbation_$(distortion_type)_phi$(phi)_d$(d).png")
        )

        # 6. Plot the accumulated chern number as a heatmap on (kx, ky) grid
        cum_chern_heatmaps = plt_accumulated_chern_heatmaps(
            data.band.kx_vals,
            data.band.ky_vals,
            data.band.cum_chern_per_band;
            title = L"\gamma = %$(gamma), \mathcal{C}(k_x, k_y), \phi = %$(phi), d = %$(d)",
            colour = :viridis
        )

        savefig(
            cum_chern_heatmaps, 
            joinpath(plot_folder, "bulk_band_accumulated_chern_A$(A)_B$(B)_m$(m)_gamma$(gamma)_perturbation_$(distortion_type)_phi$(phi)_d$(d).png")
        )

        # # 6.5 Plot the cum chern with inflection point contours
        # cum_chern_inflection = plt_accumulated_chern_2d_inflections(
        #     data.band.kx_vals,
        #     data.band.ky_vals,
        #     data.band.cum_chern_per_band;
        #     title = L"\gamma = %$(gamma), \mathcal{C}(k_x, k_y) \text{ with Inflection Contours}",
        #     colour = :viridis,
        #     smooth_window = 10, #Int(round(0.1 * (Nkx * Nky))),
        #     min_slope_tol = 1e-5
        # )

        # savefig(
        #     cum_chern_inflection, 
        #     joinpath(plot_folder, "bulk_band_accumulated_chern_2d_inflections_A$(A)_B$(B)_m$(m)_gamma$(gamma)_perturbation_$(distortion_type).png")
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

        # # 8. Plot the heatmap of accumulated berry curvature given overlapping bands at fermi energy E_cut
        # energies = data.band.plaquette_energies
        # E_cuts = collect(range(minimum(energies), stop=maximum(energies), length=50))

        # # Pre-compute global clims across ALL E_cuts for this specific gamma
        # n_bands, nkx, nky = size(energies)
        # max_abs_gamma = 0.0

        # for E_cut in E_cuts
        #     cum_2d = zeros(Float64, nkx, nky)
        #     for n in 1:n_bands, i in 1:nkx, j in 1:nky
        #         if energies[n, i, j] <= E_cut
        #             cum_2d[i, j] += data.band.berry_curvature[n, i, j] / (2π)
        #         end
        #     end
        #     max_abs_gamma = max(max_abs_gamma, maximum(abs, cum_2d))
        # end

        # # Set symmetric bounds around zero
        # global_clims = (-max_abs_gamma, max_abs_gamma)

        # for (idx, E_cut) in enumerate(E_cuts)

        #     cum_chern_E_cut_heatmap = plt_joint_band_accumulated_chern_heatmap(
        #         data.band.kx_vals,
        #         data.band.ky_vals,
        #         data.band.berry_curvature,
        #         data.band.plaquette_energies,
        #         E_cut,
        #         gamma;
        #         global_clims=global_clims,
        #         colour = :RdBu
        #     )

        #     savefig(
        #         cum_chern_E_cut_heatmap, 
        #         joinpath(plot_folder, "$(idx)_bulk_band_accumulated_chern_Ecut_$(round(E_cut, digits=2))_A$(A)_B$(B)_m$(m)_gamma$(gamma)_perturbation_$(distortion_type).png")
        #     )
        # end


        # # 9. Plot the accumulated chern number band-resolved vs energyy (no k-resolution)
        # cum_chern_separate_band_vs_energy_plot = plt_cum_chern_vs_energy(
        #     data.band.plaquette_energies,
        #     data.band.cum_chern_per_band;
        #     title = L"\gamma = %$(gamma), \mathcal{C}(E)",
        #     # ylims = (-0.2, 1.01)
        # )

        # savefig(
        #     cum_chern_separate_band_vs_energy_plot, 
        #     joinpath(plot_folder, "bulk_band_accumulated_chern_vs_energy_A$(A)_B$(B)_m$(m)_gamma$(gamma)_perturbation_$(distortion_type).png")
        # )



    end
    
end



if generate_gamma_scatter
    println("Generating gamma scatter plot...")
    gamma_data_list = Tuple{Float64, Any}[]
    gamma_data_list_with_euler = Tuple{Float64, Any}[]  # For precalc euler function
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
            push!(gamma_data_list_with_euler, (gamma, data))  # Store full data for precalc function
        end
    end

    if !isempty(gamma_data_list)
        m_val = ms[1]
        A_val = As[1]
        B_val = Bs[1]

        # # 1. Plot the energy at which the accumulated Chern number reaches 0.5 and 0.0 for each gamma value
        # band_separated_plt = plt_accumulated_chern_energy_vs_gamma(
        #     gamma_data_list; 
        #     target_p = 0.5, 
        #     target_q = 0.05,
        #     tol = 0.005
        # )

        # savefig(band_separated_plt, joinpath(plot_folder, "energy_at_accumulated_chern_0.5_vs_gamma$(minimum(gammas))-$(maximum(gammas))_A$(A_val)_B$(B_val)_m$(m_val).png"))


        # # 2. Plot the energy at which the accumulated Chern number reaches 0.5 and 0.0 for each gamma value, but only for the joint band (i.e., considering all bands together)
        # band_joint_plot = plt_joint_band_accumulated_chern_energy_vs_gamma(
        #     gamma_data_list; 
        #     target_p = 0.5, 
        #     target_q = 0.0,
        #     tol = 0.005
        # )

        # savefig(band_joint_plot, joinpath(plot_folder, "energy_at_joint_band_accumulated_chern_0.5_vs_gamma$(minimum(gammas))-$(maximum(gammas))_A$(A_val)_B$(B_val)_m$(m_val).png"))


        # # 3. Plot the energy at which the accumulated Chern number reaches its extrema (maxima and minima) for each gamma value, but only for the joint band (i.e., considering all bands together)
        # band_extrema_plot = plt_joint_band_accumulated_chern_extrema_vs_gamma(
        #     gamma_data_list; 
        #     tol = 0.005
        # )

        # savefig(band_extrema_plot, joinpath(plot_folder, "energy_at_joint_band_accumulated_chern_extrema_vs_gamma$(minimum(gammas))-$(maximum(gammas))_A$(A_val)_B$(B_val)_m$(m_val).png"))


        # # 4. Plot the energy at which the accumulated Chern number reaches specific smart chern points for each gamma value, but only for the joint band (i.e., considering all bands together), with smart detection of flat tops and bottoms, and marking of zero crossings
        # joint_tracking_plot = plt_joint_band_smart_extrema_accumulated_chern_energy_vs_gamma(
        #     gamma_data_list; 
        #     peak_tol = 0.01,
        #     e_cluster_tol = 0.1,
        #     flat_span_tol = 0.02,
        #     zero_tol = 1e-3
        # )

        # savefig(joint_tracking_plot, joinpath(plot_folder, "energy_at_joint_band_accumulated_chern_extrema_and_zero_vs_gamma$(minimum(gammas))-$(maximum(gammas))_A$(A_val)_B$(B_val)_m$(m_val).png"))


        # # 5. 
        # euler_map_plot = plt_euler_characteristic_heatmap_vs_gamma(
        #     gamma_data_list; 
        #     n_energies = 300,
        #     energy_lims = nothing, # (-2.0,-3.0), # nothing, #(-2.5, 2.5),
        #     periodic = true,
        #     add_contours = false
        # )

        # savefig(euler_map_plot, joinpath(plot_folder, "euler_characteristic_heatmap_vs_gamma$(minimum(gammas))-$(maximum(gammas))_A$(A_val)_B$(B_val)_m$(m_val).png"))

        # # 5b. Using pre-calculated Euler characteristic data
        # euler_map_plot_precalc = plt_euler_characteristic_heatmap_vs_gamma_precalc(
        #     gamma_data_list_with_euler; 
        #     add_contours = false
        # )

        # savefig(euler_map_plot_precalc, joinpath(plot_folder, "euler_characteristic_heatmap_vs_gamma_precalc$(minimum(gammas))-$(maximum(gammas))_A$(A_val)_B$(B_val)_m$(m_val).png"))



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


if generate_cum_chern_vs_d_at_EF
    println("Generating cumulative chern vs orbital displacement at fixed Fermi energy plots...")
    subfolder = joinpath(plot_folder, "cum_chern_vs_d_at_EF")
    isdir(subfolder) || mkpath(subfolder)

    data_to_unpack = joinpath(data_folder, "bandstrucutre_data", "embedding")


    E_F_value = 1.0  # Fixed Fermi energy (can be changed here)
    distortion_type = :sym_cos_sum
    cum_chern_vs_d_plot = plt_accumulated_chern_vs_d_at_EF(
        data_folder,
        ds;
        A=1.0,
        B=1.0,
        m=-1.0,
        gamma=2.0,
        distortion_type=distortion_type,
        phi=0.0,
        E_F=E_F_value,
        # ylims=(-0.2, 1.05),
        title=L"Accumulated Chern Number at $E_F = %$(E_F_value)$ vs Orbital Displacement $d$"
    )

    savefig(
        cum_chern_vs_d_plot,
        joinpath(plot_folder, "accumulated_chern_vs_d_at_EF_$(E_F_value)_IndType_$(distortion_type).png")
    )
    println("Plot saved to: $(joinpath(plot_folder, "accumulated_chern_vs_d_at_EF_$(E_F_value)_IndType_$(distortion_type).png"))")
end


if generate_energy_at_critical_chern_vs_d
    println("Generating energy at critical Chern values vs orbital displacement plots...")
    subfolder = joinpath(plot_folder, "energy_at_critical_chern_vs_d")
    isdir(subfolder) || mkpath(subfolder)

    distortion_type = :asym_sin_sum
    energy_at_crit_chern_plot = plt_energy_at_critical_chern_vs_d(
        data_folder,
        ds;
        A=1.0,
        B=1.0,
        m=-1.0,
        gamma=2.0,
        distortion_type=distortion_type,
        phi=phis[1],
        title=L"Energy at Critical Chern Values ($\mathcal{C} = \pm 0.5$) vs Orbital Displacement $d$",
        plt_vs_inverse_d = true
    )

    savefig(
        energy_at_crit_chern_plot,
        joinpath(subfolder, "energy_at_critical_chern_vs_d_IndType_$(distortion_type)_phi$(phis[1]).png")
    )
    println("Plot saved to: $(joinpath(subfolder, "energy_at_critical_chern_vs_d_IndType_$(distortion_type)_phi$(phis[1]).png"))")
end