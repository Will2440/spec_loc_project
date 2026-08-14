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

println("Julia Threads: ", Threads.nthreads())
println("BLAS Threads:  ", LinearAlgebra.BLAS.get_num_threads())


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
    B_y::Real=1.0,
    perturbation_type::Symbol=:none,
    winding_number::Int=0
)::Matrix{ComplexF64}

    winding_number >= 0 || error("winding_number must be non-negative.")

    # Base QWZ Hamiltonian in k-space
    d_x = A * sin(kx)
    d_y = A * sin(ky)
    d_z = m + 2B - B * cos(kx) - B * B_y * cos(ky)
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

## choose which directions to perturb in independently
function ribbon_perturbed_hamiltonian_qwz(
    ky::Real,
    Lx::Int;
    A::Real=1.0,
    B::Real=1.0,
    m::Real=0.0,
    gamma::Real=0.0,
    B_y::Real=0.0,
    perturbation_type::Symbol=:none,
    perturb_x::Bool=false, # standard don't perturn the OBC direction
    perturb_y::Bool=true,
    winding_number::Int=0
)
    Lx > 0 || error("Lx must be positive.")

    if winding_number > 0
        return ribbon_edge_winding_hamiltonian_qwz(
            ky, Lx;
            A=A, B=B, m=m, gamma=gamma, B_y=B_y,
            perturbation_type=perturbation_type,
            perturb_x=perturb_x,
            perturb_y=perturb_y,
            winding_number=winding_number,
        )
    end

    # 1. Construct x-hopping block (tx)
    # Base QWZ hopping along x
    tx = ComplexF64.(-0.5 * B * sigma_z - 0.5im * A * sigma_x)

    # Apply x-perturbation only if perturb_x is true
    if perturb_x
        if perturbation_type == :symmetric
            # gamma * cos(kx) -> hopping term 0.5 * gamma * identity
            tx += ComplexF64.(0.5 * gamma * identity)
        elseif perturbation_type == :tilt
            # gamma * sin(kx) -> hopping term 0.5im * gamma * identity
            tx += ComplexF64.(0.5im * gamma * identity)
        end
    end

    # 2. Construct ky-dependent on-site block
    # Base QWZ on-site term for fixed ky
    onsite_ky = (m + 2.0 * B - B * B_y * cos(ky)) * sigma_z + A * sin(ky) * sigma_y

    # Apply y-perturbation only if perturb_y is true
    if perturb_y
        if perturbation_type == :symmetric
            # gamma * cos(ky) term
            onsite_ky += gamma * cos(ky) * identity
        elseif perturbation_type == :tilt
            # Tilt is sin(kx), has no ky component
        end
    end

    # 3. Assemble 2Lx x 2Lx block tridiagonal ribbon matrix
    dim = 2 * Lx
    H = zeros(ComplexF64, dim, dim)

    for x in 1:Lx
        idx = (2x - 1):(2x)
        H[idx, idx] = onsite_ky
        
        if x < Lx
            next_idx = (2x + 1):(2x + 2)
            H[idx, next_idx] = tx
            H[next_idx, idx] = tx'
        end
    end

    return Hermitian(H)
end

# Edge winding specialized implementation
# Exact block-Fourier transform of the k_x-dependent phase rotation
function ribbon_edge_winding_hamiltonian_qwz(
    ky::Real,
    Lx::Int;
    A::Real=1.0,
    B::Real=1.0,
    m::Real=0.0,
    gamma::Real=0.0,
    B_y::Real=0.0,
    perturbation_type::Symbol=:none,
    perturb_x::Bool=false,
    perturb_y::Bool=true,
    winding_number::Int=1
)
    """
    Build ribbon Hamiltonian with edge_winding perturbation.

    The distortion is implemented exactly as
        d_x'(k_x, k_y) + i d_y'(k_x, k_y) = exp(i n k_x) * (sin(k_x) + i sin(k_y))

    For a fixed k_y, this means the off-diagonal entry is
        g_n(k_x, k_y) = exp(i n k_x) * (sin(k_x) + i sin(k_y))

    and the ribbon Hamiltonian is obtained from the exact Fourier coefficients
    of g_n and g_n^* in the x direction. This produces hoppings at ranges
    n + 1, n, and n - 1, with the n - 1 term collapsing to an on-site
    contribution when n = 1.
    """
    n = winding_number
    n >= 0 || error("winding_number must be non-negative.")

    if n == 0
        return ribbon_perturbed_hamiltonian_qwz(
            ky,
            Lx;
            A=A,
            B=B,
            m=m,
            gamma=gamma,
            B_y=B_y,
            perturbation_type=perturbation_type,
            perturb_x=perturb_x,
            perturb_y=perturb_y,
            winding_number=0,
        )
    end

    dim = 2 * Lx
    H = zeros(ComplexF64, dim, dim)

    function add_block!(block::AbstractMatrix{<:Number}, displacement::Int)
        displacement >= 0 || error("displacement must be non-negative.")

        if displacement == 0
            for x in 1:Lx
                idx = (2x - 1):(2x)
                H[idx, idx] .+= ComplexF64.(block)
            end
            return
        end

        displacement < Lx || return

        for x in 1:(Lx - displacement)
            idx = (2x - 1):(2x)
            shifted_idx = (2 * (x + displacement) - 1):(2 * (x + displacement))
            H[idx, shifted_idx] .+= ComplexF64.(block)
            H[shifted_idx, idx] .+= ComplexF64.(block')
        end
    end

    sin_ky = sin(ky)

    # On-site block from the k_x-independent part of d_z.
    onsite_ky = (m + 2.0 * B - B * B_y * cos(ky)) * sigma_z
    for x in 1:Lx
        idx = (2x - 1):(2x)
        H[idx, idx] .= onsite_ky
    end

    # -B cos(k_x) sigma_z
    add_block!(ComplexF64.(-0.5 * B * sigma_z), 1)

    # e^{i(n+1)k_x} / (2i) contribution to the lower-left entry.
    add_block!(ComplexF64.([0.0 + 0.0im 0.0 + 0.0im; -0.5im * A 0.0 + 0.0im]), n + 1)

    # -i sin(k_y) e^{i n k_x} contribution to the lower-left entry.
    if n == 0
        add_block!(ComplexF64.(A * sin_ky * sigma_y), 0)
    else
        add_block!(ComplexF64.([0.0 + 0.0im 0.0 + 0.0im; 1im * A * sin_ky 0.0 + 0.0im]), n)
    end

    # -e^{i(n-1)k_x} / (2i) contribution to the lower-left entry.
    if n == 1
        add_block!(ComplexF64.(0.5 * A * sigma_y), 0)
    elseif n > 1
        add_block!(ComplexF64.([0.0 + 0.0im 0.0 + 0.0im; 0.5im * A 0.0 + 0.0im]), n - 1)
    end

    if perturbation_type == :symmetric
        if perturb_x && gamma != 0.0
            add_block!(ComplexF64.(0.5 * gamma * identity), 1)
        end

        if perturb_y && gamma != 0.0
            for x in 1:Lx
                idx = (2x - 1):(2x)
                H[idx, idx] .+= ComplexF64.(gamma * cos(ky) * identity)
            end
        end
    elseif perturbation_type == :tilt
        if perturb_x && gamma != 0.0
            add_block!(ComplexF64.(0.5im * gamma * identity), 1)
        end
    elseif perturbation_type != :none
        error("Unknown perturbation type: $perturbation_type")
    end
    
    return Hermitian(H)
end

function ribbon_dH_dky_qwz(
    ky::Real,
    Lx::Int;
    A::Real=1.0,
    B::Real=1.0,
    m::Real=0.0,
    gamma::Real=0.0,
    B_y::Real=0.0,
    perturbation_type::Symbol=:none,
    perturb_x::Bool=false,
    perturb_y::Bool=true,
    winding_number::Int=0
)::Matrix{ComplexF64}
    Lx > 0 || error("Lx must be positive.")

    dim = 2 * Lx
    dH = zeros(ComplexF64, dim, dim)

    function add_block!(block::AbstractMatrix{<:Number}, displacement::Int)
        displacement >= 0 || error("displacement must be non-negative.")

        if displacement == 0
            for x in 1:Lx
                idx = (2x - 1):(2x)
                dH[idx, idx] .+= ComplexF64.(block)
            end
            return
        end

        displacement < Lx || return

        for x in 1:(Lx - displacement)
            idx = (2x - 1):(2x)
            shifted_idx = (2 * (x + displacement) - 1):(2 * (x + displacement))
            dH[idx, shifted_idx] .+= ComplexF64.(block)
            dH[shifted_idx, idx] .+= ComplexF64.(block')
        end
    end

    if winding_number == 0
        onsite_derivative = B * B_y * sin(ky) * sigma_z + A * cos(ky) * sigma_y

        if perturb_y
            if perturbation_type == :symmetric
                onsite_derivative += -gamma * sin(ky) * identity
            elseif perturbation_type != :tilt && perturbation_type != :none
                error("Unknown perturbation type: $perturbation_type")
            end
        elseif perturbation_type != :symmetric && perturbation_type != :tilt && perturbation_type != :none
            error("Unknown perturbation type: $perturbation_type")
        end

        add_block!(ComplexF64.(onsite_derivative), 0)
        return dH
    end

    add_block!(ComplexF64.(B * B_y * sin(ky) * sigma_z), 0)
    add_block!(ComplexF64.([0.0 + 0.0im 0.0 + 0.0im; 1im * A * cos(ky) 0.0 + 0.0im]), winding_number)

    if perturb_y
        if perturbation_type == :symmetric
            add_block!(ComplexF64.(-gamma * sin(ky) * identity), 0)
        elseif perturbation_type != :tilt && perturbation_type != :none
            error("Unknown perturbation type: $perturbation_type")
        end
    elseif perturbation_type != :symmetric && perturbation_type != :tilt && perturbation_type != :none
        error("Unknown perturbation type: $perturbation_type")
    end

    return dH
end

function compute_ribbon_spectrum(
    Lx::Int; 
    N_ky::Int=101, 
    kwargs...
)
    ky_vals = range(-π, π, length=N_ky)
    energies = zeros(Float64, 2 * Lx, N_ky)
    
    for (i, ky) in enumerate(ky_vals)
        H_ribbon = ribbon_perturbed_hamiltonian_qwz(ky, Lx; kwargs...)
        energies[:, i] = eigvals(H_ribbon)
    end
    
    return ky_vals, energies
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
    k_space_kwargs = (; [k => v for (k, v) in pairs(kwargs) if k ∈ (:A, :B, :m, :gamma, :B_y, :perturbation_type, :winding_number)]...)

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

    flat_energies = vec(plaquette_energies)
    flat_weights = vec(berry_curvature) ./ (2π)
    order = sortperm(flat_energies)

    accumulation_energies = flat_energies[order]
    accumulation_weights = flat_weights[order]
    cumulative_chern = cumsum(flat_weights[order])

    return (
        kx_vals=kx_vals,
        ky_vals=ky_vals,
        energies=energies,
        berry_curvature=berry_curvature,
        plaquette_energies=plaquette_energies,
        band_chern_numbers=band_chern_numbers,
        accumulation_energies=accumulation_energies,
        accumulation_weights=accumulation_weights,
        cumulative_chern=cumulative_chern,
    )
end

function compute_state_ipr(state::AbstractVector{<:Number})
    length(state) % 2 == 0 || error("State vector must have an even number of components.")

    probabilities = abs2.(state)
    site_probabilities = @views probabilities[1:2:end] .+ probabilities[2:2:end]

    return sum(site_probabilities .^ 2)
end

function compute_state_edge_weight(state::AbstractVector{<:Number}, edge_sites::Int=3)
    length(state) % 2 == 0 || error("State vector must have an even number of components.")

    site_probabilities = @views abs2.(state[1:2:end]) .+ abs2.(state[2:2:end])
    nsites = length(site_probabilities)
    edge_sites_clamped = clamp(edge_sites, 1, nsites)

    left_weight = sum(@view site_probabilities[1:edge_sites_clamped])
    right_weight = sum(@view site_probabilities[(nsites - edge_sites_clamped + 1):nsites])

    return left_weight + right_weight
end

function compute_gap_closure_ky_extent(
    ky_vals,
    energies,
    edge_weights;
    edge_weight_threshold::Float64=0.6,
    gap_tol::Float64=0.05,
)
    """
    Determine the periodic k_y extent where at least one edge-localized state
    remains spectrally inside the bulk gap by more than gap_tol.

    This measures how far the winding edge state stretches in k_y before it
    merges into the bulk bands.
    """
    N_ky = length(ky_vals)
    nbands = size(energies, 1)
    size(edge_weights) == size(energies) || error("edge_weights and energies must have the same shape.")

    # Edge-in-gap margin at each k_y:
    # max over edge-like states of min(distance to lower bulk, distance to upper bulk).
    edge_gap_margin = fill(-Inf, N_ky)

    for i in 1:N_ky
        edge_idx = findall(edge_weights[:, i] .>= edge_weight_threshold)
        bulk_idx = findall(edge_weights[:, i] .< edge_weight_threshold)

        isempty(edge_idx) && continue
        isempty(bulk_idx) && continue

        bulk_energies = @view energies[bulk_idx, i]
        local_best = -Inf

        for idx in edge_idx
            E = energies[idx, i]

            lower_candidates = bulk_energies[bulk_energies .< E]
            upper_candidates = bulk_energies[bulk_energies .> E]
            (isempty(lower_candidates) || isempty(upper_candidates)) && continue

            lower_bulk = maximum(lower_candidates)
            upper_bulk = minimum(upper_candidates)
            margin = min(E - lower_bulk, upper_bulk - E)
            local_best = max(local_best, margin)
        end

        edge_gap_margin[i] = local_best
    end

    edge_in_gap_mask = edge_gap_margin .> gap_tol

    if !any(edge_in_gap_mask)
        return (
            ky_closure_points=Float64[],
            ky_closure_bounds=Float64[],
            ky_extent_bounds=Float64[],
            extent_in_pi=0.0,
            extent_fraction=0.0,
            edge_in_gap_mask=edge_in_gap_mask,
            edge_gap_margin=edge_gap_margin,
        )
    end

    # Build linear contiguous runs, then merge first/last if the extent wraps the BZ seam.
    runs = Vector{Vector{Int}}()
    current = Int[]
    for i in 1:N_ky
        if edge_in_gap_mask[i]
            push!(current, i)
        elseif !isempty(current)
            push!(runs, current)
            current = Int[]
        end
    end
    if !isempty(current)
        push!(runs, current)
    end

    if length(runs) > 1 && edge_in_gap_mask[1] && edge_in_gap_mask[end]
        merged = vcat(runs[end], runs[1])
        middle = length(runs) > 2 ? runs[2:end-1] : Vector{Vector{Int}}()
        runs = vcat([merged], middle)
    end

    run_lengths = [length(r) for r in runs]
    dominant_idx = argmax(run_lengths)
    ky_indices = runs[dominant_idx]
    dom_start = first(ky_indices)
    dom_end = last(ky_indices)

    period = 2π
    ky_min, ky_max = extrema(ky_vals)

    function wrap_to_window(k::Float64)
        while k < ky_min
            k += period
        end
        while k > ky_max
            k -= period
        end
        return k
    end

    function interpolate_threshold(k1::Float64, g1::Float64, k2::Float64, g2::Float64)
        k2_adj = k2
        if k2_adj < k1
            k2_adj += period
        end

        if g1 == g2
            return wrap_to_window(0.5 * (k1 + k2_adj))
        end

        k_interp = k1 + (gap_tol - g1) * (k2_adj - k1) / (g2 - g1)
        return wrap_to_window(k_interp)
    end

    # Refine two boundaries by interpolation where the edge-in-gap margin crosses gap_tol.
    prev_idx = dom_start == 1 ? N_ky : (dom_start - 1)
    next_idx = dom_end == N_ky ? 1 : (dom_end + 1)

    left_bound = interpolate_threshold(
        Float64(ky_vals[prev_idx]),
        edge_gap_margin[prev_idx],
        Float64(ky_vals[dom_start]),
        edge_gap_margin[dom_start],
    )

    right_bound = interpolate_threshold(
        Float64(ky_vals[dom_end]),
        edge_gap_margin[dom_end],
        Float64(ky_vals[next_idx]),
        edge_gap_margin[next_idx],
    )

    extent_in_pi = mod(right_bound - left_bound, period) / π

    return (
        ky_closure_points=collect(ky_vals[ky_indices]),
        ky_closure_bounds=[left_bound, right_bound],
        ky_extent_bounds=[left_bound, right_bound],
        extent_in_pi=extent_in_pi,
        extent_fraction=length(ky_indices) / N_ky,
        edge_in_gap_mask=edge_in_gap_mask,
        edge_gap_margin=edge_gap_margin,
    )
end

function compute_edge_state_ky_extent(
    ky_vals,
    energies,
    edge_weights,
    iprs;
    edge_threshold::Float64=0.5,
    ipr_threshold::Float64=0.05,
    min_ky_length::Int=3
)
    """
    Find k_y extent of edge-localized states.
    Returns: (ky_min, ky_max, extent_fraction, num_edge_bands)
    """
    N_ky = length(ky_vals)
    nbands = size(energies, 1)
    
    # Identify edge-localized states at each k_y
    is_edge_localized = (edge_weights .> edge_threshold) .& (iprs .> ipr_threshold)
    
    # Find which k_y points have ANY edge-localized band
    ky_has_edge_states = vec(any(is_edge_localized, dims=1))
    
    if !any(ky_has_edge_states)
        return (ky_min=0.0, ky_max=0.0, extent_fraction=0.0, num_edge_bands=0, ky_indices=Int[])
    end
    
    # Find contiguous regions
    ky_indices = findall(ky_has_edge_states)
    
    # Check for largest contiguous block
    if length(ky_indices) >= min_ky_length
        ky_min = ky_vals[ky_indices[1]]
        ky_max = ky_vals[ky_indices[end]]
        extent_in_pi = (ky_max - ky_min) / π
        extent_fraction = length(ky_indices) / N_ky
        
        # Count average number of edge-localized bands
        num_edge_bands_avg = mean(vec(sum(is_edge_localized, dims=1)[ky_has_edge_states]))
        
        return (
            ky_min=ky_min,
            ky_max=ky_max,
            extent_fraction=extent_fraction,
            extent_in_pi=extent_in_pi,
            num_edge_bands=num_edge_bands_avg,
            ky_indices=ky_indices
        )
    else
        return (ky_min=0.0, ky_max=0.0, extent_fraction=0.0, num_edge_bands=0, ky_indices=Int[])
    end
end

function compute_berry_curvature_by_ky(berry_curvature::Array{Float64, 3}, ky_indices::AbstractVector{Int})
    """
    Extract integrated Berry curvature contributions for each (band, ky) pair.
    Averages the Berry curvature over all kx for each band and ky.
    
    berry_curvature: shape (nbands, Nkx, Nky)
    ky_indices: which ky indices to extract (e.g., 1:Nky for evenly spaced, or specific indices)
    
    Returns: shape (nbands, length(ky_indices))
    """
    nbands, Nkx, Nky = size(berry_curvature)
    result = zeros(Float64, nbands, length(ky_indices))
    
    for (i, iy) in enumerate(ky_indices)
        for band in 1:nbands
            result[band, i] = mean(@view berry_curvature[band, :, iy])
        end
    end
    
    return result
end

function map_cumulative_chern_to_ribbon_energies(
    ribbon_energies::AbstractMatrix{<:Real},
    accumulation_energies::AbstractVector{<:Real},
    cumulative_chern::AbstractVector{<:Real}
)
    length(accumulation_energies) == length(cumulative_chern) ||
        error("accumulation_energies and cumulative_chern must have the same length.")
    isempty(accumulation_energies) && error("accumulation_energies cannot be empty.")

    mapped = zeros(Float64, size(ribbon_energies))

    @inbounds for idx in eachindex(ribbon_energies)
        e = ribbon_energies[idx]
        pos = searchsortedlast(accumulation_energies, e)

        if pos <= 0
            mapped[idx] = cumulative_chern[1]
        elseif pos >= length(cumulative_chern)
            mapped[idx] = cumulative_chern[end]
        else
            mapped[idx] = cumulative_chern[pos]
        end
    end

    return mapped
end

function map_chern_contribution_density_to_ribbon_energies(
    ribbon_energies::AbstractMatrix{<:Real},
    contribution_energies::AbstractVector{<:Real},
    contribution_weights::AbstractVector{<:Real};
    nbins::Int=600,
)
    length(contribution_energies) == length(contribution_weights) ||
        error("contribution_energies and contribution_weights must have the same length.")
    isempty(contribution_energies) && error("contribution_energies cannot be empty.")
    nbins > 1 || error("nbins must be greater than 1.")

    emin = min(minimum(ribbon_energies), minimum(contribution_energies))
    emax = max(maximum(ribbon_energies), maximum(contribution_energies))

    if !(emax > emin)
        return zeros(Float64, size(ribbon_energies))
    end

    edges = collect(range(emin, emax; length=nbins + 1))
    bin_width = edges[2] - edges[1]

    bin_weights = zeros(Float64, nbins)
    @inbounds for (e, w) in zip(contribution_energies, contribution_weights)
        idx = clamp(searchsortedlast(edges, e), 1, nbins)
        bin_weights[idx] += w
    end

    # Approximate dC/dE by binning the Berry-curvature weights over energy.
    contribution_density = bin_weights ./ bin_width

    mapped = zeros(Float64, size(ribbon_energies))
    @inbounds for idx in eachindex(ribbon_energies)
        bin_idx = clamp(searchsortedlast(edges, ribbon_energies[idx]), 1, nbins)
        mapped[idx] = contribution_density[bin_idx]
    end

    return mapped
end

function compute_chern_density_vs_ky_and_energy(
    bulk_ky_vals::AbstractVector{<:Real},
    plaquette_energies::Array{Float64, 3},
    berry_curvature::Array{Float64, 3};
    nbins_energy::Int=600,
    energy_range::Union{Nothing, Tuple{Float64, Float64}}=nothing,
)
    size(plaquette_energies) == size(berry_curvature) ||
        error("plaquette_energies and berry_curvature must have the same shape.")
    nbands, Nkx, Nky = size(plaquette_energies)
    length(bulk_ky_vals) == Nky ||
        error("bulk_ky_vals must have length Nky = $(Nky).")
    nbands > 0 || error("Need at least one band.")
    Nkx > 0 || error("Need at least one kx point.")
    Nky > 0 || error("Need at least one ky point.")
    nbins_energy > 1 || error("nbins_energy must be greater than 1.")

    if isnothing(energy_range)
        emin = minimum(plaquette_energies)
        emax = maximum(plaquette_energies)
    else
        emin, emax = energy_range
    end

    if !(emax > emin)
        error("Invalid energy range: emax must be greater than emin.")
    end

    energy_edges = collect(range(emin, emax; length=nbins_energy + 1))
    energy_centers = 0.5 .* (energy_edges[1:end-1] .+ energy_edges[2:end])
    bin_width = energy_edges[2] - energy_edges[1]

    # dC/dE(ky, E): rows are energy bins, columns are ky points.
    density = zeros(Float64, nbins_energy, Nky)

    @inbounds for iy in 1:Nky
        for ix in 1:Nkx, band in 1:nbands
            e = plaquette_energies[band, ix, iy]
            w = berry_curvature[band, ix, iy] / (2π)
            bin_idx = clamp(searchsortedlast(energy_edges, e), 1, nbins_energy)
            density[bin_idx, iy] += w
        end
    end

    density ./= bin_width

    return (
        ky_vals=collect(bulk_ky_vals),
        energy_centers=energy_centers,
        energy_edges=energy_edges,
        density=density,
    )
end

function compute_ribbon_spectrum_and_ipr(
    Lx::Int;
    N_ky::Int=101,
    edge_sites::Int=max(2, cld(Lx, 10)),
    berry_k_grid::Int=max(41, N_ky),
    kwargs...
)
    ky_vals = range(-π, π, length=N_ky)
    nbands = 2 * Lx
    energies = zeros(Float64, nbands, N_ky)
    iprs = zeros(Float64, nbands, N_ky)
    edge_weights = zeros(Float64, nbands, N_ky)

    for (i, ky) in enumerate(ky_vals)
        H_ribbon = ribbon_perturbed_hamiltonian_qwz(ky, Lx; kwargs...)
        spectrum = eigen(H_ribbon)

        energies[:, i] .= spectrum.values

        for band in 1:nbands
            state = @view spectrum.vectors[:, band]
            iprs[band, i] = compute_state_ipr(state)
            edge_weights[band, i] = compute_state_edge_weight(state, edge_sites)
        end
    end

    bulk_berry_data = compute_bulk_band_berry_data(
        ;
        Nkx=berry_k_grid,
        Nky=berry_k_grid,
        kwargs...,
    )

    # Energy-local contribution to Chern winding projected onto ribbon eigen-energies.
    # This approximates the density dC/dE rather than the cumulative C(E).
    chern_contribution_density_on_ribbon = map_chern_contribution_density_to_ribbon_energies(
        energies,
        bulk_berry_data.accumulation_energies,
        bulk_berry_data.accumulation_weights,
    )

    return (
        ky_vals=ky_vals,
        energies=energies,
        iprs=iprs,
        edge_weights=edge_weights,
        chern_contribution_density_on_ribbon=chern_contribution_density_on_ribbon,
        bulk_kx_vals=bulk_berry_data.kx_vals,
        bulk_ky_vals=bulk_berry_data.ky_vals,
        bulk_berry_curvature=bulk_berry_data.berry_curvature,
        bulk_plaquette_energies=bulk_berry_data.plaquette_energies,
        bulk_band_chern_numbers=bulk_berry_data.band_chern_numbers,
        accumulation_energies=bulk_berry_data.accumulation_energies,
        accumulation_weights=bulk_berry_data.accumulation_weights,
        cumulative_chern=bulk_berry_data.cumulative_chern,
    )
end

function plot_ribbon_spectrum_with_localization(
    ky_vals,
    energies,
    localization_data;
    xlabel::AbstractString="k_y",
    ylabel::AbstractString="Energy",
    colorbar_title::AbstractString="IPR",
    colour::Union{Symbol, PlotUtils.ContinuousColorGradient}=:viridis,
    title::AbstractString="",
    clims::Union{Tuple{Float64, Float64}, Nothing}=nothing,
    alpha::Float64=1.0,
    gap_closure_ky::Union{Vector{<:Real}, Nothing}=nothing,
    edge_state_extent::Union{Nothing, NamedTuple}=nothing,
    edge_in_gap_mask::Union{Nothing, AbstractVector{Bool}}=nothing,
    max_ipr::Union{Nothing, Float64}=nothing
)
    nbands = size(energies, 1)

    plot_kwargs = Dict(
        :xlabel => xlabel,
        :ylabel => ylabel,
        :legend => false,
        :colorbar => true,
        :colorbar_title => colorbar_title,
        :title => title,
        :cmap => colour,
        :alpha => alpha,
    )
    if !isnothing(clims)
        plot_kwargs[:clims] = clims
    end

    plt = plot(; plot_kwargs...)

    # Format x-axis as fractions of π
    ky_min, ky_max = extrema(ky_vals)
    xticks_pos = [-π, -π/2, -π/4, 0, π/4, π/2, π]
    xticks_labels = ["-π", "-π/2", "-π/4", "0", "π/4", "π/2", "π"]
    xticks!(plt, xticks_pos, xticks_labels)
    xlims!(plt, ky_min, ky_max)

    # Lightly shade k_y regions outside the in-gap edge-state extent.
    if !isnothing(edge_in_gap_mask) && length(edge_in_gap_mask) == length(ky_vals)
        mask = edge_in_gap_mask
        N = length(mask)
        i = 1
        while i <= N
            if !mask[i]
                j = i
                while j < N && !mask[j + 1]
                    j += 1
                end

                k_left = i == 1 ? ky_vals[1] : 0.5 * (ky_vals[i - 1] + ky_vals[i])
                k_right = j == N ? ky_vals[end] : 0.5 * (ky_vals[j] + ky_vals[j + 1])
                vspan!(plt, [k_left, k_right]; color=:black, alpha=0.05, label=false)
                i = j + 1
            else
                i += 1
            end
        end
    end

    for band in 1:nbands
        band_energies = @view energies[band, :]
        band_localization = @view localization_data[band, :]

        plot!(
            plt,
            ky_vals,
            band_energies,
            line_z=band_localization,
            c=colour,
            lw=1.0,
            label=false,
            alpha=alpha
        )
    end

    # Add vertical dotted lines for gap closure
    if !isnothing(gap_closure_ky) && !isempty(gap_closure_ky)
        for ky_closure in gap_closure_ky
            vline!(plt, [ky_closure]; linestyle=:dot, color=:red, lw=1.5, alpha=0.7, label=false)
        end
    end

    # Add annotation for edge state extent
    annotation_text = ""
    if !isnothing(edge_state_extent)
        annotation_text *= @sprintf "Edge extent (gap): %.2fπ\n" edge_state_extent.extent_in_pi
    end

    # Add annotation for max IPR
    if !isnothing(max_ipr)
        annotation_text *= @sprintf "Max IPR: %.3f" max_ipr
    end

    if !isempty(annotation_text)
        # Use data coordinates: position at top-left of plot
        ky_range = extrema(ky_vals)
        e_range = extrema(energies)
        ky_pos = ky_range[1] + 0.02 * (ky_range[2] - ky_range[1])
        e_pos = e_range[2] - 0.05 * (e_range[2] - e_range[1])
        annotate!(plt, ky_pos, e_pos, text(annotation_text, 8, :left, :top, :black))
    end

    return plt
end

function plot_chern_accumulation_vs_energy(
    accumulation_energies,
    cumulative_chern;
    xlabel::AbstractString="Accumulated Chern",
    ylabel::AbstractString="Energy",
    title::AbstractString="Chern accumulation vs energy",
)
    plt = plot(
        cumulative_chern,
        accumulation_energies;
        xlabel=xlabel,
        ylabel=ylabel,
        title=title,
        lw=2.5,
        color=:darkred,
        legend=false,
        xlims=(minimum(cumulative_chern), 1),
    )

    vline!(plt, [0.0]; linestyle=:dash, color=:black, lw=1.0, label=false)
    hline!(plt, [0.0]; linestyle=:dash, color=:black, lw=1.0, label=false)

    return plt
end

function plot_chern_density_vs_ky_energy(
    ky_vals,
    energy_centers,
    density;
    ribbon_ky_vals=nothing,
    ribbon_energies=nothing,
    xlabel::AbstractString="k_y",
    ylabel::AbstractString="Energy",
    colorbar_title::AbstractString="dC/dE",
    colour::Symbol=:RdBu,
    title::AbstractString="ky-resolved Chern contribution density",
    clims::Union{Tuple{Float64, Float64}, Nothing}=nothing,
)
    plot_kwargs = Dict(
        :xlabel => xlabel,
        :ylabel => ylabel,
        :title => title,
        :colorbar => true,
        :colorbar_title => colorbar_title,
        :cmap => colour,
    )
    if !isnothing(clims)
        plot_kwargs[:clims] = clims
    end

    plt = heatmap(ky_vals, energy_centers, density; plot_kwargs...)

    if !isnothing(ribbon_ky_vals) && !isnothing(ribbon_energies)
        nbands = size(ribbon_energies, 1)
        for band in 1:nbands
            plot!(
                plt,
                ribbon_ky_vals,
                @view(ribbon_energies[band, :]),
                color=:black,
                lw=0.5,
                alpha=0.55,
                label=false,
            )
        end
    end

    return plt
end

# Example usage to plot:
Lx = 50
As = [1.0]  # collect(1.0:0.5:3.0)
Bs = [1.0] # collect(1.0:0.5:3.0)
ms = [-1.0, -0.5, -1.5] # collect(0.5:0.5:1.0) #[-1.0] #[-2.0:0.5:0.0]
B_ys = [1.0] #collect(1.0:0.1:2.0)
gammas = collect(0.0:1.0:2.0) # 0.5, 1.0, 2.0, 3.0]
winding_numbers = [0] #collect(0:1:2)
perturbation_type = :symmetric
perturb_x = false
perturb_y = true
edge_sites = max(2, cld(Lx, 10))




@showprogress for (A, B, gamma, m, B_y, winding_num) in Iterators.product(As, Bs, gammas, ms, B_ys, winding_numbers)
    kwargs = (
        A=A,
        B=B,
        m=m, 
        gamma=gamma, 
        B_y=B_y,
        perturbation_type=perturbation_type,
        perturb_x=perturb_x,
        perturb_y=perturb_y,
        winding_number=winding_num
    )
    
    ribbon_data = compute_ribbon_spectrum_and_ipr(
        Lx;
        N_ky=101,
        edge_sites=edge_sites,
        kwargs...
    )

    # Compute gap closure and edge state extent
    gap_closure_data = compute_gap_closure_ky_extent(
        collect(ribbon_data.ky_vals),
        ribbon_data.energies,
        ribbon_data.edge_weights,
        edge_weight_threshold=0.6,
        gap_tol=0.05,
    )
    
    max_ipr_val = maximum(ribbon_data.iprs)

    plt1 = plot_ribbon_spectrum_with_localization(
        ribbon_data.ky_vals,
        ribbon_data.energies,
        ribbon_data.iprs;
        colorbar_title="IPR",
        title="Ribbon spectrum with site IPR, n=$(winding_num), Lx=$(Lx)",
        clims=(0.0, 1.0),
        gap_closure_ky=gap_closure_data.ky_closure_bounds,
        edge_state_extent=gap_closure_data,
        edge_in_gap_mask=gap_closure_data.edge_in_gap_mask,
        max_ipr=max_ipr_val
    )

    plt2 = plot_ribbon_spectrum_with_localization(
        ribbon_data.ky_vals,
        ribbon_data.energies,
        ribbon_data.chern_contribution_density_on_ribbon;
        colorbar_title="Chern Contribution Density dC/dE",
        title="Ribbon spec - Chern cont density, n=$(winding_num), Lx=$(Lx)",
        colour=cgrad(:RdBu, rev=true),
        clims=(-maximum(abs.(ribbon_data.chern_contribution_density_on_ribbon)), maximum(abs.(ribbon_data.chern_contribution_density_on_ribbon))),
    )

    # plt3 = plot_ribbon_spectrum_with_localization(
    #     ribbon_data.ky_vals,
    #     ribbon_data.energies,
    #     ribbon_data.edge_weights;
    #     colorbar_title="Edge Weight",
    #     title="Ribbon spectrum with edge weight, n=$(winding_num), Lx=$(Lx)",
    # )

    plt4 = plot_chern_accumulation_vs_energy(
        ribbon_data.accumulation_energies,
        ribbon_data.cumulative_chern;
        title="Bulk Chern acc, n=$(winding_num), Lx=$(Lx), C=$(round(sum(ribbon_data.bulk_band_chern_numbers); digits=4))",
    )

    # ky_resolved_density = compute_chern_density_vs_ky_and_energy(
    #     ribbon_data.bulk_ky_vals,
    #     ribbon_data.bulk_plaquette_energies,
    #     ribbon_data.bulk_berry_curvature;
    #     nbins_energy=600,
    # )

    # ky_density_abs_max = maximum(abs.(ky_resolved_density.density))
    # plt5 = plot_chern_density_vs_ky_energy(
    #     ky_resolved_density.ky_vals,
    #     ky_resolved_density.energy_centers,
    #     ky_resolved_density.density;
    #     ribbon_ky_vals=ribbon_data.ky_vals,
    #     ribbon_energies=ribbon_data.energies,
    #     title="ky-resolved dC/dE, n=$(winding_num), Lx=$(Lx)",
    #     clims=(-ky_density_abs_max, ky_density_abs_max),
    #     colour=:plasma
    # )

    folder_name = joinpath("plots", "ribbon_spectrum", "test")
    isdir(folder_name) || mkpath(folder_name)
    savefig(plt1, joinpath(folder_name, "ribbon_spectrum_A$(A)_B$(B)_m$(m)_gamma$(gamma)_perturbation_$(perturbation_type)_B_y$(B_y)_winding_$(winding_num)_ipr.png"))
    savefig(plt2, joinpath(folder_name, "ribbon_spectrum_A$(A)_B$(B)_m$(m)_gamma$(gamma)_perturbation_$(perturbation_type)_B_y$(B_y)_winding_$(winding_num)_berry_curv.png"))
    savefig(plt4, joinpath(folder_name, "ribbon_spectrum_A$(A)_B$(B)_m$(m)_gamma$(gamma)_perturbation_$(perturbation_type)_B_y$(B_y)_winding_$(winding_num)_chern_accumulation.png"))
    # savefig(plt5, joinpath(folder_name, "ribbon_spectrum_A$(A)_B$(B)_m$(m)_gamma$(gamma)_perturbation_$(perturbation_type)_B_y$(B_y)_winding_$(winding_num)_ky_resolved_dcdE.png"))
end
