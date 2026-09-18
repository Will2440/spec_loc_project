module SpecLocSolver

using LinearAlgebra
using SparseArrays
using Statistics
using Random

export SolverCaseConfig, run_case

const SIGMA_X = ComplexF64[0 1; 1 0]
const SIGMA_Y = ComplexF64[0 -im; im 0]
const SIGMA_Z = ComplexF64[1 0; 0 -1]
const ID2 = ComplexF64[1 0; 0 1]

Base.@kwdef struct SolverCaseConfig
    A::Float64 = 1.0
    B::Float64 = 1.0
    m::Float64 = -1.0
    B_y::Float64 = 1.0
    perturbation_type::Symbol = :symmetric
    disorder_type::Symbol = :none
    winding_number::Int = 0

    Lx_ribbon::Int = 50
    Lx_obc::Int = 14
    Ly_obc::Int = 14

    gamma_vals::Vector{Float64} = [0.0, 1.0]
    W_vals::Vector{Float64} = [0.0]
    kappa_vals::Vector{Float64} = [0.2]
    E_vals::Vector{Float64} = collect(range(-5.0, 5.0; length=101))

    N_ky::Int = 101
    Nkx_bulk::Int = 81
    Nky_bulk::Int = 81
    Nkx_3d::Int = 61
    Nky_3d::Int = 61

    ldos_target_E::Float64 = 0.0
    ldos_eta::Float64 = 0.1
    dos_eta::Float64 = 0.1
    lowest_energy_count::Int = 4

    specloc_x::Int = 7
    specloc_y::Int = 7
    seed::Int = 0

    # Storage controls for expensive eigensystems
    save_hamiltonian_matrices::Bool = true
    save_full_hamiltonian_eigensystems::Bool = true
    save_full_localiser_eigensystems::Bool = false
end

site_index_qwz(x::Int, y::Int, orb::Int, Lx::Int, Ly::Int) = 2 * ((y - 1) * Lx + (x - 1)) + orb

function real_space_perturbed_qwz_blocks(; A::Real=1.0, B::Real=1.0, m::Real=0.0, gamma::Real=0.0, perturbation_type::Symbol=:none)
    onsite = ComplexF64.((m + 2.0 * B) * SIGMA_Z)
    tx = ComplexF64.(-0.5 * B * SIGMA_Z - 0.5im * A * SIGMA_X)
    ty = ComplexF64.(-0.5 * B * SIGMA_Z - 0.5im * A * SIGMA_Y)

    if perturbation_type == :symmetric
        tx += ComplexF64.(0.5 * gamma * ID2)
        ty += ComplexF64.(0.5 * gamma * ID2)
    elseif perturbation_type == :tilt
        tx += ComplexF64.(0.5im * gamma * ID2)
    elseif perturbation_type != :none
        error("Unknown perturbation type: $perturbation_type")
    end

    return onsite, tx, ty
end

function real_space_perturbed_disordered_hamiltonian_qwz(
    Lx::Int,
    Ly::Int;
    A::Real=1.0,
    B::Real=1.0,
    m::Real=0.0,
    gamma::Real=0.0,
    perturbation_type::Symbol=:none,
    disorder_type::Symbol=:none,
    W::Real=0.0,
    periodic_x::Bool=false,
    periodic_y::Bool=false,
    sparse_output::Bool=true,
    rng::AbstractRNG=Random.default_rng(),
)
    onsite_base, tx, ty = real_space_perturbed_qwz_blocks(; A=A, B=B, m=m, gamma=gamma, perturbation_type=perturbation_type)

    dim = 2 * Lx * Ly
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
        onsite_xy = copy(onsite_base)
        if W > 0.0
            random_val = W * (rand(rng) - 0.5)
            if disorder_type == :anderson
                onsite_xy += random_val * ID2
            elseif disorder_type == :mass
                onsite_xy += random_val * SIGMA_Z
            elseif disorder_type != :none
                error("Unknown disorder type: $disorder_type")
            end
        end

        add_block!(H, (x, y), (x, y), onsite_xy)

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
    winding_number::Int=0,
)::Matrix{ComplexF64}
    d_x = A * sin(kx)
    d_y = A * sin(ky)
    d_z = m + 2B - B * cos(kx) - B * B_y * cos(ky)
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
        d_x, d_y = (
            cos_nkx * d_x - sin_nkx * d_y,
            cos_nkx * d_y + sin_nkx * d_x,
        )
    end

    H_k = scalar_term * ID2 + d_x * SIGMA_X + d_y * SIGMA_Y + d_z * SIGMA_Z
    return ComplexF64.(H_k)
end

function ribbon_perturbed_hamiltonian_qwz(
    ky::Real,
    Lx::Int;
    A::Real=1.0,
    B::Real=1.0,
    m::Real=0.0,
    gamma::Real=0.0,
    B_y::Real=1.0,
    perturbation_type::Symbol=:none,
    perturb_x::Bool=false,
    perturb_y::Bool=true,
    winding_number::Int=0,
)
    tx = ComplexF64.(-0.5 * B * SIGMA_Z - 0.5im * A * SIGMA_X)

    if perturb_x
        if perturbation_type == :symmetric
            tx += ComplexF64.(0.5 * gamma * ID2)
        elseif perturbation_type == :tilt
            tx += ComplexF64.(0.5im * gamma * ID2)
        elseif perturbation_type != :none
            error("Unknown perturbation type: $perturbation_type")
        end
    end

    onsite_ky = (m + 2.0 * B - B * B_y * cos(ky)) * SIGMA_Z + A * sin(ky) * SIGMA_Y

    if perturb_y
        if perturbation_type == :symmetric
            onsite_ky += gamma * cos(ky) * ID2
        elseif perturbation_type != :tilt && perturbation_type != :none
            error("Unknown perturbation type: $perturbation_type")
        end
    end

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

function build_position_operators(Lx::Int, Ly::Int)
    dim = 2 * Lx * Ly
    xvec = zeros(Float64, dim)
    yvec = zeros(Float64, dim)
    for y in 1:Ly, x in 1:Lx
        base = 2 * ((y - 1) * Lx + (x - 1))
        xvec[base + 1] = x
        xvec[base + 2] = x
        yvec[base + 1] = y
        yvec[base + 2] = y
    end
    return Diagonal(ComplexF64.(xvec)), Diagonal(ComplexF64.(yvec))
end

function build_localiser_matrix(
    Hshift::AbstractMatrix{ComplexF64},
    Ablock::AbstractMatrix{ComplexF64},
    Bblock::AbstractMatrix{ComplexF64},
)
    return [Hshift  Ablock - im * Bblock;
            Ablock + im * Bblock  -Hshift]
end

function spectral_localiser_signature_gap_from_blocks(
    Hshift::AbstractMatrix{ComplexF64},
    Ablock::AbstractMatrix{ComplexF64},
    Bblock::AbstractMatrix{ComplexF64};
    zero_tol::Real=1e-12,
)
    L = build_localiser_matrix(Hshift, Ablock, Bblock)
    vals = eigvals(Hermitian(L))
    revals = real(vals)
    npos = count(>(zero_tol), revals)
    nneg = count(<(-zero_tol), revals)
    minabs = minimum(abs.(revals))
    signature = npos - nneg
    return signature, minabs, vals
end

function spectral_localiser_eigendecomp_from_blocks(
    Hshift::AbstractMatrix{ComplexF64},
    Ablock::AbstractMatrix{ComplexF64},
    Bblock::AbstractMatrix{ComplexF64};
    zero_tol::Real=1e-12,
)
    L = build_localiser_matrix(Hshift, Ablock, Bblock)
    F = eigen(Hermitian(L))
    vals = F.values
    revals = real(vals)
    npos = count(>(zero_tol), revals)
    nneg = count(<(-zero_tol), revals)
    minabs = minimum(abs.(revals))
    signature = npos - nneg
    return signature, minabs, vals, F.vectors
end

function spectral_localiser_signature_gap(
    H::AbstractMatrix{ComplexF64},
    X::Diagonal,
    Y::Diagonal,
    x0::Real,
    y0::Real,
    E::Real;
    kappa::Real=1.0,
    zero_tol::Real=1e-12,
)
    D = size(H, 1)
    I_D = Matrix{ComplexF64}(I, D, D)
    Xmat = Matrix(X)
    Ymat = Matrix(Y)

    Hshift = H - E * I_D
    Ablock = kappa * (Xmat - x0 * I_D)
    Bblock = kappa * (Ymat - y0 * I_D)

    return spectral_localiser_signature_gap_from_blocks(Hshift, Ablock, Bblock; zero_tol=zero_tol)
end

function compute_state_ipr(state::AbstractVector{<:Number})
    probabilities = abs2.(state)
    site_probabilities = @views probabilities[1:2:end] .+ probabilities[2:2:end]
    return sum(site_probabilities .^ 2)
end

function compute_dos(eigs::Vector{Float64}, E_range::AbstractVector{<:Real}; eta::Float64=0.1)
    dos = zeros(Float64, length(E_range))
    gaussian(E, E0, sigma) = (1 / (sigma * sqrt(2π))) * exp(-0.5 * ((E - E0) / sigma)^2)
    for (i, E_target) in enumerate(E_range)
        dos[i] = sum(gaussian.(eigs, E_target, eta))
    end
    return dos
end

function compute_ldos(F::Eigen, Lx::Int, Ly::Int, target_E::Real; eta::Float64=0.1)
    ldos = zeros(Float64, Lx, Ly)
    eigs = real(F.values)
    gaussian(E, E0, sigma) = (1 / (sigma * sqrt(2π))) * exp(-0.5 * ((E - E0) / sigma)^2)

    for i in 1:length(eigs)
        weight = gaussian(eigs[i], target_E, eta)
        if weight > 1e-6
            psi = F.vectors[:, i]
            for y in 1:Ly, x in 1:Lx
                base = 2 * ((y - 1) * Lx + (x - 1))
                prob_density = abs2(psi[base + 1]) + abs2(psi[base + 2])
                ldos[x, y] += weight * prob_density
            end
        end
    end

    return ldos
end

function lowest_energy_ldos(F::Eigen, Lx::Int, Ly::Int; nlowest::Int=4)
    eigs = real(F.values)
    order = sortperm(abs.(eigs))
    n = min(nlowest, length(order))
    ldos = zeros(Float64, Lx, Ly)
    for idx in order[1:n]
        psi = F.vectors[:, idx]
        for y in 1:Ly, x in 1:Lx
            base = 2 * ((y - 1) * Lx + (x - 1))
            ldos[x, y] += abs2(psi[base + 1]) + abs2(psi[base + 2])
        end
    end
    return ldos
end

function compute_bulk_band_berry_data(; Nkx::Int=81, Nky::Int=81, A::Real=1.0, B::Real=1.0, m::Real=0.0, gamma::Real=0.0, B_y::Real=1.0, perturbation_type::Symbol=:none, winding_number::Int=0)
    kx_vals = collect(range(-π, π; length=Nkx + 1))[1:end-1]
    ky_vals = collect(range(-π, π; length=Nky + 1))[1:end-1]
    nbands = 2

    energies = zeros(Float64, nbands, Nkx, Nky)
    eigenvectors = Array{ComplexF64}(undef, 2, nbands, Nkx, Nky)

    for (ix, kx) in enumerate(kx_vals), (iy, ky) in enumerate(ky_vals)
        spectrum = eigen(Hermitian(k_space_perturbed_qwz_hamiltonian(kx, ky; A=A, B=B, m=m, gamma=gamma, B_y=B_y, perturbation_type=perturbation_type, winding_number=winding_number)))
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

    flat_energies = vec(plaquette_energies)
    flat_weights = vec(berry_curvature) ./ (2π)
    order = sortperm(flat_energies)

    return (
        kx_vals=kx_vals,
        ky_vals=ky_vals,
        berry_curvature=berry_curvature,
        plaquette_energies=plaquette_energies,
        accumulation_energies=flat_energies[order],
        accumulation_weights=flat_weights[order],
        cumulative_chern=cumsum(flat_weights[order]),
    )
end

function map_chern_contribution_density_to_ribbon_energies(
    ribbon_energies::AbstractMatrix{<:Real},
    contribution_energies::AbstractVector{<:Real},
    contribution_weights::AbstractVector{<:Real};
    nbins::Int=600,
)
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

    contribution_density = bin_weights ./ bin_width
    mapped = zeros(Float64, size(ribbon_energies))

    @inbounds for idx in eachindex(ribbon_energies)
        bin_idx = clamp(searchsortedlast(edges, ribbon_energies[idx]), 1, nbins)
        mapped[idx] = contribution_density[bin_idx]
    end

    return mapped
end

function compute_ribbon_spectrum_bundle(
    Lx::Int;
    N_ky::Int,
    A::Real,
    B::Real,
    m::Real,
    gamma::Real,
    B_y::Real,
    perturbation_type::Symbol,
    winding_number::Int,
)
    ky_vals = range(-π, π, length=N_ky)
    nbands = 2 * Lx
    energies = zeros(Float64, nbands, N_ky)
    iprs = zeros(Float64, nbands, N_ky)

    for (i, ky) in enumerate(ky_vals)
        H_ribbon = ribbon_perturbed_hamiltonian_qwz(ky, Lx; A=A, B=B, m=m, gamma=gamma, B_y=B_y, perturbation_type=perturbation_type, winding_number=winding_number)
        spectrum = eigen(H_ribbon)
        energies[:, i] .= spectrum.values
        for band in 1:nbands
            iprs[band, i] = compute_state_ipr(@view spectrum.vectors[:, band])
        end
    end

    bulk_berry_data = compute_bulk_band_berry_data(
        Nkx=max(41, N_ky),
        Nky=max(41, N_ky),
        A=A,
        B=B,
        m=m,
        gamma=gamma,
        B_y=B_y,
        perturbation_type=perturbation_type,
        winding_number=winding_number,
    )

    chern_density = map_chern_contribution_density_to_ribbon_energies(
        energies,
        bulk_berry_data.accumulation_energies,
        bulk_berry_data.accumulation_weights,
    )

    return (
        ky_vals=collect(ky_vals),
        energies=energies,
        iprs=iprs,
        chern_contribution_density_on_ribbon=chern_density,
        accumulation_energies=bulk_berry_data.accumulation_energies,
        cumulative_chern=bulk_berry_data.cumulative_chern,
    )
end

function compute_3d_bandstructure_data(; A::Real, B::Real, m::Real, gamma::Real, B_y::Real, perturbation_type::Symbol, winding_number::Int, Nkx::Int, Nky::Int)
    kx_vals = collect(range(-π, π; length=Nkx + 1))[1:end-1]
    ky_vals = collect(range(-π, π; length=Nky + 1))[1:end-1]
    e_minus = zeros(Float64, length(kx_vals), length(ky_vals))
    e_plus = zeros(Float64, length(kx_vals), length(ky_vals))

    for (ix, kx) in enumerate(kx_vals), (iy, ky) in enumerate(ky_vals)
        vals = eigvals(Hermitian(k_space_perturbed_qwz_hamiltonian(kx, ky; A=A, B=B, m=m, gamma=gamma, B_y=B_y, perturbation_type=perturbation_type, winding_number=winding_number)))
        e_minus[ix, iy] = vals[1]
        e_plus[ix, iy] = vals[2]
    end

    return (kx_vals=kx_vals, ky_vals=ky_vals, e_minus=e_minus, e_plus=e_plus)
end

function run_case(cfg::SolverCaseConfig)
    rng = cfg.seed == 0 ? Random.default_rng() : MersenneTwister(cfg.seed)

    gammas = cfg.gamma_vals
    Ws = cfg.W_vals
    kappas = cfg.kappa_vals
    Es = cfg.E_vals

    ng, nW, nk, nE = length(gammas), length(Ws), length(kappas), length(Es)

    specloc_gap = zeros(Float64, ng, nW, nk, nE)
    specloc_signature = zeros(Int, ng, nW, nk, nE)
    specloc_index = zeros(Float64, ng, nW, nk, nE)

    spectrum_vs_gamma = Vector{Vector{Float64}}(undef, ng)
    spectrum_vs_W = Vector{Vector{Float64}}(undef, nW)
    spectrum_vs_kappa = Vector{Vector{Float64}}(undef, nk)

    # Reference indices for spectra cuts
    iW0 = cld(nW, 2)
    ik0 = cld(nk, 2)
    iE0 = cld(nE, 2)
    ig0 = cld(ng, 2)

    x0 = cfg.specloc_x
    y0 = cfg.specloc_y
    X, Y = build_position_operators(cfg.Lx_obc, cfg.Ly_obc)

    # Cache H and eigensystems per (gamma, W)
    H_cache = Matrix{ComplexF64}[]
    eig_cache = Eigen[]
    gamma_w_pairs = Tuple{Float64, Float64}[]

    for gamma in gammas, W in Ws
        H = real_space_perturbed_disordered_hamiltonian_qwz(
            cfg.Lx_obc,
            cfg.Ly_obc;
            A=cfg.A,
            B=cfg.B,
            m=cfg.m,
            gamma=gamma,
            perturbation_type=cfg.perturbation_type,
            disorder_type=cfg.disorder_type,
            W=W,
            sparse_output=false,
            rng=rng,
        )
        push!(H_cache, H)
        push!(eig_cache, eigen(Hermitian(H)))
        push!(gamma_w_pairs, (gamma, W))
    end

    pair_index(gi, wi) = (gi - 1) * nW + wi

    D = 2 * cfg.Lx_obc * cfg.Ly_obc
    I_D = Matrix{ComplexF64}(I, D, D)
    Xmat = Matrix(X)
    Ymat = Matrix(Y)
    Xshift_base = Xmat - x0 * I_D
    Yshift_base = Ymat - y0 * I_D
    Ablocks = [kappa * Xshift_base for kappa in kappas]
    Bblocks = [kappa * Yshift_base for kappa in kappas]

    localiser_eigvals = cfg.save_full_localiser_eigensystems ? Array{Vector{Float64}}(undef, ng, nW, nk, nE) : nothing
    localiser_eigvecs = cfg.save_full_localiser_eigensystems ? Array{Matrix{ComplexF64}}(undef, ng, nW, nk, nE) : nothing

    for gi in 1:ng, wi in 1:nW
        H = H_cache[pair_index(gi, wi)]
        for ei in 1:nE
            Hshift = H - Es[ei] * I_D
            for ki in 1:nk
                Ablock = Ablocks[ki]
                Bblock = Bblocks[ki]

                if cfg.save_full_localiser_eigensystems
                    sig, gap, vals, vecs = spectral_localiser_eigendecomp_from_blocks(Hshift, Ablock, Bblock)
                    localiser_eigvals[gi, wi, ki, ei] = collect(real(vals))
                    localiser_eigvecs[gi, wi, ki, ei] = vecs
                else
                    sig, gap, vals = spectral_localiser_signature_gap_from_blocks(Hshift, Ablock, Bblock)
                end

                specloc_signature[gi, wi, ki, ei] = sig
                specloc_index[gi, wi, ki, ei] = 0.5 * sig
                specloc_gap[gi, wi, ki, ei] = gap

                # Reuse spectra from the main scan for reference cuts to avoid
                # recomputing the same expensive localiser eigensolves.
                if wi == iW0 && ki == ik0 && ei == iE0
                    spectrum_vs_gamma[gi] = sort(real(vals))
                end
                if gi == ig0 && ki == ik0 && ei == iE0
                    spectrum_vs_W[wi] = sort(real(vals))
                end
                if gi == ig0 && wi == iW0 && ei == iE0
                    spectrum_vs_kappa[ki] = sort(real(vals))
                end
            end
        end
    end

    # OBC eigensystem at reference disorder strength for DOS / LDOS.
    # Also retain gamma-resolved slices at fixed W index so viewer gamma changes are meaningful.
    F_obc = eig_cache[pair_index(ig0, iW0)]
    obc_eigs = real(F_obc.values)
    dos_vals = compute_dos(obc_eigs, Es; eta=cfg.dos_eta)
    ldos_target = compute_ldos(F_obc, cfg.Lx_obc, cfg.Ly_obc, cfg.ldos_target_E; eta=cfg.ldos_eta)
    ldos_lowest = lowest_energy_ldos(F_obc, cfg.Lx_obc, cfg.Ly_obc; nlowest=cfg.lowest_energy_count)

    dos_by_gamma = Dict{Float64, Vector{Float64}}()
    ldos_target_by_gamma = Dict{Float64, Matrix{Float64}}()
    ldos_lowest_by_gamma = Dict{Float64, Matrix{Float64}}()
    for gi in 1:ng
        gamma = gammas[gi]
        Fg = eig_cache[pair_index(gi, iW0)]
        eigs_g = real(Fg.values)
        dos_by_gamma[gamma] = compute_dos(eigs_g, Es; eta=cfg.dos_eta)
        ldos_target_by_gamma[gamma] = compute_ldos(Fg, cfg.Lx_obc, cfg.Ly_obc, cfg.ldos_target_E; eta=cfg.ldos_eta)
        ldos_lowest_by_gamma[gamma] = lowest_energy_ldos(Fg, cfg.Lx_obc, cfg.Ly_obc; nlowest=cfg.lowest_energy_count)
    end

    # Ribbon and 3D band at each gamma, but only for W=0 slice semantics.
    ribbon_by_gamma = Dict{Float64, Any}()
    band3d_by_gamma = Dict{Float64, Any}()
    for gamma in gammas
        ribbon_by_gamma[gamma] = compute_ribbon_spectrum_bundle(
            cfg.Lx_ribbon;
            N_ky=cfg.N_ky,
            A=cfg.A,
            B=cfg.B,
            m=cfg.m,
            gamma=gamma,
            B_y=cfg.B_y,
            perturbation_type=cfg.perturbation_type,
            winding_number=cfg.winding_number,
        )

        band3d_by_gamma[gamma] = compute_3d_bandstructure_data(
            A=cfg.A,
            B=cfg.B,
            m=cfg.m,
            gamma=gamma,
            B_y=cfg.B_y,
            perturbation_type=cfg.perturbation_type,
            winding_number=cfg.winding_number,
            Nkx=cfg.Nkx_3d,
            Nky=cfg.Nky_3d,
        )
    end

    hamiltonian_saved = Dict(
        "gamma_w_pairs" => gamma_w_pairs,
        "matrices" => cfg.save_hamiltonian_matrices ? H_cache : nothing,
        "eigenvalues" => cfg.save_full_hamiltonian_eigensystems ? [collect(real(F.values)) for F in eig_cache] : nothing,
        "eigenvectors" => cfg.save_full_hamiltonian_eigensystems ? [F.vectors for F in eig_cache] : nothing,
    )

    localiser_saved = Dict(
        "full_saved" => cfg.save_full_localiser_eigensystems,
        "eigenvalues" => cfg.save_full_localiser_eigensystems ? localiser_eigvals : nothing,
        "eigenvectors" => cfg.save_full_localiser_eigensystems ? localiser_eigvecs : nothing,
    )

    return Dict(
        "metadata" => Dict(
            "A" => cfg.A,
            "B" => cfg.B,
            "m" => cfg.m,
            "B_y" => cfg.B_y,
            "perturbation_type" => String(cfg.perturbation_type),
            "disorder_type" => String(cfg.disorder_type),
            "winding_number" => cfg.winding_number,
            "Lx_ribbon" => cfg.Lx_ribbon,
            "Lx_obc" => cfg.Lx_obc,
            "Ly_obc" => cfg.Ly_obc,
            "specloc_x" => x0,
            "specloc_y" => y0,
        ),
        "axes" => Dict(
            "gammas" => gammas,
            "Ws" => Ws,
            "kappas" => kappas,
            "Es" => Es,
        ),
        "specloc" => Dict(
            "gap" => specloc_gap,
            "signature" => specloc_signature,
            "index" => specloc_index,
            "spectrum_vs_gamma" => spectrum_vs_gamma,
            "spectrum_vs_W" => spectrum_vs_W,
            "spectrum_vs_kappa" => spectrum_vs_kappa,
            "spectrum_ref_indices" => Dict("gamma" => ig0, "W" => iW0, "kappa" => ik0, "E" => iE0),
            "x0" => x0,
            "y0" => y0,
        ),
        "obc" => Dict(
            "eigenvalues" => obc_eigs,
            "dos" => dos_vals,
            "ldos_target" => ldos_target,
            "ldos_lowest" => ldos_lowest,
            "ldos_target_E" => cfg.ldos_target_E,
            "dos_by_gamma" => dos_by_gamma,
            "ldos_target_by_gamma" => ldos_target_by_gamma,
            "ldos_lowest_by_gamma" => ldos_lowest_by_gamma,
            "obc_ref_indices" => Dict("gamma" => ig0, "W" => iW0),
        ),
        "hamiltonian_eigensystems" => hamiltonian_saved,
        "localiser_eigensystems" => localiser_saved,
        "ribbon_by_gamma" => ribbon_by_gamma,
        "band3d_by_gamma" => band3d_by_gamma,
    )
end

end # module
