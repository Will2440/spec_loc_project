module SpecLocSolver

# Dependencies (stdlib): LinearAlgebra, SparseArrays, Random
#          (registered): KrylovKit  — install with: ]add KrylovKit
using LinearAlgebra
using SparseArrays
using Random
using KrylovKit

export SolverCaseConfig, run_case

const SIGMA_X = ComplexF64[0 1; 1 0]
const SIGMA_Y = ComplexF64[0 -im; im 0]
const SIGMA_Z = ComplexF64[1 0; 0 -1]
const ID2     = ComplexF64[1 0; 0 1]

# =====================================================================
# Config
# =====================================================================

Base.@kwdef struct SolverCaseConfig
    A::Float64 = 1.0
    B::Float64 = 1.0
    m::Float64 = -1.0
    # Available perturbation types:
    #   :none
    #   :sym_cos_sum   — +γ(cos kx + cos ky)        [legacy alias: :symmetric]
    #   :sym_cos_diff  — +γ(cos kx − cos ky)
    #   :sym_cos_add   — +γ cos(kx+ky)  (diagonal hoppings)
    #   :sym_cos_sub   — +γ cos(kx−ky)  (diagonal hoppings)
    #   :asym_sin_sum  — +γ(sin kx + sin ky)
    #   :asym_sin_diff — +γ(sin kx − sin ky)
    #   :asym_sin_add  — +γ sin(kx+ky)  (diagonal hoppings)
    #   :asym_sin_sub  — +γ sin(kx−ky)  (diagonal hoppings)
    #   :tilt          — legacy −γ sin kx
    perturbation_type::Symbol = :sym_cos_sum
    disorder_type::Symbol     = :none   # :none, :anderson, :mass

    Lx_obc::Int = 14
    Ly_obc::Int = 14

    gamma_vals::Vector{Float64} = [0.0, 1.0]
    W_vals::Vector{Float64}     = [0.0]
    kappa_vals::Vector{Float64} = [0.2]
    E_vals::Vector{Float64}     = collect(range(-5.0, 5.0; length=101))

    specloc_x::Int = 7
    specloc_y::Int = 7
    seed::Int      = 0

    # Orbital sublattice embedding for position operators
    orbital_displacement::Float64 = 0.0
    phi::Float64                   = 0.0

    # Disorder averaging (ignored when disorder_type==:none or W==0)
    n_disorder_realisations::Int = 1

    # Kappa scaling: if true, effective κ = kappa_scales[i] × Lx_obc
    scale_kappa_to_L::Bool        = false
    kappa_scales::Vector{Float64} = [0.0004]

    compute_specloc_gap::Bool  = true  # use KrylovKit shift-invert for gap
    n_gap_eigenpairs::Int      = 4
    n_spectrum_eigenpairs::Int = 20    # stored for spectrum-vs-* slices
end

# =====================================================================
# Hopping blocks
# =====================================================================
# Convention:
#   H(k) = onsite + tx e^{ikx} + tx' e^{-ikx} + ty e^{iky} + ty' e^{-iky}
#         + txy e^{i(kx+ky)} + txy' e^{-i(kx+ky)}
#         + txmy e^{i(kx-ky)} + txmy' e^{-i(kx-ky)}

function real_space_perturbed_qwz_blocks(;
    A::Real=1.0, B::Real=1.0, m::Real=0.0, gamma::Real=0.0,
    perturbation_type::Symbol=:none,
)
    onsite = ComplexF64.((m + 2.0*B) * SIGMA_Z)
    tx     = ComplexF64.(-0.5*B * SIGMA_Z - 0.5im*A * SIGMA_X)
    ty     = ComplexF64.(-0.5*B * SIGMA_Z - 0.5im*A * SIGMA_Y)
    txy    = zeros(ComplexF64, 2, 2)
    txmy   = zeros(ComplexF64, 2, 2)

    if perturbation_type == :sym_cos_sum || perturbation_type == :symmetric
        tx   .+= 0.5*gamma * ID2;   ty .+= 0.5*gamma * ID2
    elseif perturbation_type == :sym_cos_diff
        tx   .+= 0.5*gamma * ID2;   ty .-= 0.5*gamma * ID2
    elseif perturbation_type == :sym_cos_add
        txy  .+= 0.5*gamma * ID2
    elseif perturbation_type == :sym_cos_sub
        txmy .+= 0.5*gamma * ID2
    elseif perturbation_type == :asym_sin_sum
        tx   .+= -0.5im*gamma * ID2;  ty .+= -0.5im*gamma * ID2
    elseif perturbation_type == :asym_sin_diff
        tx   .+= -0.5im*gamma * ID2;  ty .+=  0.5im*gamma * ID2
    elseif perturbation_type == :asym_sin_add
        txy  .+= -0.5im*gamma * ID2
    elseif perturbation_type == :asym_sin_sub
        txmy .+= -0.5im*gamma * ID2
    elseif perturbation_type == :tilt
        tx   .+=  0.5im*gamma * ID2   # legacy: scalar −γ sin kx
    elseif perturbation_type != :none
        error("Unknown perturbation type: $perturbation_type")
    end

    return onsite, tx, ty, txy, txmy
end

# =====================================================================
# Sparse OBC Hamiltonian (COO assembly)
# =====================================================================

site_index(x, y, orb, Lx, Ly) = 2*((y-1)*Lx + (x-1)) + orb

function build_hamiltonian(
    Lx::Int, Ly::Int;
    A::Real=1.0, B::Real=1.0, m::Real=0.0, gamma::Real=0.0,
    perturbation_type::Symbol=:none,
    disorder_type::Symbol=:none, W::Real=0.0,
    rng::AbstractRNG=Random.default_rng(),
)::SparseMatrixCSC{ComplexF64,Int}
    onsite_base, tx, ty, txy, txmy =
        real_space_perturbed_qwz_blocks(; A=A, B=B, m=m, gamma=gamma,
                                          perturbation_type=perturbation_type)
    has_txy  = any(x -> abs(x) > 1e-15, txy)
    has_txmy = any(x -> abs(x) > 1e-15, txmy)
    dim = 2*Lx*Ly

    Is = Int[];        sizehint!(Is, 12*dim)
    Js = Int[];        sizehint!(Js, 12*dim)
    Vs = ComplexF64[]; sizehint!(Vs, 12*dim)

    @inline function push_block!(rx, ry, cx, cy, blk)
        rb = site_index(rx, ry, 1, Lx, Ly)
        cb = site_index(cx, cy, 1, Lx, Ly)
        @inbounds for a in 0:1, b in 0:1
            push!(Is, rb+a); push!(Js, cb+b)
            push!(Vs, ComplexF64(blk[a+1, b+1]))
        end
    end

    for y in 1:Ly, x in 1:Lx
        os = copy(onsite_base)
        if W > 0.0
            v = W*(rand(rng) - 0.5)
            if     disorder_type == :anderson; os .+= v * ID2
            elseif disorder_type == :mass;     os .+= v * SIGMA_Z
            elseif disorder_type != :none;     error("Unknown disorder type: $disorder_type")
            end
        end
        push_block!(x, y, x, y, os)

        if x < Lx
            push_block!(x+1, y, x, y, tx);   push_block!(x, y, x+1, y, tx')
        end
        if y < Ly
            push_block!(x, y+1, x, y, ty);   push_block!(x, y, x, y+1, ty')
        end
        if has_txy && x < Lx && y < Ly
            push_block!(x+1, y+1, x, y, txy);   push_block!(x, y, x+1, y+1, txy')
        end
        if has_txmy && x < Lx && y > 1
            push_block!(x+1, y-1, x, y, txmy);   push_block!(x, y, x+1, y-1, txmy')
        end
    end

    return sparse(Is, Js, Vs, dim, dim)
end

# =====================================================================
# Position operators with optional orbital sublattice displacement
# =====================================================================

function build_position_operators(
    Lx::Int, Ly::Int;
    orbital_displacement::Real=0.0, phi::Real=0.0,
)
    dim  = 2*Lx*Ly
    xvec = zeros(Float64, dim)
    yvec = zeros(Float64, dim)
    ddx  = orbital_displacement * cos(phi)
    ddy  = orbital_displacement * sin(phi)
    for y in 1:Ly, x in 1:Lx
        base = 2*((y-1)*Lx + (x-1))
        xvec[base+1] = x + ddx/2;   xvec[base+2] = x - ddx/2
        yvec[base+1] = y + ddy/2;   yvec[base+2] = y - ddy/2
    end
    return Diagonal(ComplexF64.(xvec)), Diagonal(ComplexF64.(yvec))
end

# =====================================================================
# Spectral localiser (sparse block matrix)
# =====================================================================

function build_localiser(
    Hshift::SparseMatrixCSC{ComplexF64,Int},
    A_kap::SparseMatrixCSC{ComplexF64,Int},
    B_kap::SparseMatrixCSC{ComplexF64,Int},
)::SparseMatrixCSC{ComplexF64,Int}
    return [Hshift           A_kap - im*B_kap;
            A_kap + im*B_kap  -Hshift]
end

# =====================================================================
# LDLt signature  (Sylvester's law — no full eigendecomposition needed)
# =====================================================================

function signature_ldlt(
    Lherm::Hermitian{ComplexF64,SparseMatrixCSC{ComplexF64,Int}};
    max_tries::Int=6, base_shift::Float64=1e-12,
)::Int
    shift = 0.0
    for _ in 0:max_tries
        try
            F    = ldlt(Lherm; shift=shift)
            Dv   = real.(diag(sparse(F.LD)))
            return count(>(0), Dv) - count(<(0), Dv)
        catch err
            es = string(typeof(err))
            (err isa LinearAlgebra.ZeroPivotException ||
             occursin("CHOLMOD", es) || occursin("Pivot", es)) ||
                rethrow(err)
            shift = shift == 0.0 ? base_shift : shift*10.0
        end
    end
    @warn "signature_ldlt: LDLt failed; falling back to dense eigvals"
    evs = real.(eigvals(Hermitian(Matrix(Lherm))))
    return count(>(0), evs) - count(<(0), evs)
end

# =====================================================================
# Shift-and-invert KrylovKit for lowest-|eigenvalue| spectrum
# =====================================================================

function low_lying_spectrum(
    L::SparseMatrixCSC{ComplexF64,Int};
    n::Int=10, tol::Real=1e-8, maxiter::Int=300, reg::Float64=1e-10,
)::Vector{Float64}
    sz   = size(L, 1)
    nreq = min(n, sz)
    for attempt in 0:1
        Leff = attempt == 0 ? L : L + reg*I(sz)
        try
            F     = lu(Leff)
            Linv  = v -> F \ v
            ev, _ = eigsolve(Linv, rand(ComplexF64, sz), nreq, :LM;
                             ishermitian=true, tol=tol, maxiter=maxiter)[1:2]
            evals = sort(real.(1.0 ./ ev), by=abs)
            attempt == 1 && (evals .-= reg)
            return evals[1:min(nreq, length(evals))]
        catch err
            attempt == 0 &&
                (err isa LAPACKException ||
                 err isa LinearAlgebra.SingularException ||
                 occursin("TaskFailed", string(err))) && continue
            attempt == 1 && break
        end
    end
    @warn "low_lying_spectrum: all sparse methods failed; using dense eigvals (sz=$sz)"
    evs = real.(eigvals(Hermitian(Matrix(L))))
    return sort(evs, by=abs)[1:min(nreq, length(evs))]
end

# =====================================================================
# Main entry point — spectral localiser only
# =====================================================================

function run_case(cfg::SolverCaseConfig)
    kappas = cfg.scale_kappa_to_L ?
        [s * Float64(cfg.Lx_obc) for s in cfg.kappa_scales] : cfg.kappa_vals

    gammas = cfg.gamma_vals;  Ws = cfg.W_vals;  Es = cfg.E_vals
    ng, nW, nk, nE = length(gammas), length(Ws), length(kappas), length(Es)
    n_real = max(1, cfg.n_disorder_realisations)

    x0, y0 = cfg.specloc_x, cfg.specloc_y
    D      = 2*cfg.Lx_obc*cfg.Ly_obc

    X, Y = build_position_operators(cfg.Lx_obc, cfg.Ly_obc;
                                     orbital_displacement=cfg.orbital_displacement,
                                     phi=cfg.phi)
    Xsh_sp = sparse(Diagonal(ComplexF64.(diag(X) .- x0)))
    Ysh_sp = sparse(Diagonal(ComplexF64.(diag(Y) .- y0)))
    Ak_sp  = [kappa * Xsh_sp for kappa in kappas]
    Bk_sp  = [kappa * Ysh_sp for kappa in kappas]

    iW0 = cld(nW, 2);  ik0 = cld(nk, 2);  iE0 = cld(nE, 2);  ig0 = cld(ng, 2)

    specloc_gap       = fill(NaN, ng, nW, nk, nE)
    specloc_signature = zeros(Int,     ng, nW, nk, nE)
    specloc_index     = zeros(Float64, ng, nW, nk, nE)

    spectrum_vs_gamma = Vector{Vector{Float64}}(undef, ng)
    spectrum_vs_W     = Vector{Vector{Float64}}(undef, nW)
    spectrum_vs_kappa = Vector{Vector{Float64}}(undef, nk)

    n_spec = cfg.n_spectrum_eigenpairs
    n_gap  = cfg.compute_specloc_gap ? cfg.n_gap_eigenpairs : 0
    I_sp   = sparse(I, D, D)

    for gi in 1:ng, wi in 1:nW
        gamma = gammas[gi];  W = Ws[wi]
        sig_sum = zeros(Float64, nk, nE)
        gap_sum = zeros(Float64, nk, nE)

        actual_real = (W == 0.0 || cfg.disorder_type == :none) ? 1 : n_real

        for r in 1:actual_real
            seed_r = cfg.seed == 0 ? 0 :
                     cfg.seed + ((gi-1)*nW + (wi-1))*max(n_real,1) + (r-1)
            rng_r  = seed_r == 0 ? Random.default_rng() : MersenneTwister(seed_r)

            H_r = build_hamiltonian(cfg.Lx_obc, cfg.Ly_obc;
                    A=cfg.A, B=cfg.B, m=cfg.m, gamma=gamma,
                    perturbation_type=cfg.perturbation_type,
                    disorder_type=cfg.disorder_type, W=W, rng=rng_r)

            for ei in 1:nE
                Hsh = H_r - ComplexF64(Es[ei]) * I_sp

                for ki in 1:nk
                    L     = build_localiser(Hsh, Ak_sp[ki], Bk_sp[ki])
                    Lherm = Hermitian(L)

                    need_γ = r==1 && wi==iW0 && ki==ik0 && ei==iE0 && !isassigned(spectrum_vs_gamma, gi)
                    need_W = r==1 && gi==ig0 && ki==ik0 && ei==iE0 && !isassigned(spectrum_vs_W, wi)
                    need_κ = r==1 && gi==ig0 && wi==iW0 && ei==iE0 && !isassigned(spectrum_vs_kappa, ki)
                    n_eig  = max(n_gap, (need_γ || need_W || need_κ) ? n_spec : 0)

                    sig = signature_ldlt(Lherm)
                    sig_sum[ki, ei] += sig

                    if n_eig > 0
                        evals = low_lying_spectrum(L; n=n_eig)
                        n_gap > 0 && (gap_sum[ki, ei] += minimum(abs.(evals)))
                        sl = evals[1:min(n_spec, length(evals))]
                        need_γ && (spectrum_vs_gamma[gi]  = sl)
                        need_W && (spectrum_vs_W[wi]      = sl)
                        need_κ && (spectrum_vs_kappa[ki]  = sl)
                    end
                end
            end
        end

        for ei in 1:nE, ki in 1:nk
            rsig = round(Int, sig_sum[ki, ei] / actual_real)
            specloc_signature[gi, wi, ki, ei] = rsig
            specloc_index[gi,     wi, ki, ei] = 0.5 * rsig
            specloc_gap[gi,       wi, ki, ei] = n_gap > 0 ? gap_sum[ki,ei] / actual_real : NaN
        end
    end

    return Dict(
        "metadata" => Dict(
            "A"=>cfg.A, "B"=>cfg.B, "m"=>cfg.m,
            "perturbation_type"       => String(cfg.perturbation_type),
            "disorder_type"           => String(cfg.disorder_type),
            "Lx_obc"=>cfg.Lx_obc, "Ly_obc"=>cfg.Ly_obc,
            "specloc_x"=>x0, "specloc_y"=>y0,
            "orbital_displacement"    => cfg.orbital_displacement,
            "phi"                     => cfg.phi,
            "n_disorder_realisations" => n_real,
            "scale_kappa_to_L"        => cfg.scale_kappa_to_L,
        ),
        "axes" => Dict("gammas"=>gammas, "Ws"=>Ws, "kappas"=>kappas, "Es"=>Es),
        "specloc" => Dict(
            "gap"=>specloc_gap, "signature"=>specloc_signature, "index"=>specloc_index,
            "spectrum_vs_gamma"   => spectrum_vs_gamma,
            "spectrum_vs_W"       => spectrum_vs_W,
            "spectrum_vs_kappa"   => spectrum_vs_kappa,
            "spectrum_ref_indices"=> Dict("gamma"=>ig0,"W"=>iW0,"kappa"=>ik0,"E"=>iE0),
            "x0"=>x0, "y0"=>y0,
        ),
    )
end

end # module SpecLocSolver
