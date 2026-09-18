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
using KrylovKit


println("Julia Threads: ", Threads.nthreads())
println("BLAS Threads:  ", LinearAlgebra.BLAS.get_num_threads())

## Pauli matrices
const sigma_x = [0 1; 1 0]
const sigma_y = [0 -im; im 0]
const sigma_z = [1 0; 0 -1]
const identity_2 = [1 0; 0 1]

function build_position_operators(Lx::Int, Ly::Int; 
    orbital_displacement::Real=0.0, 
    phi::Real=0.0
)
    nsites = Lx * Ly
    dim = 2 * nsites
    xvec = zeros(Float64, dim)
    yvec = zeros(Float64, dim)
    
    # Compute embedding displacement components
    dx_disp = orbital_displacement * cos(phi)
    dy_disp = orbital_displacement * sin(phi)
    
    for y in 1:Ly, x in 1:Lx
        base = 2 * ((y - 1) * Lx + (x - 1))
        # Orbital 1 (sublattice A): displaced by +d/2
        xvec[base + 1] = x + dx_disp / 2
        yvec[base + 1] = y + dy_disp / 2
        # Orbital 2 (sublattice B): displaced by -d/2
        xvec[base + 2] = x - dx_disp / 2
        yvec[base + 2] = y - dy_disp / 2
    end
    X = Diagonal(ComplexF64.(xvec))
    Y = Diagonal(ComplexF64.(yvec))
    return X, Y
end

# ## the OG
# function build_position_operators(
#     Lx::Int, 
#     Ly::Int;
#     orbital_displacement::Real=0.0,
#     phi::Real=0.0
# )

#     nsites = Lx * Ly
#     dim = 2 * nsites
#     xvec = zeros(Float64, dim)
#     yvec = zeros(Float64, dim)
#     for y in 1:Ly, x in 1:Lx
#         base = 2 * ((y - 1) * Lx + (x - 1))
#         xvec[base + 1] = x
#         xvec[base + 2] = x
#         yvec[base + 1] = y
#         yvec[base + 2] = y
#     end
#     X = Diagonal(ComplexF64.(xvec))
#     Y = Diagonal(ComplexF64.(yvec))
#     return X, Y
# end

function build_localiser_operator(
    H::AbstractMatrix{ComplexF64},
    X::Diagonal,
    Y::Diagonal,
    x0::Real,
    y0::Real,
    E::Real; 
    kappa::Real=1.0
)::SparseMatrixCSC{ComplexF64, Int}
    D = size(H, 1)
    
    # Keep operators fully sparse without dense matrix conversion
    Hshift = H - E * I
    Ablock = kappa * (X - x0 * I)
    Bblock = kappa * (Y - y0 * I)

    L = [Hshift  (Ablock - im * Bblock);
        (Ablock + im * Bblock)  -Hshift]

    return sparse(L)
end

function format_params_for_log(params::NamedTuple)::String
    """Helper to format parameters for logging"""
    param_strs = []
    for (k, v) in pairs(params)
        if isa(v, Float64)
            push!(param_strs, "$k=$(round(v, digits=4))")
        elseif isa(v, Int64)
            push!(param_strs, "$k=$v")
        else
            push!(param_strs, "$k=$(v)")
        end
    end
    return join(param_strs, ", ")
end

function compute_low_lying_localiser_shift_invert(
    L::SparseMatrixCSC{ComplexF64, Int};
    n_eigenpairs::Int=10,
    kk_tol::Real=1e-8,
    kk_maxiter::Int=300,
    params::NamedTuple=NamedTuple(),
    regularization::Float64=1e-10,
    max_tries::Int=2
)::Tuple{Vector{Float64}, Symbol}
    """
    Robust shift-and-invert eigenvalue solver with fallback strategies.
    
    Returns:
        (eigenvalues, status) where status ∈ {:success, :fallback_regularized, :fallback_dense}
    """
    n = size(L, 1)
    param_str = isempty(params) ? "" : "\n  Params: $(format_params_for_log(params))"
    
    # Strategy 1: Try standard LU factorization
    try
        F = lu(L)
        L_inv_action(x) = F \ x

        evals_inv, evecs, info = eigsolve(
            L_inv_action, 
            rand(ComplexF64, n), 
            n_eigenpairs, 
            :LM; 
            ishermitian=true, 
            tol=kk_tol, 
            maxiter=kk_maxiter
        )
        
        signed_evals = real.(1.0 ./ evals_inv)
        sort!(signed_evals, by=abs)
        return (signed_evals[1:min(n_eigenpairs, length(signed_evals))], :success)
        
    catch err
        if err isa LAPACKException || err isa LinearAlgebra.SingularException || 
           err isa TaskFailedException || occursin("TaskFailed", string(err))
            
            warn_msg = "Shift-invert LU factorization failed or matrix near-singular, attempting regularization..." * param_str
            @warn warn_msg
            
            # Strategy 2: Regularize by adding small shift diagonal
            L_reg = L + regularization * I(n)
            try
                F = lu(L_reg)
                L_inv_action(x) = F \ x

                evals_inv, evecs, info = eigsolve(
                    L_inv_action, 
                    rand(ComplexF64, n), 
                    n_eigenpairs, 
                    :LM; 
                    ishermitian=true, 
                    tol=kk_tol, 
                    maxiter=kk_maxiter
                )
                
                signed_evals = real.(1.0 ./ evals_inv) .- regularization
                sort!(signed_evals, by=abs)
                info_msg = "Regularized shift-invert succeeded." * param_str
                @info info_msg
                return (signed_evals[1:min(n_eigenpairs, length(signed_evals))], :fallback_regularized)
                
            catch err2
                warn_msg2 = "Regularized LU also failed, falling back to direct dense diagonalization." * param_str
                @warn warn_msg2
                
                # Strategy 3: Fall back to dense eigdecomposition
                try
                    L_dense = Matrix(L)
                    evals = eigvals(Hermitian(L_dense))
                    evals_sorted = sort(evals, by=abs)
                    @info "Dense eigendecomposition succeeded." * param_str
                    return (evals_sorted[1:min(n_eigenpairs, length(evals_sorted))], :fallback_dense)
                catch err3
                    @error "All eigenvalue methods failed!" * param_str * 
                           "\n  Error 1 (LU): $(typeof(err))\n  Error 2 (Regularized): $(typeof(err2))\n  Error 3 (Dense): $(typeof(err3))"
                    return (zeros(n_eigenpairs), :failed)
                end
            end
        else
            rethrow(err)
        end
    end
end

# ## OLD
# function real_space_perturbed_qwz_blocks(; 
#     A::Real=1.0, B::Real=1.0, m::Real=0.0, gamma::Real=0.0, perturbation_type::Symbol=:none
# )
#     onsite = ComplexF64.((m + 2.0 * B) * sigma_z)
#     tx = ComplexF64.(-0.5 * B * sigma_z - 0.5im * A * sigma_x)
#     ty = ComplexF64.(-0.5 * B * sigma_z - 0.5im * A * sigma_y)

#     # Symmetric perturbations (cosine-like, real scalar)
#     if perturbation_type == :sym_cos_sum || perturbation_type == :symmetric
#         # Symmetric addition: uniform real perturbation
#         tx += ComplexF64.(0.5 * gamma * identity_2)
#         ty += ComplexF64.(0.5 * gamma * identity_2)
#     elseif perturbation_type == :sym_cos_diff
#         # Symmetric difference: different signs in x and y
#         tx += ComplexF64.(0.5 * gamma * identity_2)
#         ty += ComplexF64.(-0.5 * gamma * identity_2)
#     elseif perturbation_type == :sym_cos_add
#         # Combined symmetric: average contribution
#         tx += ComplexF64.(0.25 * gamma * identity_2)
#         ty += ComplexF64.(0.25 * gamma * identity_2)
#     elseif perturbation_type == :sym_cos_sub
#         # Combined symmetric difference
#         tx += ComplexF64.(0.25 * gamma * identity_2)
#         ty += ComplexF64.(-0.25 * gamma * identity_2)
#     # Antisymmetric perturbations (sine-like, imaginary component)
#     elseif perturbation_type == :asym_sin_sum || perturbation_type == :antisymmetric
#         # Antisymmetric addition: imaginary perturbation
#         tx += ComplexF64.(0.5im * gamma * identity_2)
#         ty += ComplexF64.(0.5im * gamma * identity_2)
#     elseif perturbation_type == :asym_sin_diff
#         # Antisymmetric difference: opposite imaginary contributions
#         tx += ComplexF64.(0.5im * gamma * identity_2)
#         ty += ComplexF64.(-0.5im * gamma * identity_2)
#     elseif perturbation_type == :asym_sin_add
#         # Combined antisymmetric
#         tx += ComplexF64.(0.25im * gamma * identity_2)
#         ty += ComplexF64.(0.25im * gamma * identity_2)
#     elseif perturbation_type == :asym_sin_sub
#         # Combined antisymmetric difference
#         tx += ComplexF64.(0.25im * gamma * identity_2)
#         ty += ComplexF64.(-0.25im * gamma * identity_2)
#     elseif perturbation_type != :none
#         error("Unknown perturbation type: $perturbation_type")
#     end

#     return onsite, tx, ty
# end

function real_space_perturbed_qwz_blocks(; 
    A::Real=1.0, B::Real=1.0, m::Real=0.0, gamma::Real=0.0, perturbation_type::Symbol=:none
)
    # Convention:
    # H(k) = onsite 
    #      + tx*e^{i kx}   + tx'*e^{-i kx} 
    #      + ty*e^{i ky}   + ty'*e^{-i ky}
    #      + txy*e^{i(kx+ky)}  + txy'*e^{-i(kx+ky)}     (diagonal bond (1,1))
    #      + txmy*e^{i(kx-ky)} + txmy'*e^{-i(kx-ky)}    (diagonal bond (1,-1))
    onsite = ComplexF64.((m + 2.0 * B) * sigma_z)
    tx = ComplexF64.(-0.5 * B * sigma_z - 0.5im * A * sigma_x)
    ty = ComplexF64.(-0.5 * B * sigma_z - 0.5im * A * sigma_y)

    # Diagonal-neighbor hoppings — only populated by the _add/_sub perturbations.
    txy  = zeros(ComplexF64, 2, 2)   # bond (i,j) -> (i+1, j+1)
    txmy = zeros(ComplexF64, 2, 2)   # bond (i,j) -> (i+1, j-1)

    # General rule: adding c*I (c complex) to the hopping t_r at displacement r
    # contributes  c e^{ik.r} + c* e^{-ik.r} = 2Re(c) cos(k.r) - 2Im(c) sin(k.r)
    #   target +gamma*cos(k.r)  =>  c =  gamma/2      (real)
    #   target +gamma*sin(k.r)  =>  c = -i*gamma/2    (imaginary)

    if perturbation_type == :sym_cos_sum
        tx += ComplexF64.(0.5 * gamma * identity_2)
        ty += ComplexF64.(0.5 * gamma * identity_2)

    elseif perturbation_type == :sym_cos_diff
        tx += ComplexF64.(0.5 * gamma * identity_2)
        ty += ComplexF64.(-0.5 * gamma * identity_2)

    elseif perturbation_type == :sym_cos_add
        txy += ComplexF64.(0.5 * gamma * identity_2)

    elseif perturbation_type == :sym_cos_sub
        txmy += ComplexF64.(0.5 * gamma * identity_2)

    elseif perturbation_type == :asym_sin_sum
        tx += ComplexF64.(-0.5im * gamma * identity_2)
        ty += ComplexF64.(-0.5im * gamma * identity_2)

    elseif perturbation_type == :asym_sin_diff
        tx += ComplexF64.(-0.5im * gamma * identity_2)
        ty += ComplexF64.( 0.5im * gamma * identity_2)

    elseif perturbation_type == :asym_sin_add
        txy += ComplexF64.(-0.5im * gamma * identity_2)

    elseif perturbation_type == :asym_sin_sub
        txmy += ComplexF64.(-0.5im * gamma * identity_2)

    elseif perturbation_type != :none
        error("Unknown perturbation type: $perturbation_type")
    end

    return onsite, tx, ty, txy, txmy
end

site_index_qwz(x::Int, y::Int, orb::Int, Lx::Int, Ly::Int) = 2 * ((y - 1) * Lx + (x - 1)) + orb

# OPTIMIZED: Use COO format for sparse assembly
function real_space_perturbed_disordered_hamiltonian_qwz_coo(
    Lx::Int, Ly::Int;
    A::Real=1.0, B::Real=1.0, m::Real=0.0, gamma::Real=0.0,
    perturbation_type::Symbol=:none, disorder_type::Symbol=:none, W::Real=0.0,
    periodic_x::Bool=false, periodic_y::Bool=false
)
    onsite_base, tx, ty = real_space_perturbed_qwz_blocks(; A=A, B=B, m=m, gamma=gamma, perturbation_type=perturbation_type)
    
    nsites = Lx * Ly
    dim = 2 * nsites
    
    # Use COO format: collect all (row, col, value) triplets
    Is = Int[]
    Js = Int[]
    Vs = ComplexF64[]
    
    # Pre-allocate approximate capacity (overestimate is fine)
    sizehint!(Is, 20 * nsites)
    sizehint!(Js, 20 * nsites)
    sizehint!(Vs, 20 * nsites)
    
    function add_block_coo!(Is, Js, Vs, row_site::Tuple{Int, Int}, col_site::Tuple{Int, Int}, 
                            block::AbstractMatrix{<:Number}, Lx, Ly)
        (xr, yr) = row_site
        (xc, yc) = col_site
        row_base = site_index_qwz(xr, yr, 1, Lx, Ly)
        col_base = site_index_qwz(xc, yc, 1, Lx, Ly)
        for a in 0:1, b in 0:1
            push!(Is, row_base + a)
            push!(Js, col_base + b)
            push!(Vs, ComplexF64(block[a + 1, b + 1]))
        end
    end

    for y in 1:Ly, x in 1:Lx
        onsite_xy = copy(onsite_base)
        
        if W > 0.0
            random_val = W * (rand() - 0.5)
            if disorder_type == :anderson
                onsite_xy += random_val * identity_2
            elseif disorder_type == :mass
                onsite_xy += random_val * sigma_z
            elseif disorder_type != :none
                error("Unknown disorder type: $disorder_type")
            end
        end

        add_block_coo!(Is, Js, Vs, (x, y), (x, y), onsite_xy, Lx, Ly)

        if x < Lx
            add_block_coo!(Is, Js, Vs, (x + 1, y), (x, y), tx, Lx, Ly)
            add_block_coo!(Is, Js, Vs, (x, y), (x + 1, y), tx', Lx, Ly)
        elseif periodic_x
            add_block_coo!(Is, Js, Vs, (1, y), (x, y), tx, Lx, Ly)
            add_block_coo!(Is, Js, Vs, (x, y), (1, y), tx', Lx, Ly)
        end

        if y < Ly
            add_block_coo!(Is, Js, Vs, (x, y + 1), (x, y), ty, Lx, Ly)
            add_block_coo!(Is, Js, Vs, (x, y), (x, y + 1), ty', Lx, Ly)
        elseif periodic_y
            add_block_coo!(Is, Js, Vs, (x, 1), (x, y), ty, Lx, Ly)
            add_block_coo!(Is, Js, Vs, (x, y), (x, 1), ty', Lx, Ly)
        end
    end

    # Convert COO to sparse matrix (fast!)
    return sparse(Is, Js, Vs, dim, dim)
end

function fast_low_lying_localiser_spectrum(
    Lxs::AbstractVector{<:Integer}, 
    Lys::AbstractVector{<:Integer}; 
    Avals::AbstractVector{<:Real} = [1.0], 
    Bvals::AbstractVector{<:Real} = [1.0], 
    mvals::AbstractVector{<:Real} = [0.0],
    gammas::AbstractVector{<:Real} = [0.0],
    perturbation_type::Symbol = :none,
    disorder_type::Symbol = :none,
    Ws::AbstractVector{<:Real} = [0.0],
    n_disorder_realisations::Int=1,
    x0_sym::Union{Symbol, Real} = :centre, 
    y0_sym::Union{Symbol, Real} = :centre, 
    Es::AbstractVector{<:Real} = [0.0],
    kappas::AbstractVector{<:Real} = [1.0], 
    orbital_displacements::AbstractVector{<:Real} = [0.0],
    phis::AbstractVector{<:Real} = [0.0],
    periodic_x::Bool=false, 
    periodic_y::Bool=false, 
    n_lowest_evals::Int=10,
    keep_square::Bool=false,
    sweep_x0y0::Tuple{Bool, Int}=(false, 51),
    use_coo_assembly::Bool=true,  # NEW: Toggle optimized assembly
    scale_kappa_to_L::Bool=false,  # NEW: Whether to scale kappa inversely with system size
    kappa_scales::AbstractVector{<:Real} = [0.0004]  # NEW: Scaling factors for kappa (κ = scale/L)
)::DataFrame

    RowType = NamedTuple{(:A, :B, :m, :gamma, :W, :Lx, :Ly, :x, :y, :E, :kappa, :d, :phi, :low_lying_evals),
                         Tuple{Float64,Float64,Float64,Float64,Float64,Int64,Int64,Float64,Float64,Float64,Float64,Float64,Float64, Vector{Float64}}}
    
    all_rows = Vector{RowType}()
    do_sweep, n_points = sweep_x0y0

    for Lx in Lxs, Ly in Lys
        if keep_square && Lx != Ly
            continue
        end

        if do_sweep
            x_center = Lx / 2.0
            y_center = Ly / 2.0
            centre_range = 1.0
            x0_vals = collect(range(x_center - centre_range, x_center + centre_range, length=n_points))
            y0_vals = collect(range(y_center - centre_range, y_center + centre_range, length=n_points))
        else
            x0_fixed = (x0_sym == :centre || x0_sym == :center) ? Lx / 2.0 : Float64(x0_sym)
            y0_fixed = (y0_sym == :centre || y0_sym == :center) ? Ly / 2.0 : Float64(y0_sym)
            x0_vals = [x0_fixed]
            y0_vals = [y0_fixed]
        end

        # Determine kappas based on scaling flag
        effective_kappas = if scale_kappa_to_L
            [scale / Lx for scale in kappa_scales]
        else
            kappas
        end

        param_grid = collect(Iterators.product(Avals, Bvals, mvals, gammas, Ws, Es, effective_kappas, orbital_displacements, phis, x0_vals, y0_vals))
        n_params = length(param_grid)
        batch_rows = Vector{RowType}(undef, n_params)

        # OPTIMIZATION: 2D flat parallelization (param × realisation)
        # Pre-allocate accumulators for eigenvalues (one per parameter)
        evals_accum = [zeros(Float64, n_lowest_evals) for _ in 1:n_params]
        accum_lock = ReentrantLock()  # Thread-safe access to accumulators

        # Flatten to 2D grid: (param_idx, realisation)
        param_real_grid = collect(Iterators.product(1:n_params, 1:n_disorder_realisations))
        n_total_tasks = length(param_real_grid)

        # Progress bar with 2D parallelization (scales better with many realisations)
        @showprogress dt=1.0 Threads.@threads for idx in 1:n_total_tasks
            i, realisation = param_real_grid[idx]
            A, B, m, gamma, W, E, kappa, d, phi, x0, y0 = param_grid[i]

            # Compute Hamiltonian and eigenvalues
            if use_coo_assembly
                H = real_space_perturbed_disordered_hamiltonian_qwz_coo(Lx, Ly;
                        A=A, B=B, m=m, gamma=gamma, 
                        perturbation_type=perturbation_type, 
                        disorder_type=disorder_type, W=W,
                        periodic_x=periodic_x, periodic_y=periodic_y)
            else
                H = real_space_perturbed_disordered_hamiltonian_qwz(Lx, Ly;
                        A=A, B=B, m=m, gamma=gamma, 
                        perturbation_type=perturbation_type, 
                        disorder_type=disorder_type, W=W,
                        periodic_x=periodic_x, periodic_y=periodic_y, 
                        sparse_output=true)
            end

            # Build position operators with embedding parameters
            X_embedded, Y_embedded = build_position_operators(Lx, Ly; orbital_displacement=d, phi=phi)
            L = build_localiser_operator(H, X_embedded, Y_embedded, x0, y0, E; kappa=kappa)

            # Create parameter tuple for logging
            param_info = (A=A, B=B, m=m, gamma=gamma, W=W, E=E, kappa=kappa, d=d, phi=phi, x0=x0, y0=y0, Lx=Lx, Ly=Ly)
            
            low_lying_evals, eval_status = compute_low_lying_localiser_shift_invert(
                L; 
                n_eigenpairs=n_lowest_evals, 
                kk_tol=1e-8, 
                kk_maxiter=300,
                params=param_info
            )
            
            # Thread-safe accumulation: minimal lock time
            lock(accum_lock) do
                for j in 1:min(length(low_lying_evals), n_lowest_evals)
                    evals_accum[i][j] += low_lying_evals[j]
                end
            end
        end

        # Average results for each parameter across all realisations
        for i in 1:n_params
            A, B, m, gamma, W, E, kappa, d, phi, x0, y0 = param_grid[i]
            low_lying_evals_avg = evals_accum[i] ./ n_disorder_realisations

            batch_rows[i] = (
                A=Float64(A), B=Float64(B), m=Float64(m), gamma=Float64(gamma), W=Float64(W), 
                Lx=Int64(Lx), Ly=Int64(Ly), x=Float64(x0), y=Float64(y0), E=Float64(E), 
                kappa=Float64(kappa), d=Float64(d), phi=Float64(phi), low_lying_evals=low_lying_evals_avg
            )
        end
        
        append!(all_rows, batch_rows)
    end

    return DataFrame(all_rows)
end


#=
================================================================================
 Efficient signature / local Chern marker via LDL^dagger factorization
 for the 2D spectral localizer.

 Reference: A. Cerjan & T. A. Loring, "Tutorial: Classifying Photonic
 Topology Using the Spectral Localizer and Numerical K-Theory,"
 arXiv:2411.03515 (APL Photonics 9, 111102 (2024)).
 Relevant equations: (12)-(13) [2D localizer & local Chern marker],
 Sec. IV D, Eqs. (41)-(44) [LDL^dagger, Sylvester's law, complex embedding].

 This file assumes it is `include`-d AFTER your main script, so that
 `build_localiser_operator`, `build_position_operators`, etc. are already
 defined. Alternatively, paste the functions below directly into your
 script above `fast_low_lying_localiser_spectrum`.
================================================================================
=#

# ------------------------------------------------------------------
    # 1. Robust sparse signature of an invertible complex Hermitian
    #    matrix via LDL^dagger + Sylvester's law of inertia.
    #
    #    sig[M] = sig[D] = (#positive diagonal entries) - (#negative diagonal entries)
    #
    #    Uses Julia's native sparse LDLt (SparseArrays.CHOLMOD), which DOES
    #    support complex Hermitian input directly. CHOLMOD's LDLt uses only
    #    a fill-reducing symbolic permutation (AMD by default) and NO dynamic
    #    numerical pivoting, so on rare, non-generic (x,E) points it can hit
    #    a numerically-zero pivot and throw, even though the matrix itself is
    #    invertible. We handle that with a tiny, escalating diagonal shift
    #    (standard practice for CHOLMOD LDLt breakdowns), and fall back to a
    #    dense eigendecomposition only if that also fails (this should be a
    #    rare, isolated event, not the common path).
# ------------------------------------------------------------------
function matrix_signature_ldlt(
    Mherm::Hermitian{ComplexF64,<:SparseMatrixCSC{ComplexF64,Int}};
    max_shift_tries::Int = 6,
    base_shift::Float64 = 1e-12
)::Int
    n = size(Mherm, 1)
    shift = 0.0
    for attempt in 0:max_shift_tries
        try
            F = ldlt(Mherm; shift = shift)
            # F.D alone is not reliably convertible with `sparse()` across
            # Julia versions; F.LD (packed unit-lower-triangular L below the
            # diagonal, D ON the diagonal) is the robust way to recover D,
            # since L has an implicit unit diagonal.
            LDmat = sparse(F.LD)
            Dvals = real.(diag(LDmat))
            n_pos = count(x -> x > 0, Dvals)
            n_neg = count(x -> x < 0, Dvals)
            n_zero = n - n_pos - n_neg
            if n_zero > 0 && shift == 0.0
                @warn "matrix_signature_ldlt: $(n_zero) numerically-zero pivot(s) " *
                      "found -- (x,E) is very close to a local topological transition " *
                      "(local gap ≈ 0). Signature may be unreliable at this point."
            end
            return n_pos - n_neg
        catch err
            estr = string(typeof(err))
            if err isa LinearAlgebra.ZeroPivotException || occursin("CHOLMOD", estr) || occursin("Pivot", estr)
                shift = shift == 0.0 ? base_shift : shift * 10
                continue
            else
                rethrow(err)
            end
        end
    end
    @warn "matrix_signature_ldlt: sparse LDLt failed after $(max_shift_tries) " *
          "regularization attempts; falling back to a dense eigendecomposition " *
          "for this single (x,E) point (n=$(n)). This should be rare."
    evs = eigvals(Hermitian(Matrix(Mherm)))
    return count(x -> x > 0, evs) - count(x -> x < 0, evs)
end

# ------------------------------------------------------------------
    # 2. Portable alternative: real embedding of a complex Hermitian matrix
    #    (paper Eqs. 43-44), for use with a real-only sparse LDLt backend,
    #    or if complex CHOLMOD support proves flaky in your Julia install.
    #
    #    For Hermitian M = A + iB (A = Re(M) symmetric, B = Im(M) antisymmetric),
    #
    #        M_R = [ A  -B ]      is REAL SYMMETRIC, size 2n x 2n,
    #              [ B   A ]
    #
    #    and sig(M) = (1/2) sig(M_R). (Proof: if Mv = lambda*v with v = p+iq,
    #    then M_R*(p,q)^T = lambda*(p,q)^T, and M_R*(-q,p)^T = lambda*(-q,p)^T
    #    independently, i.e., every eigenvalue of M appears exactly TWICE in
    #    the spectrum of M_R.)
# ------------------------------------------------------------------
function real_embed_hermitian(M::SparseMatrixCSC{ComplexF64,Int})::SparseMatrixCSC{Float64,Int}
    A = real.(M)
    B = imag.(M)
    return [A  -B;
            B   A]
end

function matrix_signature_ldlt_real_embed(
    M::SparseMatrixCSC{ComplexF64,Int};
    max_shift_tries::Int = 6,
    base_shift::Float64 = 1e-12
)::Int
    Mr = Symmetric(real_embed_hermitian(M))
    n = size(Mr, 1)
    shift = 0.0
    for attempt in 0:max_shift_tries
        try
            F = ldlt(Mr; shift = shift)
            LDmat = sparse(F.LD)
            Dvals = diag(LDmat)
            n_pos = count(x -> x > 0, Dvals)
            n_neg = count(x -> x < 0, Dvals)
            sig_embed = n_pos - n_neg              # = sig(M_R)
            return div(sig_embed, 2)                # sig(M) = (1/2) sig(M_R)
        catch err
            estr = string(typeof(err))
            if err isa LinearAlgebra.ZeroPivotException || occursin("CHOLMOD", estr) || occursin("Pivot", estr)
                shift = shift == 0.0 ? base_shift : shift * 10
                continue
            else
                rethrow(err)
            end
        end
    end
    @warn "matrix_signature_ldlt_real_embed: sparse LDLt failed; falling back " *
          "to a dense eigendecomposition (n=$(n))."
    evs = eigvals(Mr)
    sig_embed = count(x -> x > 0, evs) - count(x -> x < 0, evs)
    return div(sig_embed, 2)
end

# ------------------------------------------------------------------
    # 3. Local Chern marker for the 2D spectral localizer, Eq. (13):
    #        C_(x,E)^L = (1/2) sig[L_(x,E)^(2D)]
    #    Reuses `build_localiser_operator`, already defined in your script.
# ------------------------------------------------------------------
function local_chern_marker_2d(
    H::AbstractMatrix{ComplexF64},
    X::Diagonal, Y::Diagonal,
    x0::Real, y0::Real, E::Real;
    kappa::Real = 1.0,
    method::Symbol = :complex   # :complex (native, faster) or :real_embed (portable fallback)
)::Tuple{Float64,Int}
    L = build_localiser_operator(H, X, Y, x0, y0, E; kappa = kappa)
    sig = method === :complex ?
        matrix_signature_ldlt(Hermitian(L)) :
        matrix_signature_ldlt_real_embed(L)
    return (sig / 2.0, sig)
end

# ------------------------------------------------------------------
    # 4. Drop-in extension of your `fast_low_lying_localiser_spectrum`
    #    that ALSO computes the local Chern marker at each grid point,
    #    reusing the same `L` that's already built for the low-lying
    #    spectrum calculation (so there is no duplicated matrix assembly,
    #    only one extra sparse factorization per point). Changes relative
    #    to your original function are marked "### NEW".
# ------------------------------------------------------------------
function fast_low_lying_localiser_spectrum_with_chern(
    Lxs::AbstractVector{<:Integer},
    Lys::AbstractVector{<:Integer};
    Avals::AbstractVector{<:Real} = [1.0],
    Bvals::AbstractVector{<:Real} = [1.0],
    mvals::AbstractVector{<:Real} = [0.0],
    gammas::AbstractVector{<:Real} = [0.0],
    perturbation_type::Symbol = :none,
    disorder_type::Symbol = :none,
    Ws::AbstractVector{<:Real} = [0.0],
    n_disorder_realisations::Int = 1,
    x0_sym::Union{Symbol,Real} = :centre,
    y0_sym::Union{Symbol,Real} = :centre,
    Es::AbstractVector{<:Real} = [0.0],
    kappas::AbstractVector{<:Real} = [1.0],
    orbital_displacements::AbstractVector{<:Real} = [0.0],
    phis::AbstractVector{<:Real} = [0.0],
    periodic_x::Bool = false,
    periodic_y::Bool = false,
    n_lowest_evals::Int = 10,
    keep_square::Bool = false,
    sweep_x0y0::Tuple{Bool,Int} = (false, 51),
    use_coo_assembly::Bool = true,
    scale_kappa_to_L::Bool = false,
    kappa_scales::AbstractVector{<:Real} = [0.0004],
    compute_chern::Bool = true, 
    chern_method::Symbol = :complex # (:complex or :real_embed)
)::DataFrame

    RowType = NamedTuple{
        (:A, :B, :m, :gamma, :W, :Lx, :Ly, :x, :y, :E, :kappa, :d, :phi, :low_lying_evals, :chern_marker, :signature, :eval_status), # ### NEW: eval_status
        Tuple{Float64,Float64,Float64,Float64,Float64,Int64,Int64,Float64,Float64,Float64,Float64,Float64,Float64,
              Vector{Float64},Float64,Int64,Symbol}                                                                 # ### NEW: Symbol
    }

    all_rows = Vector{RowType}()
    do_sweep, n_points = sweep_x0y0

    for Lx in Lxs, Ly in Lys
        if keep_square && Lx != Ly
            continue
        end

        X, Y = build_position_operators(Lx, Ly)

        if do_sweep
            x_center = Lx / 2.0
            y_center = Ly / 2.0
            centre_range = 1.0
            x0_vals = collect(range(x_center - centre_range, x_center + centre_range, length=n_points))
            y0_vals = collect(range(y_center - centre_range, y_center + centre_range, length=n_points))
        else
            x0_fixed = (x0_sym == :centre || x0_sym == :center) ? Lx / 2.0 : Float64(x0_sym)
            y0_fixed = (y0_sym == :centre || y0_sym == :center) ? Ly / 2.0 : Float64(y0_sym)
            x0_vals = [x0_fixed]
            y0_vals = [y0_fixed]
        end

        # Determine kappas based on scaling flag
        effective_kappas = if scale_kappa_to_L
            [scale / Lx for scale in kappa_scales]
        else
            kappas
        end

        param_grid = collect(Iterators.product(Avals, Bvals, mvals, gammas, Ws, Es, effective_kappas, orbital_displacements, phis, x0_vals, y0_vals))
        n_params = length(param_grid)
        batch_rows = Vector{RowType}(undef, n_params)

        evals_accum = [zeros(Float64, n_lowest_evals) for _ in 1:n_params]
        chern_accum = zeros(Float64, n_params)          # ### NEW
        sig_accum   = zeros(Int, n_params)               # ### NEW (last-realisation signature, for inspection)
        eval_status_accum = fill(:unknown, n_params)     # ### NEW: track eigenvalue solver status
        accum_lock = ReentrantLock()

        param_real_grid = collect(Iterators.product(1:n_params, 1:n_disorder_realisations))
        n_total_tasks = length(param_real_grid)

        @showprogress dt=1.0 Threads.@threads for idx in 1:n_total_tasks
            i, realisation = param_real_grid[idx]
            A, B, m, gamma, W, E, kappa, d, phi, x0, y0 = param_grid[i]

            H = if use_coo_assembly
                real_space_perturbed_disordered_hamiltonian_qwz_coo(Lx, Ly;
                        A=A, B=B, m=m, gamma=gamma,
                        perturbation_type=perturbation_type,
                        disorder_type=disorder_type, W=W,
                        periodic_x=periodic_x, periodic_y=periodic_y)
            else
                real_space_perturbed_disordered_hamiltonian_qwz(Lx, Ly;
                        A=A, B=B, m=m, gamma=gamma,
                        perturbation_type=perturbation_type,
                        disorder_type=disorder_type, W=W,
                        periodic_x=periodic_x, periodic_y=periodic_y,
                        sparse_output=true)
            end

            # Build position operators with embedding parameters
            X_embedded, Y_embedded = build_position_operators(Lx, Ly; orbital_displacement=d, phi=phi)
            L = build_localiser_operator(H, X_embedded, Y_embedded, x0, y0, E; kappa=kappa)

            # Create parameter tuple for logging
            param_info = (A=A, B=B, m=m, gamma=gamma, W=W, E=E, kappa=kappa, d=d, phi=phi, x0=x0, y0=y0, Lx=Lx, Ly=Ly)

            # ### NEW: Robust eigenvalue solver with parameter logging
            low_lying_evals, eval_status = compute_low_lying_localiser_shift_invert(
                L;
                n_eigenpairs=n_lowest_evals,
                kk_tol=1e-8,
                kk_maxiter=300,
                params=param_info
            )

            # ### NEW: signature / local Chern marker, reusing the same L
            chern_val = 0.0
            sig_val = 0
            if compute_chern
                sig_val = chern_method === :complex ?
                    matrix_signature_ldlt(Hermitian(L)) :
                    matrix_signature_ldlt_real_embed(L)
                chern_val = sig_val / 2.0
            end

            lock(accum_lock) do
                for j in 1:min(length(low_lying_evals), n_lowest_evals)
                    evals_accum[i][j] += low_lying_evals[j]
                end
                chern_accum[i] += chern_val
                sig_accum[i] = sig_val
                eval_status_accum[i] = eval_status
            end
        end

        for i in 1:n_params
            A, B, m, gamma, W, E, kappa, d, phi, x0, y0 = param_grid[i]
            low_lying_evals_avg = evals_accum[i] ./ n_disorder_realisations
            chern_marker_avg = chern_accum[i] / n_disorder_realisations   # ### NEW

            batch_rows[i] = (
                A=Float64(A), B=Float64(B), m=Float64(m), gamma=Float64(gamma), W=Float64(W),
                Lx=Int64(Lx), Ly=Int64(Ly), x=Float64(x0), y=Float64(y0), E=Float64(E),
                kappa=Float64(kappa), d=Float64(d), phi=Float64(phi), low_lying_evals=low_lying_evals_avg,
                chern_marker=Float64(chern_marker_avg), signature=Int64(sig_accum[i]), 
                eval_status=eval_status_accum[i]  # ### NEW: track eigenvalue solver status
            )
        end

        append!(all_rows, batch_rows)
    end

    return DataFrame(all_rows)
end

# ------------------------------------------------------------------
    # 5. Quick benchmark: LDLt signature vs. the existing shift-invert
    #    low-lying-spectrum method, at fixed system size and (x,E,kappa).
    #    Both are dominated by ONE sparse factorization of the same N x N
    #    localiser, so timings should be of the same order; LDLt should be
    #    the cheaper of the two (it does not also run Krylov iterations).
# ------------------------------------------------------------------
function benchmark_signature_vs_spectrum(Lx=30, Ly=30; n_reps=5)
    A, B, m, gamma, W = 1.0, 1.0, -1.0, 0.75, 0.05
    kappa, E = 0.02, 0.0

    X, Y = build_position_operators(Lx, Ly)
    H = real_space_perturbed_disordered_hamiltonian_qwz_coo(Lx, Ly;
            A=A, B=B, m=m, gamma=gamma, perturbation_type=:symmetric,
            disorder_type=:anderson, W=W)
    L = build_localiser_operator(H, X, Y, Lx/2, Ly/2, E; kappa=kappa)
    n = size(L, 1)
    println("Localiser size: $(n) x $(n)  (nnz = $(nnz(L)))")

    t_sig = @elapsed for _ in 1:n_reps
        matrix_signature_ldlt(Hermitian(L))
    end
    t_sig /= n_reps

    t_spec = @elapsed for _ in 1:n_reps
        compute_low_lying_localiser_shift_invert(L; n_eigenpairs=10, kk_tol=1e-8, kk_maxiter=300)
    end
    t_spec /= n_reps

    sig = matrix_signature_ldlt(Hermitian(L))
    println("Signature sig[L]        = $(sig)")
    println("Local Chern marker C    = $(sig/2)")
    println("LDLt signature time     : $(round(t_sig*1000, digits=2)) ms")
    println("Shift-invert spectrum   : $(round(t_spec*1000, digits=2)) ms")
    return (t_sig, t_spec, sig)
end







# Execution Setup
Avals = [1.0]
Bvals = [1.0]
mvals = [-1.0] #collect(-4.5:0.5:0.5)
gammavals = [0.5] #collect(-3.0:0.05:3.0) #[0.75, 0.5, 1.0, 2.0, 3.0] #collect(0.0:0.1:3.0) #collect(range(0.0, 2.0, 51))
perturbation_type = :sym_cos_sum
embedding_phis = [0.0]
embedding_ds = [0.0] #collect(0.0:0.5:5.0) #[0.0, 1.0, 10.0]
disorder_type = :none
Ws = [0.0] #[0.01, 0.1] #collect(range(0.0, 1.0, 51))
n_disorder_realisations = 1
Lxs = collect(10:10:500) #[14]
Lys = collect(10:10:500) #[14]
keep_square = true
xs = :centre
ys = :centre
sweep_x0y0 = (false, 51)
Es = [1.3] #collect(-3.0:0.05:3.0) #collect(range(-1.125, -1.075, length=101))
kappas = [0.2] #logrange(1e-3, 1e-0, 3) #collect(range(1e-3, 1e-0, 30)) #[2e-1] ## ((0.02 works well at L=50, gives kappa_scales=0.0004))
scale_kappa_to_L = true
kappa_scales = [0.0004] ## kappa = kappa_scale / L
n_lowest_evals = 4

results_df = fast_low_lying_localiser_spectrum_with_chern(
    Lxs, Lys; 
    Avals=Avals, Bvals=Bvals, mvals=mvals, gammas=gammavals,
    perturbation_type=perturbation_type, disorder_type=disorder_type, Ws=Ws,
    x0_sym=xs, y0_sym=ys, Es=Es, kappas=kappas, 
    phis=embedding_phis, orbital_displacements=embedding_ds,
    periodic_x=false, periodic_y=false,
    n_disorder_realisations=n_disorder_realisations,
    n_lowest_evals=n_lowest_evals,
    keep_square=keep_square,
    scale_kappa_to_L=scale_kappa_to_L,
    kappa_scales=kappa_scales,
    compute_chern=true,
    chern_method=:complex
)

println("Computed spectral localiser DataFrame with $(nrow(results_df)) rows and $(ncol(results_df)) columns.")

foldername = "data/qwz_disorder_lowspecloc/low_specloc_spectrum_nEvals$(n_lowest_evals)_$(perturbation_type)_disord_$(disorder_type)_nAvgReals$(n_disorder_realisations)_PARAMS_As$(Avals[1])-$(length(Avals))-$(Avals[end])_Bs$(Bvals[1])-$(Bvals[end])_ms$(mvals[1])-$(mvals[end])_gammas$(gammavals[1])-$(gammavals[end])_IndTyp$(perturbation_type)_Ws$(Ws[1])-$(Ws[end])_Lxs$(Lxs[1])-$(Lxs[end])_Lys$(Lys[1])-$(Lys[end])_keptSquare$(keep_square)_x0y0sweep$(sweep_x0y0[1])_Es$(Es[1])-$(Es[end])_kappas$(kappas[1])-$(kappas[end])_embedPhis$(embedding_phis[1])-$(embedding_phis[end])_embedDs$(embedding_ds[1])-$(embedding_ds[end])"
isdir(foldername) || mkpath(foldername)
filename = joinpath(foldername, "data.jld2")

println("Saving DataFrame to $filename ...")
t_start = time()
@save filename results_df
println("Saved DataFrame to $filename, took $(round(time() - t_start, digits=2)) seconds to save.")
