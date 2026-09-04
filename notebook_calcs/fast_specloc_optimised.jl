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

# function build_position_operators(Lx::Int, Ly::Int; 
#     orbital_displacement::Real=0.0, 
#     phi::Real=0.0
# )
#     nsites = Lx * Ly
#     dim = 2 * nsites
#     xvec = zeros(Float64, dim)
#     yvec = zeros(Float64, dim)
    
#     # Compute embedding displacement components
#     dx_disp = orbital_displacement * cos(phi)
#     dy_disp = orbital_displacement * sin(phi)
    
#     for y in 1:Ly, x in 1:Lx
#         base = 2 * ((y - 1) * Lx + (x - 1))
#         # Orbital 1 (sublattice A): displaced by +d/2
#         xvec[base + 1] = x + dx_disp / 2
#         yvec[base + 1] = y + dy_disp / 2
#         # Orbital 2 (sublattice B): displaced by -d/2
#         xvec[base + 2] = x - dx_disp / 2
#         yvec[base + 2] = y - dy_disp / 2
#     end
#     X = Diagonal(ComplexF64.(xvec))
#     Y = Diagonal(ComplexF64.(yvec))
#     return X, Y
# end

## the OG
function build_position_operators(
    Lx::Int, 
    Ly::Int;
    orbital_displacement::Real=0.0,
    phi::Real=0.0
)

    nsites = Lx * Ly
    dim = 2 * nsites
    xvec = zeros(Float64, dim)
    yvec = zeros(Float64, dim)
    for y in 1:Ly, x in 1:Lx
        base = 2 * ((y - 1) * Lx + (x - 1))
        xvec[base + 1] = x
        xvec[base + 2] = x
        yvec[base + 1] = y
        yvec[base + 2] = y
    end
    X = Diagonal(ComplexF64.(xvec))
    Y = Diagonal(ComplexF64.(yvec))
    return X, Y
end

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

function compute_low_lying_localiser_shift_invert(
    L::SparseMatrixCSC{ComplexF64, Int};
    n_eigenpairs::Int=10,
    kk_tol::Real=1e-8,
    kk_maxiter::Int=300
)::Vector{Float64}
    n = size(L, 1)
    
    # Compute LU factorization once for shift-and-invert action L^{-1} x
    F = lu(L)
    L_inv_action(x) = F \ x

    # Target largest magnitude eigenvalues of L^-1 (corresponds to smallest of L)
    # Since L is Hermitian, L^-1 is also Hermitian
    evals_inv, evecs, info = eigsolve(
        L_inv_action, 
        rand(ComplexF64, n), 
        n_eigenpairs, 
        :LM; 
        ishermitian=true, 
        tol=kk_tol, 
        maxiter=kk_maxiter
    )
    
    # Real eigenvalues of L are 1 / (eigenvalues of L^-1)
    signed_evals = real.(1.0 ./ evals_inv)
    
    sort!(signed_evals, by=abs)
    return signed_evals[1:min(n_eigenpairs, length(signed_evals))]
end

function real_space_perturbed_qwz_blocks(; 
    A::Real=1.0, B::Real=1.0, m::Real=0.0, gamma::Real=0.0, perturbation_type::Symbol=:none
)
    onsite = ComplexF64.((m + 2.0 * B) * sigma_z)
    tx = ComplexF64.(-0.5 * B * sigma_z - 0.5im * A * sigma_x)
    ty = ComplexF64.(-0.5 * B * sigma_z - 0.5im * A * sigma_y)

    # Symmetric perturbations (cosine-like, real scalar)
    if perturbation_type == :sym_cos_sum || perturbation_type == :symmetric
        # Symmetric addition: uniform real perturbation
        tx += ComplexF64.(0.5 * gamma * identity_2)
        ty += ComplexF64.(0.5 * gamma * identity_2)
    elseif perturbation_type == :sym_cos_diff
        # Symmetric difference: different signs in x and y
        tx += ComplexF64.(0.5 * gamma * identity_2)
        ty += ComplexF64.(-0.5 * gamma * identity_2)
    elseif perturbation_type == :sym_cos_add
        # Combined symmetric: average contribution
        tx += ComplexF64.(0.25 * gamma * identity_2)
        ty += ComplexF64.(0.25 * gamma * identity_2)
    elseif perturbation_type == :sym_cos_sub
        # Combined symmetric difference
        tx += ComplexF64.(0.25 * gamma * identity_2)
        ty += ComplexF64.(-0.25 * gamma * identity_2)
    # Antisymmetric perturbations (sine-like, imaginary component)
    elseif perturbation_type == :asym_sin_sum || perturbation_type == :antisymmetric
        # Antisymmetric addition: imaginary perturbation
        tx += ComplexF64.(0.5im * gamma * identity_2)
        ty += ComplexF64.(0.5im * gamma * identity_2)
    elseif perturbation_type == :asym_sin_diff
        # Antisymmetric difference: opposite imaginary contributions
        tx += ComplexF64.(0.5im * gamma * identity_2)
        ty += ComplexF64.(-0.5im * gamma * identity_2)
    elseif perturbation_type == :asym_sin_add
        # Combined antisymmetric
        tx += ComplexF64.(0.25im * gamma * identity_2)
        ty += ComplexF64.(0.25im * gamma * identity_2)
    elseif perturbation_type == :asym_sin_sub
        # Combined antisymmetric difference
        tx += ComplexF64.(0.25im * gamma * identity_2)
        ty += ComplexF64.(-0.25im * gamma * identity_2)
    elseif perturbation_type != :none
        error("Unknown perturbation type: $perturbation_type")
    end

    return onsite, tx, ty
end

site_index_qwz(x::Int, y::Int, orb::Int, Lx::Int, Ly::Int) = 2 * ((y - 1) * Lx + (x - 1)) + orb

# OPTIMIZED: Use COO format for sparse assembly (5-10x faster than element-wise +=)
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
    scale_kappa_to_L::Bool=false,  # NEW: Whether to scale kappa by system size
    kappa_scales::AbstractVector{<:Real} = [0.0004]  # NEW: Scaling factors for kappa
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
            [scale * Lx for scale in kappa_scales]
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

            low_lying_evals = compute_low_lying_localiser_shift_invert(
                L; 
                n_eigenpairs=n_lowest_evals, 
                kk_tol=1e-8, 
                kk_maxiter=300
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




# Execution Setup
Avals = [1.0]
Bvals = [1.0]
mvals = [-1.0]
gammavals = collect(-3.0:0.05:3.0) #[0.75, 0.5, 1.0, 2.0, 3.0] #collect(0.0:0.1:3.0) #collect(range(0.0, 2.0, 51))
perturbation_type = :asym_sin_sum
embedding_phis = [0.0]
embedding_ds = [0.0]#, 1.0, 10.0]
disorder_type = :none
Ws = [0.0] #[0.01, 0.1] #collect(range(0.0, 1.0, 51))
n_disorder_realisations = 1
Lxs = [20] #collect(10:10:200) #[14]
Lys = [20] #collect(10:10:200) #[14]
keep_square = true
xs = :centre
ys = :centre
sweep_x0y0 = (false, 51)
Es = collect(-3.0:0.05:3.0) #collect(range(-1.125, -1.075, length=101))
kappas = [0.05] #logrange(1e-3, 1e-0, 3) #collect(range(1e-3, 1e-0, 30)) #[2e-1] ## ((0.02 works well at L=50, gives kappa_scales=0.0004))
scale_kappa_to_L = false
kappa_scales = [0.0004] ## kappa = kappa_scale * L
n_lowest_evals = 4

results_df = fast_low_lying_localiser_spectrum(
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
    kappa_scales=kappa_scales
)

println("Computed spectral localiser DataFrame with $(nrow(results_df)) rows and $(ncol(results_df)) columns.")

foldername = "data/qwz_disorder_lowspecloc/low_specloc_spectrum_nEvals$(n_lowest_evals)_$(perturbation_type)_disord_$(disorder_type)_nAvgReals$(n_disorder_realisations)_PARAMS_As$(Avals[1])-$(length(Avals))-$(Avals[end])_Bs$(Bvals[1])-$(Bvals[end])_ms$(mvals[1])-$(mvals[end])_gammas$(gammavals[1])-$(gammavals[end])_Ws$(Ws[1])-$(Ws[end])_Lxs$(Lxs[1])-$(Lxs[end])_Lys$(Lys[1])-$(Lys[end])_keptSquare$(keep_square)_x0y0sweep$(sweep_x0y0[1])_Es$(Es[1])-$(Es[end])_kappas$(kappas[1])-$(kappas[end])_embedPhis$(embedding_phis[1])-$(embedding_phis[end])_embedDs$(embedding_ds[1])-$(embedding_ds[end])"
isdir(foldername) || mkpath(foldername)
filename = joinpath(foldername, "data.jld2")

println("Saving DataFrame to $filename ...")
t_start = time()
@save filename results_df
println("Saved DataFrame to $filename, took $(round(time() - t_start, digits=2)) seconds to save.")
