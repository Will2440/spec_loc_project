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
using LDLFactorizations

println("Julia Threads: ", Threads.nthreads())
println("BLAS Threads:  ", LinearAlgebra.BLAS.get_num_threads())


## Pauli matrices
sigma_x = [0 1; 1 0]
sigma_y = [0 -im; im 0]
sigma_z = [1 0; 0 -1]
identity = [1 0; 0 1]

# Build position operators (x and y) that act on the full 2-orbital per site basis.
function build_position_operators(
    Lx::Int, 
    Ly::Int
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

# Compute the signature (n_pos - n_neg) and also return the minimum absolute eigenvalue.
function spectral_localiser_signature(
    H::AbstractMatrix{ComplexF64},
    X::Diagonal,
    Y::Diagonal,
    x0::Real,
    y0::Real,
    E::Real; 
    kappa::Real=1.0,
    zero_tol::Real=1e-12
)::Tuple{Int, Float64}

    D = size(H, 1)
    I_D = Matrix{ComplexF64}(I, D, D)
    Xmat = Matrix(X)
    Ymat = Matrix(Y)

    Hshift = H - E * I_D
    Ablock = kappa * (Xmat - x0 * I_D)
    Bblock = kappa * (Ymat - y0 * I_D)

    L = [Hshift  Ablock - im * Bblock;
         Ablock + im * Bblock  -Hshift]

    vals = eigvals(Hermitian(L))
    revals = real(vals)
    npos = count(>(zero_tol), revals)
    nneg = count(<(-zero_tol), revals)
    minabs = minimum(abs.(revals))
    signature = npos - nneg
    return signature, minabs
end

# Compute the Chern number from the localiser signature: C = 1/2 * signature.
function chern_from_localiser_signature(
    signature::Int, 
    minabs::Real; 
    zero_tol::Real=1e-12, 
    chern_tol::Real=1e-6
)
    if minabs <= zero_tol
        return missing
    end
    chern_d = signature / 2
    chern_round = round(Int, chern_d)
    if abs(chern_d - chern_round) > chern_tol
        return missing
    end
    return Int8(chern_round)
end


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

# # 3. Updated Spectral Localizer Wrapper
## Threaded
function slow_compute_perturbed_spectral_localiser_df(
    Lx::Int, 
    Ly::Int; 
    Avals::AbstractVector{<:Real} = [1.0], 
    Bvals::AbstractVector{<:Real} = [1.0], 
    mvals::AbstractVector{<:Real} = [0.0],
    gammas::AbstractVector{<:Real} = [0.0],
    perturbation_type::Symbol = :none,
    xs::AbstractVector = collect(1:Lx), 
    ys::AbstractVector = collect(1:Ly), 
    Es::AbstractVector = [0.0],
    kappas::AbstractVector{<:Real} = [1.0], 
    periodic_x::Bool=false, 
    periodic_y::Bool=false, 
    sparse_output::Bool=false
)::DataFrame

    # 1. Pre-build position operators
    X, Y = build_position_operators(Lx, Ly)

    # 2. Collect the parameter grid into a Vector so Threads.@threads can index into it
    grid = collect(Iterators.product(Avals, Bvals, mvals, gammas, xs, ys, Es, kappas))
    n_iters = length(grid)

    # 3. Define element type and pre-allocate output array (CRITICAL for thread safety)
    RowType = NamedTuple{(:A, :B, :m, :gamma, :Lx, :Ly, :x, :y, :E, :kappa, :localiser_gap, :signature, :chern_number),
                         Tuple{Float64,Float64,Float64,Float64,Int64,Int64,Float64,Float64,Float64,Float64,Float64,Int64,Union{Missing, Int8}}}
    
    rows = Vector{RowType}(undef, n_iters)

    # 4. Multi-threaded loop with ProgressMeter
    @showprogress dt=0.1 Threads.@threads for i in 1:n_iters
        A, B, m, gamma, x0, y0, E, kappa = grid[i]

        # Build perturbed H for this parameter set
        H = real_space_perturbed_hamiltonian_qwz(Lx, Ly; 
                A=A, B=B, m=m, gamma=gamma, 
                perturbation_type=perturbation_type, 
                periodic_x=periodic_x, periodic_y=periodic_y, sparse_output=sparse_output)
        
        sig, minabs = spectral_localiser_signature(H, X, Y, x0, y0, E; kappa=kappa)
        chern = chern_from_localiser_signature(sig, minabs; zero_tol=1e-12, chern_tol=1e-6)
        
        # Thread-safe assignment by index
        rows[i] = (A=Float64(A), B=Float64(B), m=Float64(m), gamma=Float64(gamma), Lx=Lx, Ly=Ly, 
                   x=Float64(x0), y=Float64(y0), E=Float64(E), kappa=Float64(kappa), 
                   localiser_gap=Float64(minabs), signature=Int64(sig), chern_number=chern)
    end

    return DataFrame(rows)
end



## search energy and space at good kappa

Avals = [1.0]
Bvals = [1.0]
mvals = [1.0, -1.0]
gammavals = [0.0, 2.0]
perturbation_type=:symmetric
Lx = 20
Ly = 20
xs = collect(range(-1, Lx+1, length=25))
ys = collect(range(-1, Ly+1, length=25)) #collect(1:Ly)
Es = collect(range(0.0, 0.0, length=1))
kappas = [2e-1]#logrange(1e-3, 5e1, 25)

println("--------- Computation Parameters: ----------")
println("Avals: ", minimum(Avals), " - ", length(Avals), " - ", maximum(Avals))
println("Bvals: ", minimum(Bvals), " - ", length(Bvals), " - ", maximum(Bvals))
println("mvals: ", minimum(mvals), " - ", length(mvals), " - ", maximum(mvals))
println("gammavals: ", minimum(gammavals), " - ", length(gammavals), " - ", maximum(gammavals))
println("kappas: ", minimum(kappas), " - ", length(kappas), " - ", maximum(kappas))
println("xs: ", minimum(xs), " - ", length(xs), " - ", maximum(xs))
println("ys: ", minimum(ys), " - ", length(ys), " - ", maximum(ys))
println("Es: ", minimum(Es), " - ", length(Es), " - ", maximum(Es))
println("Lx: ", Lx)
println("Ly: ", Ly)
println("perturbation_type: ", perturbation_type)
println("--------------------------------------------\n")
println("n combs: $(length(Avals) * length(Bvals) * length(mvals) * length(gammavals) * length(xs) * length(ys) * length(Es) * length(kappas)) for $(Lx)x$(Ly) system size\n")
println("--------------------------------------------\n")

results_df = slow_compute_perturbed_spectral_localiser_df(
    perturbation_type=perturbation_type,
    Lx, 
    Ly; 
    Avals=Avals, 
    Bvals=Bvals, 
    mvals=mvals, 
    gammas=gammavals,
    xs=xs, 
    ys=ys, 
    Es=Es, 
    kappas=kappas, 
    periodic_x=false, 
    periodic_y=false, 
    sparse_output=true
)

println("Computed spectral localiser DataFrame with $(nrow(results_df)) rows and $(ncol(results_df)) columns.")

foldername = "data/qwz_model"
isdir(foldername) || mkpath(foldername)
filename = joinpath(foldername, "slow_perturbed_$(perturbation_type)_localiser_scan_As$(Avals[1])-$(length(Avals))-$(Avals[end])_Bs$(Bvals[1])-$(Bvals[end])_ms$(mvals[1])-$(mvals[end])_gammas$(gammavals[1])-$(gammavals[end]).jld2")

println("Saving DataFrame to $filename ...")
t_start = time()
@save filename results_df
println("Saved DataFrame to $filename, took $(round(time() - t_start, digits=2)) seconds to save.")