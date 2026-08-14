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
    perturbation_type::Symbol=:none
)::Matrix{ComplexF64}

    # Base QWZ Hamiltonian in k-space
    H_k = (m + 2B - B * (cos(kx) + cos(ky))) * sigma_z +
          A * sin(kx) * sigma_x +
          A * sin(ky) * sigma_y

    # Apply perturbations
    if perturbation_type == :symmetric
        H_k += gamma * (cos(kx) + cos(ky)) * identity
    elseif perturbation_type == :tilt
        H_k += gamma * sin(kx) * identity
    elseif perturbation_type == :none
        # Do nothing
    else
        error("Unknown perturbation type: $perturbation_type")
    end

    return ComplexF64.(H_k)
end

# 1. Total Density of States (DOS)
function compute_dos(
    eigs::Vector{Float64}, 
    E_range::AbstractVector{<:Real}; 
    eta=0.05
)
    dos = zeros(Float64, length(E_range))
    gaussian(E, E0, sigma) = (1 / (sigma * sqrt(2π))) * exp(-0.5 * ((E - E0) / sigma)^2)
    
    for (i, E_target) in enumerate(E_range)
        dos[i] = sum(gaussian.(eigs, E_target, eta))
    end
    return dos
end

# 2. Local Density of States (LDOS)
function compute_ldos(
    F::Eigen, 
    Lx::Int, 
    Ly::Int, 
    target_E::Real; 
    eta::Real=0.05
)
    ldos = zeros(Float64, Lx, Ly)
    eigs = real(F.values)
    
    # Gaussian broadening function
    gaussian(E, E0, sigma) = (1 / (sigma * sqrt(2π))) * exp(-0.5 * ((E - E0) / sigma)^2)
    
    for i in 1:length(eigs)
        weight = gaussian(eigs[i], target_E, eta)
        
        # Optimization: only sum states that are actually close to target_E
        if weight > 1e-5 
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



# Example usage to plot:
A = 1.0
B = 1.0
m = -1.0
Lxs = [14]
Lys = [14]
gammas = collect(range(0.0, 3.0, 30)) #[0.0, 0.5, 1.0, 2.0]
perturbation_type = :symmetric

periodic_x = false
periodic_y = false

# # OBC spectrum and min|E| eigenstate density
# @showprogress for (Lx, Ly) in zip(Lxs, Lys)
#     for gamma in gammas
#         H_obc = real_space_perturbed_hamiltonian_qwz(
#             Lx, Ly; 
#             A=A, 
#             B=B, 
#             m=m, 
#             gamma=gamma, 
#             perturbation_type=perturbation_type, 
#             periodic_x=periodic_x, 
#             periodic_y=periodic_y, 
#             sparse_output=false
#         )

#         F = eigen(H_obc)
#         eigs = real(F.values)
#         order = sortperm(abs.(eigs))
#         eigs_sorted = eigs[sortperm(eigs)]

#         edge_idx = order[1]
#         psi = F.vectors[:, edge_idx]
#         edge_density = zeros(Float64, Lx, Ly)
#         for y in 1:Ly, x in 1:Lx
#             base = 2 * ((y - 1) * Lx + (x - 1))
#             edge_density[x, y] = abs2(psi[base + 1]) + abs2(psi[base + 2])
#         end

#         p1 = scatter(
#             1:length(eigs_sorted), eigs_sorted;
#             label="OBC eigenvalues",
#             markerstrokewidth=0,
#             markersize=3,
#             color=:slateblue,
#             xlabel="sorted state index",
#             ylabel="energy",
#             title="Open-boundary spectrum"
#         )
#         hline!(p1, [0.0]; linestyle=:dash, color=:black, label=false)

#         p2 = heatmap(
#             1:Lx, 1:Ly, edge_density';
#             xlabel="x", ylabel="y",
#             title="Lowest-|E| eigenstate density",
#             color=:viridis,
#             aspect_ratio=1,
#             colorbar_title="|ψ|²"
#         )

#         plot(p1, p2; layout=(1, 2), size=(800, 600))

#         foldername = joinpath("plots", "pert_obc_spec", "m$(m)")
#         isdir(foldername) || mkpath(foldername)

#         savefig(joinpath(foldername, "qwz_obc_edge_state_A$(A)_B$(B)_m$(m)_gamma$(gamma)_perturbation_$(perturbation_type)_Lx$(Lx)_Ly$(Ly).png"))
#     end
# end


# Calculate LDOS at E = 0.0 with a broadening of 0.1
target_E = 2.0
broadening_width = 0.1

# @showprogress for (Lx, Ly) in zip(Lxs, Lys)
#     for gamma in gammas
#         H_obc = real_space_perturbed_hamiltonian_qwz(
#             Lx, Ly; 
#             A=A, 
#             B=B, 
#             m=m, 
#             gamma=gamma, 
#             perturbation_type=perturbation_type, 
#             periodic_x=periodic_x, 
#             periodic_y=periodic_y, 
#             sparse_output=false
#         )

#         F = eigen(H_obc)

#         ldos_heatmap = compute_ldos(
#             F, 
#             Lx, 
#             Ly, 
#             target_E, 
#             eta=broadening_width
#         )

#         p3 = heatmap(
#             1:Lx, 1:Ly, ldos_heatmap';
#             xlabel="x", ylabel="y",
#             title="LDOS at E = $(target_E), eta = $(broadening_width)",
#             color=:viridis,
#             aspect_ratio=1,
#             colorbar_title="ρ(x,y,E)"
#         )

#         plot(p3; size=(400, 400))

#         foldername = joinpath("plots", "pert_obc_spec", "m$(m)")
#         isdir(foldername) || mkpath(foldername)

#         savefig(joinpath(foldername, "qwz_obc_ldos_A$(A)_B$(B)_m$(m)_gamma$(gamma)_perturbation_$(perturbation_type)_Lx$(Lx)_Ly$(Ly)_E$(target_E)_eta$(broadening_width).png"))
#     end
# end


energy_range = collect(range(-5.0, 5.0, length=101))
Lx = 14
Ly = 14

# df to save results of DOS calc
df = DataFrame(gamma=Float64[], E=Vector{Float64}[], dos=Vector{Float64}[])

@showprogress for gamma in gammas
    H_obc = real_space_perturbed_hamiltonian_qwz(
        Lx, Ly; 
        A=A, 
        B=B, 
        m=m, 
        gamma=gamma, 
        perturbation_type=perturbation_type, 
        periodic_x=periodic_x, 
        periodic_y=periodic_y, 
        sparse_output=false
    )

    F = eigen(H_obc)
    dos = compute_dos(real(F.values), energy_range; eta=0.1)

    push!(df, (gamma=gamma, E=energy_range, dos=dos))
    
    plt = plot(energy_range, dos;
        xlabel="energy",
        ylabel="DOS",
        title="Density of States vs Energy for gamma = $(gamma)",
        legend=false
    )

    foldername = joinpath("plots", "pert_dos")
    isdir(foldername) || mkpath(foldername)
    savefig(plt, joinpath(foldername, "qwz_dos_A$(A)_B$(B)_m$(m)_gamma$(gamma)_perturbation_$(perturbation_type)_Lx$(Lx)_Ly$(Ly).png"))
end

function heatmap_gammaVsE_dos(
    df::DataFrame;
    filename::String="plots/pert_dos/heatmap_gammaVsE_dos.png",
    logscale::Bool=false
)
    
    # Extract unique gamma values and ensure they're sorted
    gammas_unique = sort(unique(df.gamma))
    
    # Get energy range from the first row (assuming all rows have same E values)
    E_range = df.E[1]
    
    # Create matrix: rows = E values, columns = gamma values
    dos_mat = fill(NaN, length(E_range), length(gammas_unique))
    
    for row in eachrow(df)
        gamma_idx = findfirst(==(row.gamma), gammas_unique)
        dos_mat[:, gamma_idx] = row.dos
    end
    
    if logscale
        # Filter out non-positive values before taking log to avoid -Inf/NaN
        # add small value to avoid log(0)
        dos_mat[dos_mat .<= 0] .= NaN
        dos_mat = dos_mat .+ 1e-6  # Add a small value to avoid log(0)
        dos_mat = log10.(dos_mat)
        colorbar_title = "log10(DOS)"
    else
        colorbar_title = "DOS"
    end
    # dos_mat_log = fill(NaN, size(dos_mat))
    # for i in eachindex(dos_mat)
    #     if dos_mat[i] > 0
    #         dos_mat_log[i] = log10(dos_mat[i])
    #     end
    # end
    # dos_mat = dos_mat_log
    # colorbar_title = "log10(DOS)"
    # else
    #     colorbar_title = "DOS"
    # end
    
    p = heatmap(gammas_unique, E_range, dos_mat;
        xlabel="gamma",
        ylabel="Energy",
        title="Density of States: gamma vs Energy",
        colorbar_title=colorbar_title,
        aspect_ratio=:auto)
    
    isdir(dirname(filename)) || mkpath(dirname(filename))
    savefig(p, filename)
    
end

# Generate the heatmap
heatmap_gammaVsE_dos(df; 
    filename=joinpath("plots", "pert_dos", "heatmap_gammaVsE_dos_A$(A)_B$(B)_m$(m)_perturbation_$(perturbation_type)_Lx$(Lx)_Ly$(Ly).png"),
    logscale=true
)