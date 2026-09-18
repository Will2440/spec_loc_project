using Plots
using LaTeXStrings
using DataFrames
using JLD2
using Glob
using ProgressMeter
using Peaks
using Statistics


# specloc_data_path = "/Users/Will/Documents/spec_loc_project/notebook_calcs/data/test_pm_gamma_m-1.0"
specloc_data_path = "/Users/Will/Documents/spec_loc_project/notebook_calcs/data/qwZ_fast_low_spec/clean_low_specloc_spectrum_nEvals10_symmetric_disord_none_PARAMS_As1.0-1-1.0_Bs1.0-1.0_ms-1.0--1.0_gammas0.0-2.0_Ws0.0-0.0_Lxs50-50_Lys50-50_keptSquaretrue_x0y0sweepfalse_Es0.0-3.0_kappas0.001-1.0"




function process_jld2_files_fluid(
    folder_path::String;
    dir_type::Symbol=:folder # :folder or :file
)
    if dir_type == :folder
        file_paths = glob("*.jld2", folder_path)
        
        # Pre-allocate a vector to hold all the dataframes
        dfs = DataFrame[]
        failed_files = String[]

        @showprogress for file_path in file_paths
            jld_data = nothing
            try
                jld_data = load(file_path)
            catch err
                push!(failed_files, file_path)
                @warn "Failed to load JLD2 file (skipping)" file=file_path error=err
                continue
            end

            if haskey(jld_data, "results_df")
                push!(dfs, DataFrame(jld_data["results_df"]))
            else
                @warn "JLD2 file missing results_df key (skipping)" file=file_path
            end
        end

        println("Number of files read: $(length(dfs)) / $(length(file_paths))")
        if !isempty(failed_files)
            @warn "Some JLD2 files could not be processed" failed_count=length(failed_files) failed_files=failed_files
        end
        
        # Combine everything dynamically
        combined_dataframe = isempty(dfs) ? DataFrame() : vcat(dfs..., cols=:union)
    
    elseif dir_type == :file
        # If a single file is provided, load it directly
        jld_data = load(folder_path)
        if haskey(jld_data, "results_df")
            combined_dataframe = DataFrame(jld_data["results_df"])
        else
            error("JLD2 file missing results_df key")
        end
    else
        error("Invalid dir_type. Use :folder or :file.")
    end

    # Print the resulting column names dynamically
    println("\n=== DataFrame Loading Complete ===")
    println("DataFrame column names: | Type | min-length-max")
    for col in names(combined_dataframe)
        print("  - :$col", " (", eltype(combined_dataframe[!, col]), ")", "")
        if col in ["A", "B", "m", "gamma", "kappa", "E", "Lx", "Ly", "x", "y"]
            println(" ", minimum(combined_dataframe[!, col]), " - ", length(unique(combined_dataframe[!, col])), " - ", maximum(combined_dataframe[!, col]))
        end
    end
    println("==================================\n")

    return combined_dataframe
end

specloc_df = process_jld2_files_fluid(specloc_data_path)

specloc_Avals = sort(unique(specloc_df.A))
specloc_Bvals = sort(unique(specloc_df.B))
specloc_mvals = sort(unique(specloc_df.m))
specloc_gammavals = sort(unique(specloc_df.gamma))
specloc_kappavals = sort(unique(specloc_df.kappa))
specloc_Evals = sort(unique(specloc_df.E))

# Keep kappa filtering lazy: pass `kappa = target_kappa` into plotting calls.
target_kappa = 0.03162277660168379
# target_kappa = nothing


######################################################################################################
######################################################################################################

crit_chern_data_path = "/Users/Will/Documents/spec_loc_project/notebook_calcs/data/bulk_band_berry_data/bandstrucutre_data"
target_perturbation = :symmetric

file_paths = glob("*.jld2", crit_chern_data_path)

# Store full tuple: (A, B, m, gamma, band_data)
crit_chern_data_list = Tuple{Float64, Float64, Float64, Float64, Any, Any}[]


for filepath in file_paths
    fname = basename(filepath)

    perturbation_match = match(r"_perturbation_([A-Za-z0-9_]+)\.jld2$", fname)
    if target_perturbation !== nothing
        if isnothing(perturbation_match) || Symbol(perturbation_match.captures[1]) != target_perturbation
            continue
        end
    end
    
    # Parse parameters directly from filename
    m_match = match(r"A([0-9\.\-]+)_B([0-9\.\-]+)_m([0-9\.\-]+)_gamma([0-9\.\-]+)", fname)
    if isnothing(m_match)
        @warn "Filename pattern mismatch, skipping file:" fname
        continue
    end

    A_val     = parse(Float64, m_match.captures[1])
    B_val     = parse(Float64, m_match.captures[2])
    m_val     = parse(Float64, m_match.captures[3])
    gamma_val = parse(Float64, m_match.captures[4])

    jld_data = try
        load(filepath)
    catch err
        @warn "Failed to load JLD2 file" file=filepath error=err
        continue
    end

    if !haskey(jld_data, "data")
        continue
    end
    data = jld_data["data"]

    # Extract band data structure
    band_data = if hasproperty(data, :band)
        data.band
    elseif data isa AbstractDict && haskey(data, :band)
        data[:band]
    else
        data
    end

    # Extract computed euler character
    euler_char_data = if hasproperty(data, :euler_char)
        data.euler_char
    elseif data isa AbstractDict && haskey(data, :euler_char)
        data[:euler_char]
    else
        nothing
    end

    push!(crit_chern_data_list, (A_val, B_val, m_val, gamma_val, band_data, euler_char_data))
end

# Extract unique parameter ranges safely from tuple indices
crit_chern_Avals     = sort(unique([item[1] for item in crit_chern_data_list]))
crit_chern_Bvals     = sort(unique([item[2] for item in crit_chern_data_list]))
crit_chern_mvals     = sort(unique([item[3] for item in crit_chern_data_list]))
crit_chern_gammavals = sort(unique([item[4] for item in crit_chern_data_list]))



######################################################################################################
######################################################################################################

## Check data ranges against each other, set valid matching data ranges
function approx_intersect(v1, v2; tol=1e-6)
    return filter(x -> any(y -> isapprox(x, y; atol=tol), v2), v1)
end

valid_Avals     = approx_intersect(specloc_Avals, crit_chern_Avals)
valid_Bvals     = approx_intersect(specloc_Bvals, crit_chern_Bvals)
valid_mvals     = approx_intersect(specloc_mvals, crit_chern_mvals)
# valid_gammavals = approx_intersect(specloc_gammavals, crit_chern_gammavals)

all_gammavals = sort(unique(vcat(specloc_gammavals, crit_chern_gammavals)))

if any(isempty, (valid_Avals, valid_Bvals, valid_mvals))#, valid_gammavals))
    error("Parameter Mismatch: No overlapping ranges found across (A, B, m).")
end

# Filter specloc DataFrame ONLY by (A, B, m)
specloc_valid_df = filter(row -> 
    any(a -> isapprox(row.A, a; atol=1e-6), valid_Avals) &&
    any(b -> isapprox(row.B, b; atol=1e-6), valid_Bvals) &&
    any(m -> isapprox(row.m, m; atol=1e-6), valid_mvals),
    specloc_df
)

# Filter critical Chern data list ONLY by (A, B, m)
crit_chern_valid_list = filter(item -> 
    any(a -> isapprox(item[1], a; atol=1e-6), valid_Avals) &&
    any(b -> isapprox(item[2], b; atol=1e-6), valid_Bvals) &&
    any(m -> isapprox(item[3], m; atol=1e-6), valid_mvals),
    crit_chern_data_list
)

## Generate plot output path from valid parameter ranges

function fmt_param(vals)
    return length(vals) == 1 ? "$(vals[1])" : "$(minimum(vals))_to_$(maximum(vals))"
end

folder_name = "combined_A$(fmt_param(valid_Avals))_B$(fmt_param(valid_Bvals))_m$(fmt_param(valid_mvals))_gamma$(fmt_param(all_gammavals))"

base_results_path = "/Users/Will/Documents/spec_loc_project/notebook_calcs/results"
output_path = joinpath(base_results_path, folder_name)

isdir(output_path) || mkpath(output_path)
println("Plot output directory set to:\n$output_path")


function filter_df_by_fixed_variables(df::DataFrame, fixed_variables; atol::Real=1e-8)
    nrow(df) == 0 && return df

    col_names = Set(string.(names(df)))
    mask = trues(nrow(df))

    for (k, v) in fixed_variables
        string(k) in col_names || error("Filter key $(k) is not a DataFrame column.")
        col = df[!, k]

        if v isa Real && eltype(col) <: Real
            mask .&= abs.(Float64.(col) .- Float64(v)) .<= Float64(atol)
        else
            mask .&= col .== v
        end
    end

    return df[mask, :]
end

function compute_euler_characteristic(mask::AbstractMatrix{Bool}; periodic::Bool=true)
    nx, ny = size(mask)
    vertices = count(mask)
    edges_h, edges_v, faces = 0, 0, 0

    for j in 1:ny, i in 1:nx
        if mask[i, j]
            i_next = (periodic && i == nx) ? 1 : i + 1
            j_next = (periodic && j == ny) ? 1 : j + 1

            if j_next <= ny && mask[i, j_next]
                edges_h += 1
            end
            if i_next <= nx && mask[i_next, j]
                edges_v += 1
            end
            if i_next <= nx && j_next <= ny && mask[i_next, j] && mask[i, j_next] && mask[i_next, j_next]
                faces += 1
            end
        end
    end

    return vertices - (edges_h + edges_v) + faces
end

function extract_euler_payload(euler_char_data)
    if isnothing(euler_char_data) || euler_char_data isa Missing
        return nothing
    end

    e_mesh = nothing
    chi_vec = nothing

    if hasproperty(euler_char_data, :E_mesh)
        e_mesh = euler_char_data.E_mesh
    elseif euler_char_data isa AbstractDict && haskey(euler_char_data, :E_mesh)
        e_mesh = euler_char_data[:E_mesh]
    elseif euler_char_data isa AbstractDict && haskey(euler_char_data, "E_mesh")
        e_mesh = euler_char_data["E_mesh"]
    end

    if hasproperty(euler_char_data, :euler_char)
        chi_vec = euler_char_data.euler_char
    elseif euler_char_data isa AbstractDict && haskey(euler_char_data, :euler_char)
        chi_vec = euler_char_data[:euler_char]
    elseif euler_char_data isa AbstractDict && haskey(euler_char_data, "euler_char")
        chi_vec = euler_char_data["euler_char"]
    end

    if isnothing(e_mesh) || isnothing(chi_vec)
        return nothing
    end

    e_vals = Float64.(collect(e_mesh))
    chi_vals = Float64.(collect(chi_vec))
    return length(e_vals) == length(chi_vals) ? (e_vals, chi_vals) : nothing
end

# Count the number of connected components (pockets) in a boolean mask
function count_connected_components(mask::AbstractMatrix{Bool}; periodic::Bool=true)
    nx, ny = size(mask)
    visited = falses(nx, ny)
    num_components = 0
    
    # BFS helper to mark connected component
    function bfs!(start_i::Int, start_j::Int)
        queue = [(start_i, start_j)]
        visited[start_i, start_j] = true
        
        while !isempty(queue)
            (i, j) = popfirst!(queue)
            
            # Check 4-neighbors (with periodic wrapping)
            for (di, dj) in [(0, 1), (0, -1), (1, 0), (-1, 0)]
                ni = i + di
                nj = j + dj
                
                if periodic
                    # Use mod() for true modulo wrapping (% is remainder with sign of dividend)
                    ni = 1 + mod(ni - 1, nx)
                    nj = 1 + mod(nj - 1, ny)
                else
                    (1 <= ni <= nx && 1 <= nj <= ny) || continue
                end
                
                if !visited[ni, nj] && mask[ni, nj]
                    visited[ni, nj] = true
                    push!(queue, (ni, nj))
                end
            end
        end
    end
    
    for j in 1:ny, i in 1:nx
        if mask[i, j] && !visited[i, j]
            bfs!(i, j)
            num_components += 1
        end
    end
    
    return num_components
end

# Compute the number of holes and pockets from a boolean occupation mask
# pockets = number of connected components of occupied region
# The relation is: holes + pockets = euler_characteristic
function compute_holes_and_pockets(mask::AbstractMatrix{Bool}; periodic::Bool=true)
    pockets = count_connected_components(mask; periodic=periodic)
    euler_char = compute_euler_characteristic(mask; periodic=periodic)
    holes = euler_char - pockets
    return (pockets, holes)
end

# Create payload for pockets or holes: (e_vals, quantity_vals)
function compute_holes_pockets_payload_from_band_data(
    band_data; 
    quantity::Symbol=:pockets,
    E_mesh=range(0.0, 3.0, length=201), 
    periodic::Bool=true
)
    if !hasproperty(band_data, :energies)
        return nothing
    end

    # Consider both bands: a state is occupied if occupied in either band 1 or band 2
    band1_energies = band_data.energies[1, :, :]
    band2_energies = band_data.energies[2, :, :]
    e_vals = Float64.(collect(E_mesh))
    quantity_vals = Vector{Float64}(undef, length(e_vals))

    for (idx, ef) in enumerate(e_vals)
        # Union of occupied states from both bands
        occupied_mask = (band1_energies .<= ef) .| (band2_energies .<= ef)
        if quantity == :pockets
            quantity_vals[idx] = Float64(count_connected_components(occupied_mask; periodic=periodic))
        elseif quantity == :holes
            pockets, holes = compute_holes_and_pockets(occupied_mask; periodic=periodic)
            quantity_vals[idx] = Float64(holes)
        else
            error("Unknown quantity: $quantity. Use :pockets or :holes.")
        end
    end

    return (e_vals, quantity_vals)
end

function compute_euler_payload_from_band_data(band_data; E_mesh=range(0.0, 3.0, length=201), periodic::Bool=true)
    if !hasproperty(band_data, :energies)
        return nothing
    end

    # Consider both bands: a state is occupied if occupied in either band 1 or band 2
    band1_energies = band_data.energies[1, :, :]
    band2_energies = band_data.energies[2, :, :]
    e_vals = Float64.(collect(E_mesh))
    chi_vals = Vector{Float64}(undef, length(e_vals))

    for (idx, ef) in enumerate(e_vals)
        # Union of occupied states from both bands
        occupied_mask = (band1_energies .<= ef) .| (band2_energies .<= ef)
        chi_vals[idx] = compute_euler_characteristic(occupied_mask; periodic=periodic)
    end

    return (e_vals, chi_vals)
end

function extract_pockets_payload(euler_char_data)
    if isnothing(euler_char_data) || euler_char_data isa Missing
        return nothing
    end

    e_mesh = nothing
    pockets_vec = nothing

    if hasproperty(euler_char_data, :E_mesh)
        e_mesh = euler_char_data.E_mesh
    elseif euler_char_data isa AbstractDict && haskey(euler_char_data, :E_mesh)
        e_mesh = euler_char_data[:E_mesh]
    elseif euler_char_data isa AbstractDict && haskey(euler_char_data, "E_mesh")
        e_mesh = euler_char_data["E_mesh"]
    end

    if hasproperty(euler_char_data, :n_pockets)
        pockets_vec = euler_char_data.n_pockets
    elseif euler_char_data isa AbstractDict && haskey(euler_char_data, :n_pockets)
        pockets_vec = euler_char_data[:n_pockets]
    elseif euler_char_data isa AbstractDict && haskey(euler_char_data, "n_pockets")
        pockets_vec = euler_char_data["n_pockets"]
    end

    if isnothing(e_mesh) || isnothing(pockets_vec)
        return nothing
    end

    e_vals = Float64.(collect(e_mesh))
    pockets_vals = Float64.(collect(pockets_vec))
    return length(e_vals) == length(pockets_vals) ? (e_vals, pockets_vals) : nothing
end

function extract_holes_payload(euler_char_data)
    if isnothing(euler_char_data) || euler_char_data isa Missing
        return nothing
    end

    e_mesh = nothing
    holes_vec = nothing

    if hasproperty(euler_char_data, :E_mesh)
        e_mesh = euler_char_data.E_mesh
    elseif euler_char_data isa AbstractDict && haskey(euler_char_data, :E_mesh)
        e_mesh = euler_char_data[:E_mesh]
    elseif euler_char_data isa AbstractDict && haskey(euler_char_data, "E_mesh")
        e_mesh = euler_char_data["E_mesh"]
    end

    if hasproperty(euler_char_data, :n_holes)
        holes_vec = euler_char_data.n_holes
    elseif euler_char_data isa AbstractDict && haskey(euler_char_data, :n_holes)
        holes_vec = euler_char_data[:n_holes]
    elseif euler_char_data isa AbstractDict && haskey(euler_char_data, "n_holes")
        holes_vec = euler_char_data["n_holes"]
    end

    if isnothing(e_mesh) || isnothing(holes_vec)
        return nothing
    end

    e_vals = Float64.(collect(e_mesh))
    holes_vals = Float64.(collect(holes_vec))
    return length(e_vals) == length(holes_vals) ? (e_vals, holes_vals) : nothing
end



######################################################################################################
######################################################################################################


function plt_specloc_gamma_vs_E_overlaid_crit_cum_chern(
    df::DataFrame,
    gamma_data_pairs::Vector{<:Tuple{Real, Any}}; 
    atol::Real=1e-8,
    logscale::Bool=false,
    target_p::Float64=0.5,
    target_q::Float64=0.0,
    target_r::Float64=1.0,
    tol::Float64=0.01,
    filename::String="plots/localiser_gamma_E_heatmap.png",
    xlims::Union{Nothing, Tuple{Real, Real}}=nothing,
    ylims::Union{Nothing, Tuple{Real, Real}}=nothing,
    fixed_variables...
)
    # --------------------------------------------------------------------------
    # 1. Filter DataFrame by fixed parameter values
    # --------------------------------------------------------------------------
    col_names_str = string.(names(df))
    subdf = df
    for (k, v) in fixed_variables
        string(k) in col_names_str || error("Filter key $(k) is not a DataFrame column.")
        col = subdf[!, k]
        if v isa Real && eltype(col) <: Real
            mask = abs.(Float64.(col) .- Float64(v)) .<= Float64(atol)
            subdf = subdf[mask, :]
        else
            mask = col .== v
            subdf = subdf[mask, :]
        end
    end

    # --------------------------------------------------------------------------
    # 2. Extract Heatmap Data (Localiser Gap & Chern)
    # --------------------------------------------------------------------------
    gammas = sort(unique(subdf.gamma))
    Es = sort(unique(subdf.E))
    ngamma, nE = length(gammas), length(Es)

    gap_mat = fill(NaN, nE, ngamma)
    chern_mat = fill(NaN, nE, ngamma)

    for row in eachrow(subdf)
        xi = findfirst(row.gamma .== gammas)
        yi = findfirst(row.E .== Es)
        if !isnothing(xi) && !isnothing(yi)
            gap_mat[yi, xi] = row.localiser_gap
            chern_mat[yi, xi] = ismissing(row.chern_number) ? NaN : Float64(row.chern_number)
        end
    end

    if logscale
        gap_mat = log10.(gap_mat)
        gap_title = "log10(min|λ|)"
    else
        gap_title = "min|λ|"
    end

    # --------------------------------------------------------------------------
    # 3. Extract Scatter Contour Points from Critical Chern Data
    # --------------------------------------------------------------------------
    g_b1_pos_p, e_b1_pos_p = Float64[], Float64[]
    g_b1_neg_p, e_b1_neg_p = Float64[], Float64[]
    g_b1_pos_q, e_b1_pos_q = Float64[], Float64[]
    g_b1_neg_q, e_b1_neg_q = Float64[], Float64[]

    g_b2_pos_p, e_b2_pos_p = Float64[], Float64[]
    g_b2_neg_p, e_b2_neg_p = Float64[], Float64[]
    g_b2_pos_q, e_b2_pos_q = Float64[], Float64[]
    g_b2_neg_q, e_b2_neg_q = Float64[], Float64[]

    q_is_zero = abs(target_q) < 1e-8

    for (gamma, band_data) in gamma_data_pairs
        g = Float64(gamma)

        cum_1 = band_data.cum_chern_per_band[1, :, :]
        E_1   = band_data.plaquette_energies[1, :, :]
        cum_2 = band_data.cum_chern_per_band[2, :, :]
        E_2   = band_data.plaquette_energies[2, :, :]

        for i in findall(v -> abs(v - target_p) <= tol, cum_1); push!(g_b1_pos_p, g); push!(e_b1_pos_p, E_1[i]); end
        for i in findall(v -> abs(v + target_p) <= tol, cum_1); push!(g_b1_neg_p, g); push!(e_b1_neg_p, E_1[i]); end
        for i in findall(v -> abs(v - target_q) <= tol, cum_1); push!(g_b1_pos_q, g); push!(e_b1_pos_q, E_1[i]); end
        if !q_is_zero
            for i in findall(v -> abs(v + target_q) <= tol, cum_1); push!(g_b1_neg_q, g); push!(e_b1_neg_q, E_1[i]); end
        end

        for i in findall(v -> abs(v - target_p) <= tol, cum_2); push!(g_b2_pos_p, g); push!(e_b2_pos_p, E_2[i]); end
        for i in findall(v -> abs(v + target_p) <= tol, cum_2); push!(g_b2_neg_p, g); push!(e_b2_neg_p, E_2[i]); end
        for i in findall(v -> abs(v - target_q) <= tol, cum_2); push!(g_b2_pos_q, g); push!(e_b2_pos_q, E_2[i]); end
        if !q_is_zero
            for i in findall(v -> abs(v + target_q) <= tol, cum_2); push!(g_b2_neg_q, g); push!(e_b2_neg_q, E_2[i]); end
        end
    end

    lbl_q_b1_pos = q_is_zero ? "Band 1, " * L" %$(target_q)" : "Band 1, " * L" +%$(target_q)"
    lbl_q_b2_pos = q_is_zero ? "Band 2, " * L" %$(target_q)" : "Band 2, " * L" +%$(target_q)"

    series_configs = [
        (g_b1_pos_p, e_b1_pos_p, "Band 1, " * L" +%$(target_p)", :cyan,        :circle),
        (g_b1_neg_p, e_b1_neg_p, "Band 1, " * L" -%$(target_p)", :dodgerblue,  :diamond),
        (g_b1_pos_q, e_b1_pos_q, lbl_q_b1_pos,                  :springgreen, :utriangle),
        (g_b1_neg_q, e_b1_neg_q, "Band 1, " * L" -%$(target_q)", :lime,        :dtriangle),
        (g_b2_pos_p, e_b2_pos_p, "Band 2, " * L" +%$(target_p)", :magenta,     :circle),
        (g_b2_neg_p, e_b2_neg_p, "Band 2, " * L" -%$(target_p)", :hotpink,     :diamond),
        (g_b2_pos_q, e_b2_pos_q, lbl_q_b2_pos,                  :yellow,      :utriangle),
        (g_b2_neg_q, e_b2_neg_q, "Band 2, " * L" -%$(target_q)", :orange,      :dtriangle)
    ]

    # --------------------------------------------------------------------------
    # 4. Generate Heatmaps & Overlay Scatter Points (Equal & Square Layout)
    # --------------------------------------------------------------------------
    fixed_title_str = join(["$k=$v" for (k, v) in fixed_variables], ", ")

    # Explicitly enable colorbars on BOTH subplots for symmetrical layout width
    p1 = heatmap(gammas, Es, gap_mat;
        xlabel=L"\gamma", ylabel=L"E",
        title="Localiser Gap - ($(fixed_title_str))",
        colorbar_title=gap_title,
        colorbar=true,
        colormap=:plasma)
        # aspect_ratio=:equal)
    
    p2 = heatmap(gammas, Es, chern_mat;
        xlabel=L"\gamma", ylabel=L"E",
        title="Localiser Chern Number",
        colorbar_title="Chern",
        colorbar=true,
        colormap=:RdBu)
        # aspect_ratio=: equal)

    # Overlay scatter series on both heatmaps
    for (g_pts, e_pts, lbl, clr, shp) in series_configs
        if !isempty(g_pts)
            scatter!(p1, g_pts, e_pts;
                label=lbl,
                color=clr,
                markershape=shp,
                markersize=3.5,
                markerstrokewidth=0.4,
                markerstrokecolor=:black,
                alpha=0.85,
                xlims=xlims,
                ylims=ylims
            )

            scatter!(p2, g_pts, e_pts;
                label=lbl,
                color=clr,
                markershape=shp,
                markersize=3.5,
                markerstrokewidth=0.4,
                markerstrokecolor=:black,
                alpha=0.85,
                xlims=xlims,
                ylims=ylims
            )
        end
    end

    # Set inside legend on p2 to prevent shrinking p2 relative to p1
    plot!(p2, legend=:topright)

    # 1100x550 size provides 2 side-by-side square panels with colorbars
    plt = plot(p1, p2; layout=(1, 2), size=(1600, 600), margin=5Plots.mm)
    
    isdir(dirname(filename)) || mkpath(dirname(filename))
    savefig(plt, filename)

    return plt
end

function plt_specloc_gamma_vs_E_overlaid_crit_cum_chern_three_targets(
    df::DataFrame,
    gamma_data_pairs::Union{Vector{<:Tuple{Real, Any}}, Nothing}=nothing; 
    atol::Real=1e-8,
    logscale::Bool=false,
    target_p::Float64=0.5,
    target_q::Float64=0.0,
    target_r::Float64=1.0,
    tol::Float64=0.01,
    filename::String="plots/localiser_gamma_E_heatmap.png",
    xlims::Union{Nothing, Tuple{Real, Real}}=nothing,
    ylims::Union{Nothing, Tuple{Real, Real}}=nothing,
    fixed_variables...
)
    # --------------------------------------------------------------------------
    # 1. Filter DataFrame by fixed parameter values
    # --------------------------------------------------------------------------
    col_names_str = string.(names(df))
    subdf = df
    for (k, v) in fixed_variables
        string(k) in col_names_str || error("Filter key $(k) is not a DataFrame column.")
        col = subdf[!, k]
        if v isa Real && eltype(col) <: Real
            mask = abs.(Float64.(col) .- Float64(v)) .<= Float64(atol)
            subdf = subdf[mask, :]
        else
            mask = col .== v
            subdf = subdf[mask, :]
        end
    end

    if nrow(subdf) == 0
        @warn "No data found matching parameters: $fixed_variables"
        return nothing
    end

    # --------------------------------------------------------------------------
    # 2. Extract Heatmap Data (Inferred Gap & Chern)
    # --------------------------------------------------------------------------
    gammas = sort(unique(subdf.gamma))
    Es = sort(unique(subdf.E))
    ngamma, nE = length(gammas), length(Es)

    gap_mat = fill(NaN, nE, ngamma)
    chern_mat = fill(NaN, nE, ngamma)

    for row in eachrow(subdf)
        xi = findfirst(==(row.gamma), gammas)
        yi = findfirst(==(row.E), Es)
        
        if !isnothing(xi) && !isnothing(yi)
            # Infer localiser gap from low_lying_evals if localiser_gap is absent/missing
            gap_val = if hasproperty(row, :localiser_gap) && !ismissing(row.localiser_gap)
                Float64(row.localiser_gap)
            elseif hasproperty(row, :low_lying_evals) && !ismissing(row.low_lying_evals)
                minimum(abs.(vec(row.low_lying_evals)))
            else
                NaN
            end

            gap_mat[yi, xi] = gap_val
            chern_mat[yi, xi] = hasproperty(row, :chern_number) && !ismissing(row.chern_number) ? Float64(row.chern_number) : NaN
        end
    end

    # Log-scale transformation logic
    if logscale
        gap_mat = log10.(max.(gap_mat, 1e-12))
        gap_title = L"\log_{10}(\text{min}|\lambda|)"
    else
        gap_title = L"\text{min}|\lambda|"
    end

    # Explicit clims safety bounds
    valid_gaps = filter(!isnan, vec(gap_mat))
    c_min, c_max = isempty(valid_gaps) ? (0.0, 1.0) : (minimum(valid_gaps), maximum(valid_gaps))
    if c_min == c_max
        c_min -= 1e-5
        c_max += 1e-5
    end

    fixed_title_str = join(["$k=$v" for (k, v) in fixed_variables], ", ")

    # --------------------------------------------------------------------------
    # 3. Base Heatmaps
    # --------------------------------------------------------------------------
    p1 = heatmap(gammas, Es, gap_mat;
        clims = (c_min, c_max),
        xlabel = L"\gamma", ylabel = L"E",
        title = "Localiser Gap - ($fixed_title_str)",
        colorbar_title = gap_title,
        colorbar = true,
        colormap = :plasma,
        aspect_ratio = :auto)
    
    p2 = heatmap(gammas, Es, chern_mat;
        xlabel = L"\gamma", ylabel = L"E",
        title = "Localiser Chern Number",
        colorbar_title = "Chern",
        colorbar = true,
        colormap = :RdBu,
        aspect_ratio = :auto)

    # --------------------------------------------------------------------------
    # 4. Extract & Overlay Scatter Contour Points (3 Targets: p, q, r)
    # --------------------------------------------------------------------------
    if !isnothing(gamma_data_pairs)
        g_b1_pos_p, e_b1_pos_p = Float64[], Float64[]
        g_b1_neg_p, e_b1_neg_p = Float64[], Float64[]
        g_b1_pos_q, e_b1_pos_q = Float64[], Float64[]
        g_b1_neg_q, e_b1_neg_q = Float64[], Float64[]
        g_b1_pos_r, e_b1_pos_r = Float64[], Float64[]
        g_b1_neg_r, e_b1_neg_r = Float64[], Float64[]

        g_b2_pos_p, e_b2_pos_p = Float64[], Float64[]
        g_b2_neg_p, e_b2_neg_p = Float64[], Float64[]
        g_b2_pos_q, e_b2_pos_q = Float64[], Float64[]
        g_b2_neg_q, e_b2_neg_q = Float64[], Float64[]
        g_b2_pos_r, e_b2_pos_r = Float64[], Float64[]
        g_b2_neg_r, e_b2_neg_r = Float64[], Float64[]

        q_is_zero = abs(target_q) < 1e-8
        r_is_zero = abs(target_r) < 1e-8

        for (gamma, band_data) in gamma_data_pairs
            g = Float64(gamma)

            cum_1 = band_data.cum_chern_per_band[1, :, :]
            E_1   = band_data.plaquette_energies[1, :, :]
            cum_2 = band_data.cum_chern_per_band[2, :, :]
            E_2   = band_data.plaquette_energies[2, :, :]

            # --- Band 1 ---
            for i in findall(v -> abs(v - target_p) <= tol, cum_1); push!(g_b1_pos_p, g); push!(e_b1_pos_p, E_1[i]); end
            for i in findall(v -> abs(v + target_p) <= tol, cum_1); push!(g_b1_neg_p, g); push!(e_b1_neg_p, E_1[i]); end
            for i in findall(v -> abs(v - target_q) <= tol, cum_1); push!(g_b1_pos_q, g); push!(e_b1_pos_q, E_1[i]); end
            if !q_is_zero
                for i in findall(v -> abs(v + target_q) <= tol, cum_1); push!(g_b1_neg_q, g); push!(e_b1_neg_q, E_1[i]); end
            end
            for i in findall(v -> abs(v - target_r) <= tol, cum_1); push!(g_b1_pos_r, g); push!(e_b1_pos_r, E_1[i]); end
            if !r_is_zero
                for i in findall(v -> abs(v + target_r) <= tol, cum_1); push!(g_b1_neg_r, g); push!(e_b1_neg_r, E_1[i]); end
            end

            # --- Band 2 ---
            for i in findall(v -> abs(v - target_p) <= tol, cum_2); push!(g_b2_pos_p, g); push!(e_b2_pos_p, E_2[i]); end
            for i in findall(v -> abs(v + target_p) <= tol, cum_2); push!(g_b2_neg_p, g); push!(e_b2_neg_p, E_2[i]); end
            for i in findall(v -> abs(v - target_q) <= tol, cum_2); push!(g_b2_pos_q, g); push!(e_b2_pos_q, E_2[i]); end
            if !q_is_zero
                for i in findall(v -> abs(v + target_q) <= tol, cum_2); push!(g_b2_neg_q, g); push!(e_b2_neg_q, E_2[i]); end
            end
            for i in findall(v -> abs(v - target_r) <= tol, cum_2); push!(g_b2_pos_r, g); push!(e_b2_pos_r, E_2[i]); end
            if !r_is_zero
                for i in findall(v -> abs(v + target_r) <= tol, cum_2); push!(g_b2_neg_r, g); push!(e_b2_neg_r, E_2[i]); end
            end
        end

        lbl_q_b1_pos = q_is_zero ? "Band 1, " * L" %$(target_q)" : "Band 1, " * L" +%$(target_q)"
        lbl_q_b2_pos = q_is_zero ? "Band 2, " * L" %$(target_q)" : "Band 2, " * L" +%$(target_q)"

        lbl_r_b1_pos = r_is_zero ? "Band 1, " * L" %$(target_r)" : "Band 1, " * L" +%$(target_r)"
        lbl_r_b2_pos = r_is_zero ? "Band 2, " * L" %$(target_r)" : "Band 2, " * L" +%$(target_r)"

        series_configs = [
            # Band 1 series
            (g_b1_pos_p, e_b1_pos_p, "Band 1, " * L" +%$(target_p)", :cyan,         :circle),
            (g_b1_neg_p, e_b1_neg_p, "Band 1, " * L" -%$(target_p)", :dodgerblue,   :diamond),
            (g_b1_pos_q, e_b1_pos_q, lbl_q_b1_pos,                  :springgreen,  :utriangle),
            (g_b1_neg_q, e_b1_neg_q, "Band 1, " * L" -%$(target_q)", :lime,         :dtriangle),
            (g_b1_pos_r, e_b1_pos_r, lbl_r_b1_pos,                  :mediumpurple, :star5),
            (g_b1_neg_r, e_b1_neg_r, "Band 1, " * L" -%$(target_r)", :purple,       :pentagon),

            # Band 2 series
            (g_b2_pos_p, e_b2_pos_p, "Band 2, " * L" +%$(target_p)", :magenta,      :circle),
            (g_b2_neg_p, e_b2_neg_p, "Band 2, " * L" -%$(target_p)", :hotpink,      :diamond),
            (g_b2_pos_q, e_b2_pos_q, lbl_q_b2_pos,                  :yellow,       :utriangle),
            (g_b2_neg_q, e_b2_neg_q, "Band 2, " * L" -%$(target_q)", :orange,       :dtriangle),
            (g_b2_pos_r, e_b2_pos_r, lbl_r_b2_pos,                  :crimson,      :star5),
            (g_b2_neg_r, e_b2_neg_r, "Band 2, " * L" -%$(target_r)", :darkred,      :pentagon)
        ]

        # Overlay scatter series on both heatmaps
        for (g_pts, e_pts, lbl, clr, shp) in series_configs
            if !isempty(g_pts)
                scatter!(p1, g_pts, e_pts;
                    label=lbl,
                    color=clr,
                    markershape=shp,
                    markersize=3.5,
                    markerstrokewidth=0.4,
                    markerstrokecolor=:black,
                    alpha=0.85,
                    xlims=xlims,
                    ylims=ylims
                )

                scatter!(p2, g_pts, e_pts;
                    label=lbl,
                    color=clr,
                    markershape=shp,
                    markersize=3.5,
                    markerstrokewidth=0.4,
                    markerstrokecolor=:black,
                    alpha=0.85,
                    xlims=xlims,
                    ylims=ylims
                )
            end
        end
    end

    plot!(p2, legend=:topright)

    plt = plot(p1, p2; layout=(1, 2), size=(1600, 600), margin=5Plots.mm)
    
    if !isempty(filename)
        isdir(dirname(filename)) || mkpath(dirname(filename))
        savefig(plt, filename)
    end

    return plt
end




function plt_specloc_gamma_vs_E_overlaid_joint_band_crit_cum_chern(
    df::DataFrame,
    gamma_data_pairs::Vector{<:Tuple{Real, Any}}; 
    atol::Real=1e-8,
    logscale::Bool=false,
    target_p::Float64=0.5,
    target_q::Float64=0.0,
    tol::Float64=0.01,
    filename::String="plots/localiser_gamma_E_joint_heatmap.png",
    xlims::Union{Nothing, Tuple{Real, Real}}=nothing,
    ylims::Union{Nothing, Tuple{Real, Real}}=nothing,
    fixed_variables...
)
    # --------------------------------------------------------------------------
    # 1. Filter DataFrame by fixed parameter values
    # --------------------------------------------------------------------------
    col_names_str = string.(names(df))
    subdf = df
    for (k, v) in fixed_variables
        string(k) in col_names_str || error("Filter key $(k) is not a DataFrame column.")
        col = subdf[!, k]
        if v isa Real && eltype(col) <: Real
            mask = abs.(Float64.(col) .- Float64(v)) .<= Float64(atol)
            subdf = subdf[mask, :]
        else
            mask = col .== v
            subdf = subdf[mask, :]
        end
    end

    # --------------------------------------------------------------------------
    # 2. Extract Heatmap Data (Localiser Gap & Chern)
    # --------------------------------------------------------------------------
    gammas = sort(unique(subdf.gamma))
    Es = sort(unique(subdf.E))
    ngamma, nE = length(gammas), length(Es)

    gap_mat = fill(NaN, nE, ngamma)
    chern_mat = fill(NaN, nE, ngamma)

    for row in eachrow(subdf)
        xi = findfirst(row.gamma .== gammas)
        yi = findfirst(row.E .== Es)
        if !isnothing(xi) && !isnothing(yi)
            gap_mat[yi, xi] = row.localiser_gap
            chern_mat[yi, xi] = ismissing(row.chern_number) ? NaN : Float64(row.chern_number)
        end
    end

    if logscale
        gap_mat = log10.(gap_mat)
        gap_title = "log10(min|λ|)"
    else
        gap_title = "min|λ|"
    end

    # --------------------------------------------------------------------------
    # 3. Extract Scatter Points from Joint Band Accumulated Chern Data
    # --------------------------------------------------------------------------
    g_pos_p, e_pos_p = Float64[], Float64[]
    g_neg_p, e_neg_p = Float64[], Float64[]
    g_pos_q, e_pos_q = Float64[], Float64[]
    g_neg_q, e_neg_q = Float64[], Float64[]

    q_is_zero = abs(target_q) < 1e-8

    for (gamma, band_data) in gamma_data_pairs
        g = Float64(gamma)

        # 1. Flatten energy across all bands
        E_flat = vec(band_data.plaquette_energies)

        # 2. Compute global energy-accumulated Chern number across all overlapping bands
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

        # 3. Filter points matching target contours
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

    lbl_q_pos = q_is_zero ? L"C(E) = %$(target_q)" : L"C(E) = +%$(target_q)"

    series_configs = [
        (g_pos_p, e_pos_p, L"C(E) = +%$(target_p)", :cyan,        :circle),
        (g_neg_p, e_neg_p, L"C(E) = -%$(target_p)", :dodgerblue,  :diamond),
        (g_pos_q, e_pos_q, lbl_q_pos,               :springgreen, :utriangle),
        (g_neg_q, e_neg_q, L"C(E) = -%$(target_q)", :lime,        :dtriangle)
    ]

    # --------------------------------------------------------------------------
    # 4. Generate Heatmaps & Overlay Scatter Points
    # --------------------------------------------------------------------------
    fixed_title_str = join(["$k=$v" for (k, v) in fixed_variables], ", ")

    p1 = heatmap(gammas, Es, gap_mat;
        xlabel=L"\gamma", ylabel=L"E",
        title="Localiser Gap - ($(fixed_title_str))",
        colorbar_title=gap_title,
        colorbar=true,
        colormap=:plasma)
    
    p2 = heatmap(gammas, Es, chern_mat;
        xlabel=L"\gamma", ylabel=L"E",
        title="Localiser Chern Number",
        colorbar_title="Chern",
        colorbar=true,
        colormap=:RdBu)

    # Overlay global joint scatter series on both heatmaps
    for (g_pts, e_pts, lbl, clr, shp) in series_configs
        if !isempty(g_pts)
            scatter!(p1, g_pts, e_pts;
                label=lbl,
                color=clr,
                markershape=shp,
                markersize=3.5,
                markerstrokewidth=0.4,
                markerstrokecolor=:black,
                alpha=0.85,
                xlims=xlims,
                ylims=ylims
            )

            scatter!(p2, g_pts, e_pts;
                label=lbl,
                color=clr,
                markershape=shp,
                markersize=3.5,
                markerstrokewidth=0.4,
                markerstrokecolor=:black,
                alpha=0.85,
                xlims=xlims,
                ylims=ylims
            )
        end
    end

    plot!(p2, legend=:topright)

    plt = plot(p1, p2; layout=(1, 2), size=(1600, 600), margin=5Plots.mm)
    
    isdir(dirname(filename)) || mkpath(dirname(filename))
    savefig(plt, filename)

    return plt
end

function plt_specloc_gamma_vs_E_overlaid_joint_band_max_cum_chern(
    df::DataFrame,
    gamma_data_pairs::Vector{<:Tuple{Real, Any}}; 
    atol::Real=1e-8,
    logscale::Bool=false,
    target_zero::Float64=0.0,
    tol::Float64=0.01,
    filename::String="plots/localiser_gamma_E_max_cum_chern_heatmap.png",
    xlims::Union{Nothing, Tuple{Real, Real}}=nothing,
    ylims::Union{Nothing, Tuple{Real, Real}}=nothing,
    max_colormap::Symbol=:viridis,
    fixed_variables...
)
    # --------------------------------------------------------------------------
    # 1. Filter DataFrame by fixed parameter values
    # --------------------------------------------------------------------------
    col_names_str = string.(names(df))
    subdf = df
    for (k, v) in fixed_variables
        string(k) in col_names_str || error("Filter key $(k) is not a DataFrame column.")
        col = subdf[!, k]
        if v isa Real && eltype(col) <: Real
            mask = abs.(Float64.(col) .- Float64(v)) .<= Float64(atol)
            subdf = subdf[mask, :]
        else
            mask = col .== v
            subdf = subdf[mask, :]
        end
    end

    # --------------------------------------------------------------------------
    # 2. Extract Heatmap Data (Localiser Gap & Chern)
    # --------------------------------------------------------------------------
    gammas = sort(unique(subdf.gamma))
    Es = sort(unique(subdf.E))
    ngamma, nE = length(gammas), length(Es)

    gap_mat = fill(NaN, nE, ngamma)
    chern_mat = fill(NaN, nE, ngamma)

    for row in eachrow(subdf)
        xi = findfirst(row.gamma .== gammas)
        yi = findfirst(row.E .== Es)
        if !isnothing(xi) && !isnothing(yi)
            gap_mat[yi, xi] = row.localiser_gap
            chern_mat[yi, xi] = ismissing(row.chern_number) ? NaN : Float64(row.chern_number)
        end
    end

    if logscale
        gap_mat = log10.(gap_mat)
        gap_title = "log10(min|λ|)"
    else
        gap_title = "min|λ|"
    end

    # --------------------------------------------------------------------------
    # 3. Extract Scatter Data (C = 0 contours & Max Cum Chern per gamma)
    # --------------------------------------------------------------------------
    g_zero, e_zero = Float64[], Float64[]
    g_max, e_max, c_max = Float64[], Float64[], Float64[]

    for (gamma, band_data) in gamma_data_pairs
        g = Float64(gamma)

        # Flatten energy across all bands
        E_flat = vec(band_data.plaquette_energies)

        # Compute global energy-accumulated Chern number
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

        # (a) Track critical C = 0 points
        for i in findall(v -> abs(v - target_zero) <= tol, cum_flat)
            push!(g_zero, g)
            push!(e_zero, E_flat[i])
        end

        # (b) Track Maximum Cumulative Chern value and its energy locations for this gamma
        if !isempty(cum_flat)
            max_val = maximum(cum_flat)
            for i in findall(v -> abs(v - max_val) <= tol, cum_flat)
                push!(g_max, g)
                push!(e_max, E_flat[i])
                push!(c_max, max_val)
            end
        end
    end

    # Compute color limits for maximum cumulative Chern scatter points
    clims_max = if !isempty(c_max)
        cmin, cmax_val = minimum(c_max), maximum(c_max)
        cmin == cmax_val ? (cmin - 0.1, cmax_val + 0.1) : (cmin, cmax_val)
    else
        (0.0, 1.0)
    end

    # --------------------------------------------------------------------------
    # 4. Generate Subplots with Independent Colorbars
    # --------------------------------------------------------------------------
    fixed_title_str = join(["$k=$v" for (k, v) in fixed_variables], ", ")

    p1 = heatmap(gammas, Es, gap_mat;
        xlabel=L"\gamma", ylabel=L"E",
        title="Localiser Gap - ($(fixed_title_str))",
        colorbar_title=gap_title,
        colorbar=true,
        colormap=:plasma)

    p2 = heatmap(gammas, Es, chern_mat;
        xlabel=L"\gamma", ylabel=L"E",
        title="Localiser Chern Number",
        colorbar_title="Chern",
        colorbar=true,
        colormap=:RdBu)

    # Panel 3: Dedicated plot providing an independent color bar scale for Max C(E)
    p3 = scatter(
        xlabel=L"\gamma", ylabel=L"E",
        title=L"Max\ Cumulative\ Chern\ C_{\max}(E)",
        xlims=xlims, ylims=ylims
    )

    # Overlay C = 0 scatter points on all panels
    if !isempty(g_zero)
        for p in (p1, p2, p3)
            scatter!(p, g_zero, e_zero;
                label=L"C = %$(target_zero)",
                color=:white,
                markershape=:diamond,
                markersize=3.5,
                markerstrokewidth=0.6,
                markerstrokecolor=:black,
                alpha=0.9
            )
        end
    end

    # Overlay Max Cum Chern scatter points color-coded by c_max
    if !isempty(g_max)
        # On p1 and p2: overlay colored markers using identical colormap without duplicating colorbars
        scatter!(p1, g_max, e_max;
            marker_z=c_max,
            label=L"C_{\max}",
            colormap=max_colormap,
            clims=clims_max,
            colorbar=false,
            markershape=:circle,
            markersize=4.0,
            markerstrokewidth=0.4,
            markerstrokecolor=:black,
            alpha=0.9
        )

        scatter!(p2, g_max, e_max;
            marker_z=c_max,
            label=L"C_{\max}",
            colormap=max_colormap,
            clims=clims_max,
            colorbar=false,
            markershape=:circle,
            markersize=4.0,
            markerstrokewidth=0.4,
            markerstrokecolor=:black,
            alpha=0.9
        )

        # On p3: render with independent colorbar scale enabled
        scatter!(p3, g_max, e_max;
            marker_z=c_max,
            label=L"C_{\max}",
            colorbar_title=L"C_{\max}",
            colormap=max_colormap,
            clims=clims_max,
            colorbar=true,
            markershape=:circle,
            markersize=4.0,
            markerstrokewidth=0.4,
            markerstrokecolor=:black,
            alpha=0.9
        )
    end

    plot!(p1, legend=:topright, xlims=xlims, ylims=ylims)
    plot!(p2, legend=:topright, xlims=xlims, ylims=ylims)
    plot!(p3, legend=:topright, xlims=xlims, ylims=ylims)

    # 3-panel layout (2400x600 px) to cleanly fit all colorbars side-by-side
    plt = plot(p1, p2, p3; layout=(1, 3), size=(2400, 600), margin=5Plots.mm)
    
    isdir(dirname(filename)) || mkpath(dirname(filename))
    savefig(plt, filename)

    return plt
end

function plt_specloc_gamma_vs_E_overlaid_joint_band_smart_extrema(
    df::DataFrame,
    gamma_data_pairs::Vector{<:Tuple{Real, Any}}; 
    atol::Real=1e-8,
    logscale::Bool=false,
    peak_tol::Float64=0.01,
    e_cluster_tol::Float64=0.1,
    flat_span_tol::Float64=0.02,
    zero_tol::Float64=1e-4,
    filename::String="plots/localiser_gamma_E_smart_extrema_heatmap.png",
    xlims::Union{Nothing, Tuple{Real, Real}}=nothing,
    ylims::Union{Nothing, Tuple{Real, Real}}=nothing,
    fixed_variables...
)
    # --------------------------------------------------------------------------
    # 1. Filter DataFrame by fixed parameter values
    # --------------------------------------------------------------------------
    col_names_str = string.(names(df))
    subdf = df
    for (k, v) in fixed_variables
        string(k) in col_names_str || error("Filter key $(k) is not a DataFrame column.")
        col = subdf[!, k]
        if v isa Real && eltype(col) <: Real
            mask = abs.(Float64.(col) .- Float64(v)) .<= Float64(atol)
            subdf = subdf[mask, :]
        else
            mask = col .== v
            subdf = subdf[mask, :]
        end
    end

    # --------------------------------------------------------------------------
    # 2. Extract Heatmap Data (Localiser Gap & Chern)
    # --------------------------------------------------------------------------
    gammas = sort(unique(subdf.gamma))
    Es = sort(unique(subdf.E))
    ngamma, nE = length(gammas), length(Es)

    gap_mat = fill(NaN, nE, ngamma)
    chern_mat = fill(NaN, nE, ngamma)

    for row in eachrow(subdf)
        xi = findfirst(row.gamma .== gammas)
        yi = findfirst(row.E .== Es)
        if !isnothing(xi) && !isnothing(yi)
            gap_mat[yi, xi] = row.localiser_gap
            chern_mat[yi, xi] = ismissing(row.chern_number) ? NaN : Float64(row.chern_number)
        end
    end

    if logscale
        gap_mat = log10.(gap_mat)
        gap_title = "log10(min|λ|)"
    else
        gap_title = "min|λ|"
    end

    # --------------------------------------------------------------------------
    # 3. Extract Smart Extrema Scatter Points from Joint Band Chern Data
    # --------------------------------------------------------------------------
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

    series_configs = [
        (g_max, e_max, L"E(\mathcal{C}_{max})", :crimson,    :utriangle),
        (g_2nd, e_2nd, L"E(\mathcal{C}_{2nd max})", :darkorange, :diamond),
        (g_min, e_min, L"E(\mathcal{C}_{min})", :royalblue,  :dtriangle),
        (g_zero, e_zero, L"E(\mathcal{C}=0)",           :lime,       :circle)
    ]

    # --------------------------------------------------------------------------
    # 4. Generate Heatmaps & Overlay Extrema Points
    # --------------------------------------------------------------------------
    fixed_title_str = join(["$k=$v" for (k, v) in fixed_variables], ", ")

    p1 = heatmap(gammas, Es, gap_mat;
        xlabel=L"\gamma", ylabel=L"E",
        title="Localiser Gap - ($(fixed_title_str))",
        colorbar_title=gap_title,
        colorbar=true,
        colormap=:plasma)
    
    p2 = heatmap(gammas, Es, chern_mat;
        xlabel=L"\gamma", ylabel=L"E",
        title="Localiser Chern Number",
        colorbar_title="Chern",
        colorbar=true,
        colormap=:RdBu)

    # Overlay global smart extrema points on both heatmaps
    for (g_pts, e_pts, lbl, clr, shp) in series_configs
        if !isempty(g_pts)
            scatter!(p1, g_pts, e_pts;
                label=lbl,
                color=clr,
                markershape=shp,
                markersize=4.0,
                markerstrokewidth=0.4,
                markerstrokecolor=:black,
                alpha=0.85,
                xlims=xlims,
                ylims=ylims
            )

            scatter!(p2, g_pts, e_pts;
                label=lbl,
                color=clr,
                markershape=shp,
                markersize=4.0,
                markerstrokewidth=0.4,
                markerstrokecolor=:black,
                alpha=0.85,
                xlims=xlims,
                ylims=ylims
            )
        end
    end

    plot!(p2, legend=:topright)

    plt = plot(p1, p2; layout=(1, 2), size=(1600, 600), margin=5Plots.mm)
    
    isdir(dirname(filename)) || mkpath(dirname(filename))
    savefig(plt, filename)

    return plt
end




function plt_specloc_gamma_vs_E_overlaid_inflections(
    df::DataFrame,
    gamma_data_pairs::Vector{<:Tuple{Real, Any}}; 
    atol::Real=1e-8,
    logscale::Bool=false,
    smooth_window::Int=15,       # 2D moving-average window size (odd integer)
    min_slope_tol::Float64=1e-3, # Gradient magnitude threshold to mask plateau noise
    filename::String="plots/localiser_gamma_E_heatmap.png",
    xlims::Union{Nothing, Tuple{Real, Real}}=nothing,
    ylims::Union{Nothing, Tuple{Real, Real}}=nothing,
    fixed_variables...
)
    # --------------------------------------------------------------------------
    # 1. Filter DataFrame by fixed parameter values
    # --------------------------------------------------------------------------
    col_names_str = string.(names(df))
    subdf = df
    for (k, v) in fixed_variables
        string(k) in col_names_str || error("Filter key $(k) is not a DataFrame column.")
        col = subdf[!, k]
        if v isa Real && eltype(col) <: Real
            mask = abs.(Float64.(col) .- Float64(v)) .<= Float64(atol)
            subdf = subdf[mask, :]
        else
            mask = col .== v
            subdf = subdf[mask, :]
        end
    end

    # --------------------------------------------------------------------------
    # 2. Extract Heatmap Data (Localiser Gap & Chern)
    # --------------------------------------------------------------------------
    gammas = sort(unique(subdf.gamma))
    Es = sort(unique(subdf.E))
    ngamma, nE = length(gammas), length(Es)

    gap_mat = fill(NaN, nE, ngamma)
    chern_mat = fill(NaN, nE, ngamma)

    for row in eachrow(subdf)
        xi = findfirst(row.gamma .== gammas)
        yi = findfirst(row.E .== Es)
        if !isnothing(xi) && !isnothing(yi)
            gap_mat[yi, xi] = row.localiser_gap
            chern_mat[yi, xi] = ismissing(row.chern_number) ? NaN : Float64(row.chern_number)
        end
    end

    if logscale
        gap_mat = log10.(gap_mat)
        gap_title = "log10(min|λ|)"
    else
        gap_title = "min|λ|"
    end

    # --------------------------------------------------------------------------
    # 3. Extract 2D Inflection Energies (∇²C = 0 with Gradient Mask)
    # --------------------------------------------------------------------------
    g_b1_inf, e_b1_inf = Float64[], Float64[]
    g_b2_inf, e_b2_inf = Float64[], Float64[]

    w = max(1, smooth_window)
    r = div(w, 2)

    for (gamma, band_data) in gamma_data_pairs
        g = Float64(gamma)

        kx_vals = band_data.kx_vals
        ky_vals = band_data.ky_vals
        dkx = kx_vals[2] - kx_vals[1]
        dky = ky_vals[2] - ky_vals[1]

        # Extract band arrays (2D grids)
        cum_1 = band_data.cum_chern_per_band[1, :, :]'
        E_1   = band_data.plaquette_energies[1, :, :]'
        cum_2 = band_data.cum_chern_per_band[2, :, :]'
        E_2   = band_data.plaquette_energies[2, :, :]'

        ny, nx = size(cum_1)

        for (band_idx, cum_mat, E_mat, g_buf, e_buf) in [
            (1, cum_1, E_1, g_b1_inf, e_b1_inf),
            (2, cum_2, E_2, g_b2_inf, e_b2_inf)
        ]
            # A. 2D Periodic Moving-Average Smooth
            Z_smooth = copy(cum_mat)
            if r > 0
                for iy in 1:ny, ix in 1:nx
                    y_idxs = mod1.((iy - r):(iy + r), ny)
                    x_idxs = mod1.((ix - r):(ix + r), nx)
                    Z_smooth[iy, ix] = mean(cum_mat[y_idxs, x_idxs])
                end
            end

            # B. Periodic Finite Differences
            Z_left  = circshift(Z_smooth, (0, 1))
            Z_right = circshift(Z_smooth, (0, -1))
            Z_up    = circshift(Z_smooth, (1, 0))
            Z_down  = circshift(Z_smooth, (-1, 0))

            # First Derivatives (Gradient Magnitude)
            C_x = (Z_right .- Z_left) ./ (2 * dkx)
            C_y = (Z_down .- Z_up) ./ (2 * dky)
            grad_mag = sqrt.(C_x.^2 .+ C_y.^2)

            # Second Derivatives & Laplacian
            C_xx = (Z_right .- 2 .* Z_smooth .+ Z_left) ./ (dkx^2)
            C_yy = (Z_down  .- 2 .* Z_smooth .+ Z_up)   ./ (dky^2)
            laplacian = C_xx .+ C_yy

            # C. Detect Zero Crossings with Gradient Masking
            # Check horizontal and vertical neighbor pairs for sign changes in Laplacian
            valid_mask = grad_mag .>= min_slope_tol

            # Horizontal zero-crossings
            L_right = circshift(laplacian, (0, -1))
            V_right = circshift(valid_mask, (0, -1))
            h_crossings = (laplacian .* L_right .<= 0) .& valid_mask .& V_right

            # Vertical zero-crossings
            L_down = circshift(laplacian, (-1, 0))
            V_down = circshift(valid_mask, (-1, 0))
            v_crossings = (laplacian .* L_down .<= 0) .& valid_mask .& V_down

            inflection_indices = findall(h_crossings .| v_crossings)

            for idx in inflection_indices
                push!(g_buf, g)
                push!(e_buf, E_mat[idx])
            end
        end
    end

    series_configs = [
        (g_b1_inf, e_b1_inf, "Band 1 Inflections", :cyan,    :circle),
        (g_b2_inf, e_b2_inf, "Band 2 Inflections", :magenta, :diamond)
    ]

    # --------------------------------------------------------------------------
    # 4. Generate Heatmaps & Overlay Scatter Points
    # --------------------------------------------------------------------------
    fixed_title_str = join(["$k=$v" for (k, v) in fixed_variables], ", ")

    p1 = heatmap(gammas, Es, gap_mat;
        xlabel=L"\gamma", ylabel=L"E",
        title="Localiser Gap - ($(fixed_title_str))",
        colorbar_title=gap_title,
        colorbar=true,
        colormap=:plasma)
    
    p2 = heatmap(gammas, Es, chern_mat;
        xlabel=L"\gamma", ylabel=L"E",
        title="Localiser Chern Number",
        colorbar_title="Chern",
        colorbar=true,
        colormap=:RdBu)

    for (g_pts, e_pts, lbl, clr, shp) in series_configs
        if !isempty(g_pts)
            scatter!(p1, g_pts, e_pts;
                label=lbl,
                color=clr,
                markershape=shp,
                markersize=2.5,
                markerstrokewidth=0.2,
                markerstrokecolor=:black,
                alpha=0.7,
                xlims=xlims,
                ylims=ylims
            )

            scatter!(p2, g_pts, e_pts;
                label=lbl,
                color=clr,
                markershape=shp,
                markersize=2.5,
                markerstrokewidth=0.2,
                markerstrokecolor=:black,
                alpha=0.7,
                xlims=xlims,
                ylims=ylims
            )
        end
    end

    plot!(p2, legend=:topright)

    plt = plot(p1, p2; layout=(1, 2), size=(1600, 600), margin=5Plots.mm)
    
    isdir(dirname(filename)) || mkpath(dirname(filename))
    savefig(plt, filename)

    return plt
end


function plt_specloc_gamma_vs_E_overlaid_euler_char_(
    df::DataFrame,
    crit_chern_data_list::Union{Vector{<:Tuple{Real, Real, Real, Real, Any, Any}}, Vector{<:Tuple{Real, Real, Real, Real, Any}}, Nothing}=nothing; 
    atol::Real=1e-8,
    logscale::Bool=false,
    filename::String="plots/localiser_gamma_E_heatmap_euler.png",
    xlims::Union{Nothing, Tuple{Real, Real}}=nothing,
    ylims::Union{Nothing, Tuple{Real, Real}}=nothing,
    clims_overide::Union{Nothing, Tuple{Real, Real}}=nothing,
    fixed_variables...
)
    # --------------------------------------------------------------------------
    # 1. Filter DataFrame by fixed parameter values
    # --------------------------------------------------------------------------
    subdf = filter_df_by_fixed_variables(df, fixed_variables; atol=atol)

    if nrow(subdf) == 0
        @warn "No data found matching parameters: $fixed_variables"
        return nothing
    end

    # --------------------------------------------------------------------------
    # 2. Extract Spectral Localiser Gap Data
    # --------------------------------------------------------------------------
    gammas = sort(unique(subdf.gamma))
    Es = sort(unique(subdf.E))
    ngamma, nE = length(gammas), length(Es)

    gap_mat = fill(NaN, nE, ngamma)

    for row in eachrow(subdf)
        xi = findfirst(==(row.gamma), gammas)
        yi = findfirst(==(row.E), Es)
        
        if !isnothing(xi) && !isnothing(yi)
            gap_val = if hasproperty(row, :localiser_gap) && !ismissing(row.localiser_gap)
                Float64(row.localiser_gap)
            elseif hasproperty(row, :low_lying_evals) && !ismissing(row.low_lying_evals)
                minimum(abs.(vec(row.low_lying_evals)))
            else
                NaN
            end

            gap_mat[yi, xi] = gap_val
        end
    end

    if logscale
        gap_mat = log10.(max.(gap_mat, 1e-12))
        gap_title = L"\log_{10}(\text{min}|\lambda|)"
    else
        gap_title = L"\text{min}|\lambda|"
    end

    valid_gaps = filter(!isnan, vec(gap_mat))
    c_min, c_max = isempty(valid_gaps) ? (0.0, 1.0) : (minimum(valid_gaps), maximum(valid_gaps))
    if c_min == c_max
        c_min -= 1e-5
        c_max += 1e-5
    end

    if !isnothing(clims_overide)
        c_min, c_max = clims_overide
    end

    fixed_title_str = join(["$k=$v" for (k, v) in fixed_variables], ", ")

    # --------------------------------------------------------------------------
    # 3. Create Left Panel: Spectral Localiser Gap Heatmap
    # --------------------------------------------------------------------------
    p1 = heatmap(gammas, Es, gap_mat;
        clims = (c_min, c_max),
        xlabel = L"\gamma", ylabel = L"E",
        title = "Spectral Localiser Gap - ($fixed_title_str)",
        colorbar_title = gap_title,
        colorbar = true,
        colormap = :plasma,
        aspect_ratio = :auto,
        xlims = xlims,
        ylims = ylims)
    
    # --------------------------------------------------------------------------
    # 4. Extract Euler Characteristic Data & Create Right Panel
    #    - Build Euler heatmap (right) from crit_chern_data_list
    #    - Map Euler grid onto the spectral-localiser grid (nearest average)
    #      so contours can be drawn exactly on the left panel's axes.
    # --------------------------------------------------------------------------
    p2 = plot(title="Euler Characteristic heatmap unavailable")

    if !isnothing(crit_chern_data_list) && !isempty(crit_chern_data_list)
        filtered_list = filter(item -> begin
            match_all = true
            param_map = Dict(:A => item[1], :B => item[2], :m => item[3])
            for (k, v) in fixed_variables
                if haskey(param_map, k) && v isa Real
                    if !isapprox(Float64(param_map[k]), Float64(v); atol=atol)
                        match_all = false
                        break
                    end
                end
            end
            return match_all
        end, crit_chern_data_list)

        if !isempty(filtered_list)
            euler_payloads = Dict{Float64, Tuple{Vector{Float64}, Vector{Float64}}}()
            for item in filtered_list
                gamma_val = Float64(item[4])
                e_char_data = length(item) >= 6 ? item[6] : nothing
                payload = extract_euler_payload(e_char_data)

                if isnothing(payload)
                    payload = compute_euler_payload_from_band_data(item[5])
                end

                if !isnothing(payload)
                    euler_payloads[gamma_val] = payload
                end
            end

            if !isempty(euler_payloads)
                # Canonical Euler grid (may differ from `gammas`, `Es`)
                euler_gammas = sort(collect(keys(euler_payloads)))
                euler_Es = euler_payloads[euler_gammas[1]][1]
                n_e_g = length(euler_gammas)
                n_e_E = length(euler_Es)
                euler_mat = fill(NaN, n_e_E, n_e_g)

                for (gi, gamma_val) in enumerate(euler_gammas)
                    e_vals, chi_vals = euler_payloads[gamma_val]
                    # If the energy mesh matches (within tolerance) use directly
                    if length(e_vals) == n_e_E && all(isapprox.(e_vals, euler_Es; atol=1e-8))
                        euler_mat[:, gi] .= chi_vals
                    else
                        # Attempt to resample/interpolate chi_vals onto canonical euler_Es
                        # Ensure arrays are Float64 and sorted by energy
                        ev = Float64.(e_vals)
                        ch = Float64.(chi_vals)
                        sort_idx = sortperm(ev)
                        evs = ev[sort_idx]
                        chs = ch[sort_idx]

                        interp_col = fill(NaN, n_e_E)
                        for (ti, Et) in enumerate(euler_Es)
                            # find interval [evs[j], evs[j+1]] containing Et
                            j = findlast(x -> x <= Et, evs)
                            if !isnothing(j) && j < length(evs)
                                x0, x1 = evs[j], evs[j+1]
                                y0, y1 = chs[j], chs[j+1]
                                if x1 != x0
                                    t = (Et - x0) / (x1 - x0)
                                    interp_col[ti] = (1 - t) * y0 + t * y1
                                else
                                    interp_col[ti] = y0
                                end
                            elseif !isnothing(j) && j == length(evs) && isapprox(Et, evs[end]; atol=1e-8)
                                interp_col[ti] = chs[end]
                            else
                                # If Et is below the first evs value, try nearest neighbor
                                j2 = findfirst(x -> x >= Et, evs)
                                if !isnothing(j2)
                                    interp_col[ti] = chs[j2]
                                else
                                    interp_col[ti] = NaN
                                end
                            end
                        end

                        euler_mat[:, gi] .= interp_col
                        # Silently interpolate - no warning needed
                    end
                end

                valid_count = count(!isnan, euler_mat)
                if valid_count > 0
                    plot_ylims = isnothing(ylims) ? (minimum(euler_Es), maximum(euler_Es)) : ylims

                    # Right panel: bare Euler heatmap (no contours)
                    p2 = heatmap(euler_gammas, euler_Es, euler_mat;
                        xlabel = L"\gamma", ylabel = L"E",
                        title = "Euler Characteristic " * L"\chi",
                        colorbar_title = L"\chi",
                        colorbar = true,
                        colormap = :viridis,
                        aspect_ratio = :auto,
                        xlims = xlims,
                        ylims = plot_ylims)

                    # Map Euler data onto the spectral-localiser grid so contours
                    # are drawn in the exact same plotting axes as `p1`.
                    nx = length(gammas)
                    ny = length(Es)
                    sum_mat = zeros(Float64, ny, nx)
                    count_mat = zeros(Int, ny, nx)

                    for (egj, eg) in enumerate(euler_gammas)
                        # find nearest x index in gammas
                        xdist = abs.(gammas .- eg)
                        xidx = findmin(xdist)[2]
                        for (eyi, eE) in enumerate(euler_Es)
                            ydist = abs.(Es .- eE)
                            yidx = findmin(ydist)[2]
                            v = euler_mat[eyi, egj]
                            if !isnan(v)
                                sum_mat[yidx, xidx] += Float64(v)
                                count_mat[yidx, xidx] += 1
                            end
                        end
                    end

                    mapped_euler = fill(NaN, ny, nx)
                    for j in 1:ny, i in 1:nx
                        if count_mat[j, i] > 0
                            mapped_euler[j, i] = sum_mat[j, i] / count_mat[j, i]
                        end
                    end

                    # Simple edge detection: find where Euler characteristic changes value
                    edge_gamma_pts = Float64[]
                    edge_E_pts = Float64[]

                    nE_grid = length(Es)
                    ngamma_grid = length(gammas)

                    # Scan the grid and detect edges (transitions between different integer values)
                    for j in 1:(nE_grid-1)
                        for i in 1:(ngamma_grid-1)
                            curr_val = mapped_euler[j, i]
                            if !isnan(curr_val)
                                # Check horizontal neighbor (right)
                                right_val = mapped_euler[j, i+1]
                                if !isnan(right_val) && round(Int, curr_val) != round(Int, right_val)
                                    # Edge between cells i and i+1
                                    push!(edge_gamma_pts, (gammas[i] + gammas[i+1]) / 2)
                                    push!(edge_E_pts, Es[j])
                                end
                                
                                # Check vertical neighbor (down)
                                down_val = mapped_euler[j+1, i]
                                if !isnan(down_val) && round(Int, curr_val) != round(Int, down_val)
                                    # Edge between cells j and j+1
                                    push!(edge_gamma_pts, gammas[i])
                                    push!(edge_E_pts, (Es[j] + Es[j+1]) / 2)
                                end
                            end
                        end
                    end

                    # Plot edge points as white scatter on the spectral localiser plot
                    if !isempty(edge_gamma_pts)
                        scatter!(p1, edge_gamma_pts, edge_E_pts;
                            marker = :circle,
                            markersize = 2.0,
                            markercolor = :red,
                            markerstrokewidth = 0,
                            alpha = 0.8,
                            label = false)
                    end
                else
                    @warn "Euler matrix has no valid entries after assembly"
                end
            else
                @warn "No valid Euler payloads found for selected parameters"
            end
        end
    end

    plt = plot(p1, p2; layout=(1, 2), size=(1600, 600), margin=5Plots.mm)
    
    if !isempty(filename)
        isdir(dirname(filename)) || mkpath(dirname(filename))
        savefig(plt, filename)
    end

    return plt
end

function plt_specloc_gamma_vs_E_overlaid_topological_quantity(
    df::DataFrame,
    crit_chern_data_list::Union{Vector{<:Tuple{Real, Real, Real, Real, Any, Any}}, Vector{<:Tuple{Real, Real, Real, Real, Any}}, Nothing}=nothing; 
    plot_type::Symbol=:holes,
    atol::Real=1e-8,
    logscale::Bool=false,
    filename::String="plots/localiser_gamma_E_heatmap_topology.png",
    xlims::Union{Nothing, Tuple{Real, Real}}=nothing,
    ylims::Union{Nothing, Tuple{Real, Real}}=nothing,
    clims_overide::Union{Nothing, Tuple{Real, Real}}=nothing,
    fixed_variables...
)
    # Validate plot_type parameter
    plot_type in (:holes, :pockets) || error("plot_type must be :holes or :pockets, got :$plot_type")
    
    # --------------------------------------------------------------------------
    # 1. Filter DataFrame by fixed parameter values
    # --------------------------------------------------------------------------
    subdf = filter_df_by_fixed_variables(df, fixed_variables; atol=atol)

    if nrow(subdf) == 0
        @warn "No data found matching parameters: $fixed_variables"
        return nothing
    end

    # --------------------------------------------------------------------------
    # 2. Extract Spectral Localiser Gap Data
    # --------------------------------------------------------------------------
    gammas = sort(unique(subdf.gamma))
    Es = sort(unique(subdf.E))
    ngamma, nE = length(gammas), length(Es)

    gap_mat = fill(NaN, nE, ngamma)

    for row in eachrow(subdf)
        xi = findfirst(==(row.gamma), gammas)
        yi = findfirst(==(row.E), Es)
        
        if !isnothing(xi) && !isnothing(yi)
            gap_val = if hasproperty(row, :localiser_gap) && !ismissing(row.localiser_gap)
                Float64(row.localiser_gap)
            elseif hasproperty(row, :low_lying_evals) && !ismissing(row.low_lying_evals)
                minimum(abs.(vec(row.low_lying_evals)))
            else
                NaN
            end

            gap_mat[yi, xi] = gap_val
        end
    end

    if logscale
        gap_mat = log10.(max.(gap_mat, 1e-12))
        gap_title = L"\log_{10}(\text{min}|\lambda|)"
    else
        gap_title = L"\text{min}|\lambda|"
    end

    valid_gaps = filter(!isnan, vec(gap_mat))
    c_min, c_max = isempty(valid_gaps) ? (0.0, 1.0) : (minimum(valid_gaps), maximum(valid_gaps))
    if c_min == c_max
        c_min -= 1e-5
        c_max += 1e-5
    end

    if !isnothing(clims_overide)
        c_min, c_max = clims_overide
    end

    fixed_title_str = join(["$k=$v" for (k, v) in fixed_variables], ", ")

    # --------------------------------------------------------------------------
    # 3. Create Left Panel: Spectral Localiser Gap Heatmap
    # --------------------------------------------------------------------------
    p1 = heatmap(gammas, Es, gap_mat;
        clims = (c_min, c_max),
        xlabel = L"\gamma", ylabel = L"E",
        title = "Spectral Localiser Gap - ($fixed_title_str)",
        colorbar_title = gap_title,
        colorbar = true,
        colormap = :plasma,
        aspect_ratio = :auto,
        xlims = xlims,
        ylims = ylims)
    
    # --------------------------------------------------------------------------
    # 4. Extract Topological Quantity Data (Holes or Pockets) & Create Right Panel
    # --------------------------------------------------------------------------
    p2 = plot(title="$(plot_type |> string |> titlecase) heatmap unavailable")
    
    quantity_label = plot_type == :holes ? "Number of Holes" : "Number of Pockets"
    quantity_short = plot_type == :holes ? L"N_h" : L"N_p"

    if !isnothing(crit_chern_data_list) && !isempty(crit_chern_data_list)
        filtered_list = filter(item -> begin
            match_all = true
            param_map = Dict(:A => item[1], :B => item[2], :m => item[3])
            for (k, v) in fixed_variables
                if haskey(param_map, k) && v isa Real
                    if !isapprox(Float64(param_map[k]), Float64(v); atol=atol)
                        match_all = false
                        break
                    end
                end
            end
            return match_all
        end, crit_chern_data_list)

        if !isempty(filtered_list)
            quantity_payloads = Dict{Float64, Tuple{Vector{Float64}, Vector{Float64}}}()
            for item in filtered_list
                gamma_val = Float64(item[4])
                
                # Try to compute from band data
                payload = compute_holes_pockets_payload_from_band_data(item[5]; quantity=plot_type)

                if !isnothing(payload)
                    quantity_payloads[gamma_val] = payload
                end
            end

            if !isempty(quantity_payloads)
                # Canonical grid (may differ from `gammas`, `Es`)
                quantity_gammas = sort(collect(keys(quantity_payloads)))
                quantity_Es = quantity_payloads[quantity_gammas[1]][1]
                n_q_g = length(quantity_gammas)
                n_q_E = length(quantity_Es)
                quantity_mat = fill(NaN, n_q_E, n_q_g)

                for (gi, gamma_val) in enumerate(quantity_gammas)
                    e_vals, qty_vals = quantity_payloads[gamma_val]
                    # If the energy mesh matches (within tolerance) use directly
                    if length(e_vals) == n_q_E && all(isapprox.(e_vals, quantity_Es; atol=1e-8))
                        quantity_mat[:, gi] .= qty_vals
                    else
                        # Attempt to resample/interpolate qty_vals onto canonical quantity_Es
                        ev = Float64.(e_vals)
                        qty = Float64.(qty_vals)
                        sort_idx = sortperm(ev)
                        evs = ev[sort_idx]
                        qtys = qty[sort_idx]

                        interp_col = fill(NaN, n_q_E)
                        for (ti, Et) in enumerate(quantity_Es)
                            # find interval [evs[j], evs[j+1]] containing Et
                            j = findlast(x -> x <= Et, evs)
                            if !isnothing(j) && j < length(evs)
                                x0, x1 = evs[j], evs[j+1]
                                y0, y1 = qtys[j], qtys[j+1]
                                if x1 != x0
                                    t = (Et - x0) / (x1 - x0)
                                    interp_col[ti] = (1 - t) * y0 + t * y1
                                else
                                    interp_col[ti] = y0
                                end
                            elseif !isnothing(j) && j == length(evs) && isapprox(Et, evs[end]; atol=1e-8)
                                interp_col[ti] = qtys[end]
                            else
                                # If Et is below the first evs value, try nearest neighbor
                                j2 = findfirst(x -> x >= Et, evs)
                                if !isnothing(j2)
                                    interp_col[ti] = qtys[j2]
                                else
                                    interp_col[ti] = NaN
                                end
                            end
                        end

                        quantity_mat[:, gi] .= interp_col
                    end
                end

                valid_count = count(!isnan, quantity_mat)
                if valid_count > 0
                    plot_ylims = isnothing(ylims) ? (minimum(quantity_Es), maximum(quantity_Es)) : ylims

                    # Right panel: bare topological quantity heatmap
                    p2 = heatmap(quantity_gammas, quantity_Es, quantity_mat;
                        xlabel = L"\gamma", ylabel = L"E",
                        title = quantity_label,
                        colorbar_title = quantity_short,
                        colorbar = true,
                        colormap = :viridis,
                        aspect_ratio = :auto,
                        xlims = xlims,
                        ylims = plot_ylims)

                    # Map topological data onto the spectral-localiser grid so contours
                    # are drawn in the exact same plotting axes as `p1`.
                    nx = length(gammas)
                    ny = length(Es)
                    sum_mat = zeros(Float64, ny, nx)
                    count_mat = zeros(Int, ny, nx)

                    for (qgj, qg) in enumerate(quantity_gammas)
                        # find nearest x index in gammas
                        xdist = abs.(gammas .- qg)
                        xidx = findmin(xdist)[2]
                        for (qyi, qE) in enumerate(quantity_Es)
                            ydist = abs.(Es .- qE)
                            yidx = findmin(ydist)[2]
                            v = quantity_mat[qyi, qgj]
                            if !isnan(v)
                                sum_mat[yidx, xidx] += Float64(v)
                                count_mat[yidx, xidx] += 1
                            end
                        end
                    end

                    mapped_quantity = fill(NaN, ny, nx)
                    for j in 1:ny, i in 1:nx
                        if count_mat[j, i] > 0
                            mapped_quantity[j, i] = sum_mat[j, i] / count_mat[j, i]
                        end
                    end

                    # Simple edge detection: find where topological quantity changes value
                    edge_gamma_pts = Float64[]
                    edge_E_pts = Float64[]

                    nE_grid = length(Es)
                    ngamma_grid = length(gammas)

                    # Scan the grid and detect edges (transitions between different integer values)
                    for j in 1:(nE_grid-1)
                        for i in 1:(ngamma_grid-1)
                            curr_val = mapped_quantity[j, i]
                            if !isnan(curr_val)
                                # Check horizontal neighbor (right)
                                right_val = mapped_quantity[j, i+1]
                                if !isnan(right_val) && round(Int, curr_val) != round(Int, right_val)
                                    # Edge between cells i and i+1
                                    push!(edge_gamma_pts, (gammas[i] + gammas[i+1]) / 2)
                                    push!(edge_E_pts, Es[j])
                                end
                                
                                # Check vertical neighbor (down)
                                down_val = mapped_quantity[j+1, i]
                                if !isnan(down_val) && round(Int, curr_val) != round(Int, down_val)
                                    # Edge between cells j and j+1
                                    push!(edge_gamma_pts, gammas[i])
                                    push!(edge_E_pts, (Es[j] + Es[j+1]) / 2)
                                end
                            end
                        end
                    end

                    # Plot edge points as red scatter on the spectral localiser plot
                    if !isempty(edge_gamma_pts)
                        scatter!(p1, edge_gamma_pts, edge_E_pts;
                            marker = :circle,
                            markersize = 2.0,
                            markercolor = :red,
                            markerstrokewidth = 0,
                            alpha = 0.8,
                            label = false)
                    end
                else
                    @warn "$(plot_type |> string |> titlecase) matrix has no valid entries after assembly"
                end
            else
                @warn "No valid $(plot_type) payloads found for selected parameters"
            end
        end
    end

    plt = plot(p1, p2; layout=(1, 2), size=(1600, 600), margin=5Plots.mm)
    
    if !isempty(filename)
        isdir(dirname(filename)) || mkpath(dirname(filename))
        savefig(plt, filename)
    end

    return plt
end

function plt_specloc_gamma_vs_E_overlaid_holes_pockets_2x2(
    df::DataFrame,
    crit_chern_data_list::Union{Vector{<:Tuple{Real, Real, Real, Real, Any, Any}}, Vector{<:Tuple{Real, Real, Real, Real, Any}}, Nothing}=nothing; 
    atol::Real=1e-8,
    logscale::Bool=false,
    filename::String="plots/localiser_gamma_E_heatmap_holes_pockets.png",
    xlims::Union{Nothing, Tuple{Real, Real}}=nothing,
    ylims::Union{Nothing, Tuple{Real, Real}}=nothing,
    clims_overide::Union{Nothing, Tuple{Real, Real}}=nothing,
    fixed_variables...
)
    # --------------------------------------------------------------------------
    # 1. Filter DataFrame by fixed parameter values
    # --------------------------------------------------------------------------
    subdf = filter_df_by_fixed_variables(df, fixed_variables; atol=atol)

    if nrow(subdf) == 0
        @warn "No data found matching parameters: $fixed_variables"
        return nothing
    end

    # --------------------------------------------------------------------------
    # 2. Extract Spectral Localiser Gap Data
    # --------------------------------------------------------------------------
    gammas = sort(unique(subdf.gamma))
    Es = sort(unique(subdf.E))
    ngamma, nE = length(gammas), length(Es)

    gap_mat = fill(NaN, nE, ngamma)

    for row in eachrow(subdf)
        xi = findfirst(==(row.gamma), gammas)
        yi = findfirst(==(row.E), Es)
        
        if !isnothing(xi) && !isnothing(yi)
            gap_val = if hasproperty(row, :localiser_gap) && !ismissing(row.localiser_gap)
                Float64(row.localiser_gap)
            elseif hasproperty(row, :low_lying_evals) && !ismissing(row.low_lying_evals)
                minimum(abs.(vec(row.low_lying_evals)))
            else
                NaN
            end

            gap_mat[yi, xi] = gap_val
        end
    end

    if logscale
        gap_mat = log10.(max.(gap_mat, 1e-12))
        gap_title = L"\log_{10}(\text{min}|\lambda|)"
    else
        gap_title = L"\text{min}|\lambda|"
    end

    valid_gaps = filter(!isnan, vec(gap_mat))
    c_min, c_max = isempty(valid_gaps) ? (0.0, 1.0) : (minimum(valid_gaps), maximum(valid_gaps))
    if c_min == c_max
        c_min -= 1e-5
        c_max += 1e-5
    end

    if !isnothing(clims_overide)
        c_min, c_max = clims_overide
    end

    fixed_title_str = join(["$k=$v" for (k, v) in fixed_variables], ", ")

    # --------------------------------------------------------------------------
    # Helper function to create heatmap with contours for a given quantity
    # --------------------------------------------------------------------------
    function create_quantity_panels(quantity::Symbol, quantity_name::String, quantity_short::Union{String, LaTeXString})
        # Create LEFT panel: spectral localiser gap heatmap
        p_left = heatmap(gammas, Es, gap_mat;
            clims = (c_min, c_max),
            xlabel = L"\gamma", ylabel = L"E",
            title = "Spectral Localiser Gap - ($fixed_title_str)",
            colorbar_title = gap_title,
            colorbar = true,
            colormap = :plasma,
            aspect_ratio = :auto,
            xlims = xlims,
            ylims = ylims)
        
        # Create RIGHT panel: topological quantity heatmap
        p_right = plot(title="$(quantity_name) heatmap unavailable")

        if !isnothing(crit_chern_data_list) && !isempty(crit_chern_data_list)
            filtered_list = filter(item -> begin
                match_all = true
                param_map = Dict(:A => item[1], :B => item[2], :m => item[3])
                for (k, v) in fixed_variables
                    if haskey(param_map, k) && v isa Real
                        if !isapprox(Float64(param_map[k]), Float64(v); atol=atol)
                            match_all = false
                            break
                        end
                    end
                end
                return match_all
            end, crit_chern_data_list)

            if !isempty(filtered_list)
                quantity_payloads = Dict{Float64, Tuple{Vector{Float64}, Vector{Float64}}}()
                for item in filtered_list
                    gamma_val = Float64(item[4])
                    e_char_data = length(item) >= 6 ? item[6] : nothing
                    
                    # Extract the appropriate payload based on quantity
                    payload = if quantity == :holes
                        extract_holes_payload(e_char_data)
                    elseif quantity == :pockets
                        extract_pockets_payload(e_char_data)
                    else
                        nothing
                    end

                    if !isnothing(payload)
                        quantity_payloads[gamma_val] = payload
                    end
                end

                if !isempty(quantity_payloads)
                    # Canonical grid
                    quantity_gammas = sort(collect(keys(quantity_payloads)))
                    quantity_Es = quantity_payloads[quantity_gammas[1]][1]
                    n_q_g = length(quantity_gammas)
                    n_q_E = length(quantity_Es)
                    quantity_mat = fill(NaN, n_q_E, n_q_g)

                    for (gi, gamma_val) in enumerate(quantity_gammas)
                        e_vals, qty_vals = quantity_payloads[gamma_val]
                        # If energy mesh matches, use directly
                        if length(e_vals) == n_q_E && all(isapprox.(e_vals, quantity_Es; atol=1e-8))
                            quantity_mat[:, gi] .= qty_vals
                        else
                            # Resample/interpolate onto canonical grid
                            ev = Float64.(e_vals)
                            qty = Float64.(qty_vals)
                            sort_idx = sortperm(ev)
                            evs = ev[sort_idx]
                            qtys = qty[sort_idx]

                            interp_col = fill(NaN, n_q_E)
                            for (ti, Et) in enumerate(quantity_Es)
                                j = findlast(x -> x <= Et, evs)
                                if !isnothing(j) && j < length(evs)
                                    x0, x1 = evs[j], evs[j+1]
                                    y0, y1 = qtys[j], qtys[j+1]
                                    if x1 != x0
                                        t = (Et - x0) / (x1 - x0)
                                        interp_col[ti] = (1 - t) * y0 + t * y1
                                    else
                                        interp_col[ti] = y0
                                    end
                                elseif !isnothing(j) && j == length(evs) && isapprox(Et, evs[end]; atol=1e-8)
                                    interp_col[ti] = qtys[end]
                                else
                                    j2 = findfirst(x -> x >= Et, evs)
                                    if !isnothing(j2)
                                        interp_col[ti] = qtys[j2]
                                    else
                                        interp_col[ti] = NaN
                                    end
                                end
                            end

                            quantity_mat[:, gi] .= interp_col
                        end
                    end

                    valid_count = count(!isnan, quantity_mat)
                    if valid_count > 0
                        plot_ylims = isnothing(ylims) ? (minimum(quantity_Es), maximum(quantity_Es)) : ylims

                        # Right panel: topological quantity heatmap
                        p_right = heatmap(quantity_gammas, quantity_Es, quantity_mat;
                            xlabel = L"\gamma", ylabel = L"E",
                            title = quantity_name,
                            colorbar_title = quantity_short,
                            colorbar = true,
                            colormap = :viridis,
                            aspect_ratio = :auto,
                            xlims = xlims,
                            ylims = plot_ylims)

                        # Map topological data onto spectral-localiser grid for contours
                        nx = length(gammas)
                        ny = length(Es)
                        sum_mat = zeros(Float64, ny, nx)
                        count_mat = zeros(Int, ny, nx)

                        for (qgj, qg) in enumerate(quantity_gammas)
                            xdist = abs.(gammas .- qg)
                            xidx = findmin(xdist)[2]
                            for (qyi, qE) in enumerate(quantity_Es)
                                ydist = abs.(Es .- qE)
                                yidx = findmin(ydist)[2]
                                v = quantity_mat[qyi, qgj]
                                if !isnan(v)
                                    sum_mat[yidx, xidx] += Float64(v)
                                    count_mat[yidx, xidx] += 1
                                end
                            end
                        end

                        mapped_quantity = fill(NaN, ny, nx)
                        for j in 1:ny, i in 1:nx
                            if count_mat[j, i] > 0
                                mapped_quantity[j, i] = sum_mat[j, i] / count_mat[j, i]
                            end
                        end

                        # Edge detection: find where quantity changes value
                        edge_gamma_pts = Float64[]
                        edge_E_pts = Float64[]

                        nE_grid = length(Es)
                        ngamma_grid = length(gammas)

                        for j in 1:(nE_grid-1)
                            for i in 1:(ngamma_grid-1)
                                curr_val = mapped_quantity[j, i]
                                if !isnan(curr_val)
                                    # Check right neighbor
                                    right_val = mapped_quantity[j, i+1]
                                    if !isnan(right_val) && round(Int, curr_val) != round(Int, right_val)
                                        push!(edge_gamma_pts, (gammas[i] + gammas[i+1]) / 2)
                                        push!(edge_E_pts, Es[j])
                                    end
                                    
                                    # Check down neighbor
                                    down_val = mapped_quantity[j+1, i]
                                    if !isnan(down_val) && round(Int, curr_val) != round(Int, down_val)
                                        push!(edge_gamma_pts, gammas[i])
                                        push!(edge_E_pts, (Es[j] + Es[j+1]) / 2)
                                    end
                                end
                            end
                        end

                        # Overlay contours on left panel
                        if !isempty(edge_gamma_pts)
                            scatter!(p_left, edge_gamma_pts, edge_E_pts;
                                marker = :circle,
                                markersize = 2.0,
                                markercolor = :red,
                                markerstrokewidth = 0,
                                alpha = 0.8,
                                label = false)
                        end
                    end
                end
            end
        end

        return (p_left, p_right)
    end

    # --------------------------------------------------------------------------
    # 3. Create 2x2 grid: holes (top), pockets (bottom)
    # --------------------------------------------------------------------------
    p_holes_left, p_holes_right = create_quantity_panels(:holes, "Number of Holes", L"N_h")
    p_pockets_left, p_pockets_right = create_quantity_panels(:pockets, "Number of Pockets", L"N_p")

    # Create 2x2 layout
    plt = plot(
        p_holes_left, p_holes_right,
        p_pockets_left, p_pockets_right;
        layout = (2, 2),
        size = (1600, 1200),
        margin = 5Plots.mm
    )

    if !isempty(filename)
        isdir(dirname(filename)) || mkpath(dirname(filename))
        savefig(plt, filename)
    end

    return plt
end



for A in valid_Avals, B in valid_Bvals, m in valid_mvals
    
    # Extract (gamma, band_data) pairs for current slice
    # sub_chern_pairs = [(item[4], item[5]) for item in crit_chern_valid_list if 
    #     isapprox(item[1], A; atol=1e-6) &&
    #     isapprox(item[2], B; atol=1e-6) &&
    #     isapprox(item[3], m; atol=1e-6)
    # ]

    # Extract 2-element tuples (gamma, band_data) for functions expecting 2 items
    sub_chern_pairs = [(item[4], item[5]) for item in crit_chern_valid_list if 
        isapprox(item[1], A; atol=1e-6) &&
        isapprox(item[2], B; atol=1e-6) &&
        isapprox(item[3], m; atol=1e-6)
    ]

    # Extract full 6-element tuples (A, B, m, gamma, band_data, euler_char_data)
    sub_euler_tuples = [item for item in crit_chern_valid_list if 
        isapprox(item[1], A; atol=1e-6) &&
        isapprox(item[2], B; atol=1e-6) &&
        isapprox(item[3], m; atol=1e-6)
    ]

    isempty(sub_chern_pairs) && continue

    # # 1. Generate Heatmap of Localiser Gap & Chern Number with Critical Chern Contours
    # fn = joinpath(output_path, "localiser_heatmap_A$(A)_B$(B)_m$(m).png")

    # plt_specloc_gamma_vs_E_overlaid_crit_cum_chern(
    #     specloc_valid_df,
    #     sub_chern_pairs;
    #     logscale = true,
    #     target_p = 0.5,
    #     target_q = 0.0,
    #     tol = 0.005,
    #     xlims = (-3.0,3.0),
    #     ylims = (-6.0, 5.0),
    #     filename = fn,
    #     A = A, B = B, m = m
    # )

    # # 2. with three contours (C = 0, C = +1, C = -1) overlaid on the same heatmap
    # fn = joinpath(output_path, "zoomed_localiser_rrrr_heatmap_A$(A)_B$(B)_m$(m).png")

    # plt_specloc_gamma_vs_E_overlaid_crit_cum_chern_three_targets(
    #     specloc_valid_df,
    #     sub_chern_pairs;
    #     logscale = true,
    #     target_p = 0.5,
    #     target_q = 0.0,
    #     target_r=1.0,
    #     tol = 0.005,
    #     xlims = (0.0,2.0),
    #     ylims = (0.0, 3.0),
    #     filename = fn,
    #     A = A, B = B, m = m
    # )

    # # 3. Generate Heatmap of Localiser Gap & Chern Number with Joint Band Critical Chern Contours
    # fn = joinpath(output_path, "localiser_joint_heatmap_A$(A)_B$(B)_m$(m).png")

    # plt_specloc_gamma_vs_E_overlaid_joint_band_crit_cum_chern(
    #     specloc_valid_df,
    #     sub_chern_pairs;
    #     logscale = true,
    #     target_p = 0.5,
    #     target_q = 0.0,
    #     tol = 0.005,
    #     xlims = (-3.0,3.0),
    #     ylims = (-6.0, 5.0),
    #     filename = fn,
    #     A = A, B = B, m = m
    # )

    # # 4. max only 
    # fn = joinpath(output_path, "localiser_max_cum_chern_heatmap_A$(A)_B$(B)_m$(m).png")

    # plt_specloc_gamma_vs_E_overlaid_joint_band_max_cum_chern(
    #     specloc_valid_df,
    #     sub_chern_pairs;
    #     logscale = true,
    #     target_zero = 0.0,
    #     tol = 0.005,
    #     xlims = (-3.0,3.0),
    #     ylims = (-6.0, 5.0),
    #     filename = fn,
    #     A = A, B = B, m = m
    # )

    # # 5. all smart extrema (max, 2nd max, min, zero) overlaid on the same heatmap
    # fn = joinpath(output_path, "localiser_smart_extrema_heatmap_A$(A)_B$(B)_m$(m).png")

    # plt_specloc_gamma_vs_E_overlaid_joint_band_smart_extrema(
    #     specloc_valid_df,
    #     sub_chern_pairs;
    #     logscale = true,
    #     peak_tol = 0.01,
    #     e_cluster_tol = 0.1,
    #     flat_span_tol = 0.02,
    #     zero_tol = 1e-3,
    #     xlims = (-3.0,3.0),
    #     ylims = (-6.0, 5.0),
    #     filename = fn,
    #     A = A, B = B, m = m
    # )


    # # 6. inflection points (∇²C = 0) overlaid on the same heatmap
    # fn = joinpath(output_path, "localiser_inflection_heatmap_A$(A)_B$(B)_m$(m).png")

    # plt_specloc_gamma_vs_E_overlaid_inflections(
    #     specloc_valid_df,
    #     sub_chern_pairs;
    #     logscale = true,
    #     smooth_window = 10,
    #     min_slope_tol = 1e-5,
    #     xlims = (-3.0,3.0),
    #     ylims = (-6.0, 5.0),
    #     filename = fn,
    #     A = A, B = B, m = m
    # )

    # # # 7. Euler characteristic contours overlaid on the same heatmap
    # fn = joinpath(output_path, "zoomed_nolog_localiser_overlaid_euler_char_heatmap_A$(A)_B$(B)_m$(m).png")
    
    # plt_specloc_gamma_vs_E_overlaid_euler_char_(
    #     specloc_valid_df,
    #     sub_euler_tuples;
    #     logscale = false,
    #     xlims = (0.0, 2.0), #(-3.0,3.0),
    #     ylims = (0.0, 3.0), #(-3.0,3.0),
    #     clims_overide = (0.0, 0.01),
    #     filename = fn,
    #     A = A, B = B, m = m, kappa = target_kappa
    # )

    # # # 8. plt either N_holes or N_pockets (instead of full Euler characteristic) overlaid on the same heatmap
    # type=:holes
    # fn = joinpath(output_path, "zoomed_localiser_overlaid_$(type)_heatmap_A$(A)_B$(B)_m$(m).png")
    # println("Saving $(type) heatmap to: $fn")

    # plt_specloc_gamma_vs_E_overlaid_topological_quantity(
    #     specloc_valid_df,
    #     sub_euler_tuples;
    #     plot_type = type,
    #     logscale = true,
    #     xlims = (0.0, 2.0), #(-3.0,3.0),
    #     ylims = (0.0, 3.0), #(-3.0,3.0),
    #     clims_overide = nothing, #(0.0, 0.01),
    #     filename = fn,
    #     A = A, B = B, m = m, kappa = target_kappa
    # )


    # 8b.
    fn = joinpath(output_path, "zoomed_localiser_overlaid_holes_pockets_heatmap_A$(A)_B$(B)_m$(m).png")

    plt_specloc_gamma_vs_E_overlaid_holes_pockets_2x2(
        specloc_valid_df,
        crit_chern_valid_list;
        logscale = true,
        xlims = (0.0, 2.0),
        ylims = (0.0, 3.0),
        filename = fn,
        A = 1.0, B = 1.0, m = -1.0, kappa = target_kappa
    )

end
