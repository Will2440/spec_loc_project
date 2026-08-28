using Plots
using LaTeXStrings
using DataFrames
using JLD2
using Glob
using ProgressMeter
using Peaks


specloc_data_path = "/Users/Will/Documents/spec_loc_project/notebook_calcs/data/test_pm_gamma_m-1.5"



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


######################################################################################################
######################################################################################################

crit_chern_data_path = "/Users/Will/Documents/spec_loc_project/notebook_calcs/data/bulk_band_berry_data/bandstrucutre_data"

file_paths = glob("*.jld2", crit_chern_data_path)

# Store full tuple: (A, B, m, gamma, band_data)
crit_chern_data_list = Tuple{Float64, Float64, Float64, Float64, Any}[]

for filepath in file_paths
    fname = basename(filepath)
    
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

    push!(crit_chern_data_list, (A_val, B_val, m_val, gamma_val, band_data))
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








function plt_specloc_gamma_vs_E_overlaid_crit_cum_chern_rrrr(
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

    # --------------------------------------------------------------------------
    # 4. Generate Heatmaps & Overlay Scatter Points (Equal & Square Layout)
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

    plot!(p2, legend=:topright)

    plt = plot(p1, p2; layout=(1, 2), size=(1600, 600), margin=5Plots.mm)
    
    isdir(dirname(filename)) || mkpath(dirname(filename))
    savefig(plt, filename)

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

# function plt_specloc_gamma_vs_E_overlaid_joint_band_2nd_max_cum_chern(
#     df::DataFrame,
#     gamma_data_pairs::Vector{<:Tuple{Real, Any}}; 
#     atol::Real=1e-8,
#     logscale::Bool=false,
#     target_zero::Float64=0.0,
#     tol::Float64=0.01,
#     filename::String="plots/localiser_gamma_E_2nd_max_cum_chern_heatmap.png",
#     xlims::Union{Nothing, Tuple{Real, Real}}=nothing,
#     ylims::Union{Nothing, Tuple{Real, Real}}=nothing,
#     max_colormap::Symbol=:viridis,
#     fixed_variables...
# )
#     # --------------------------------------------------------------------------
#     # 1. Filter DataFrame by fixed parameter values
#     # --------------------------------------------------------------------------
#     col_names_str = string.(names(df))
#     subdf = df
#     for (k, v) in fixed_variables
#         string(k) in col_names_str || error("Filter key $(k) is not a DataFrame column.")
#         col = subdf[!, k]
#         if v isa Real && eltype(col) <: Real
#             mask = abs.(Float64.(col) .- Float64(v)) .<= Float64(atol)
#             subdf = subdf[mask, :]
#         else
#             mask = col .== v
#             subdf = subdf[mask, :]
#         end
#     end

#     # --------------------------------------------------------------------------
#     # 2. Extract Heatmap Data (Localiser Gap & Chern)
#     # --------------------------------------------------------------------------
#     gammas = sort(unique(subdf.gamma))
#     Es = sort(unique(subdf.E))
#     ngamma, nE = length(gammas), length(Es)

#     gap_mat = fill(NaN, nE, ngamma)
#     chern_mat = fill(NaN, nE, ngamma)

#     for row in eachrow(subdf)
#         xi = findfirst(row.gamma .== gammas)
#         yi = findfirst(row.E .== Es)
#         if !isnothing(xi) && !isnothing(yi)
#             gap_mat[yi, xi] = row.localiser_gap
#             chern_mat[yi, xi] = ismissing(row.chern_number) ? NaN : Float64(row.chern_number)
#         end
#     end

#     if logscale
#         gap_mat = log10.(gap_mat)
#         gap_title = "log10(min|λ|)"
#     else
#         gap_title = "min|λ|"
#     end

#     # --------------------------------------------------------------------------
#     # 3. Extract Scatter Data (C = 0 contours & 2nd Highest Local Max Cum Chern)
#     # --------------------------------------------------------------------------
#     g_zero, e_zero = Float64[], Float64[]
#     g_2nd_max, e_2nd_max, c_2nd_max = Float64[], Float64[], Float64[]

#     for (gamma, band_data) in gamma_data_pairs
#         g = Float64(gamma)

#         # Flatten energy across all bands
#         E_flat = vec(band_data.plaquette_energies)

#         # Compute global energy-accumulated Chern number
#         cum_flat = if hasproperty(band_data, :cum_chern_global)
#             vec(band_data.cum_chern_global)
#         elseif hasproperty(band_data, :berry_curvature)
#             b_flat = vec(band_data.berry_curvature) ./ (2π)
#             sort_idx = sortperm(E_flat)
#             cum_sorted = cumsum(b_flat[sort_idx])
#             cum_resorted = similar(cum_sorted)
#             cum_resorted[sort_idx] = cum_sorted
#             cum_resorted
#         else
#             vec(band_data.cum_chern_per_band)
#         end

#         # (a) Track critical C = 0 points
#         for i in findall(v -> abs(v - target_zero) <= tol, cum_flat)
#             push!(g_zero, g)
#             push!(e_zero, E_flat[i])
#         end

#         # (b) Sort C(E) along ascending energy to identify 1D local maxima peaks
#         sort_idx = sortperm(E_flat)
#         E_sorted = E_flat[sort_idx]
#         cum_sorted = cum_flat[sort_idx]

#         N = length(cum_sorted)
#         peak_indices = Int[]

#         if N >= 3
#             if cum_sorted[1] > cum_sorted[2]
#                 push!(peak_indices, 1)
#             end
#             for i in 2:(N-1)
#                 if cum_sorted[i] >= cum_sorted[i-1] && cum_sorted[i] >= cum_sorted[i+1]
#                     # Avoid duplicate counting on flat plateaus
#                     if cum_sorted[i] > cum_sorted[i-1] || cum_sorted[i] > cum_sorted[i+1]
#                         push!(peak_indices, i)
#                     end
#                 end
#             end
#             if cum_sorted[N] > cum_sorted[N-1]
#                 push!(peak_indices, N)
#             end
#         elseif N > 0
#             push!(peak_indices, argmax(cum_sorted))
#         end

#         # Extract distinct local peak values (descending)
#         peak_vals = cum_sorted[peak_indices]
#         unique_peak_vals = Float64[]
#         for v in sort(peak_vals, rev=true)
#             if isempty(unique_peak_vals) || all(abs.(unique_peak_vals .- v) .> 1e-4)
#                 push!(unique_peak_vals, v)
#             end
#         end

#         # Collect points corresponding to the 2nd highest local maximum (if present)
#         if length(unique_peak_vals) >= 2
#             second_max_val = unique_peak_vals[2]
#             for idx in peak_indices
#                 if abs(cum_sorted[idx] - second_max_val) <= tol
#                     push!(g_2nd_max, g)
#                     push!(e_2nd_max, E_sorted[idx])
#                     push!(c_2nd_max, second_max_val)
#                 end
#             end
#         end
#     end

#     # Color limits for the 2nd maximum cumulative Chern scale
#     clims_2nd = if !isempty(c_2nd_max)
#         cmin, cmax_val = minimum(c_2nd_max), maximum(c_2nd_max)
#         cmin == cmax_val ? (cmin - 0.1, cmax_val + 0.1) : (cmin, cmax_val)
#     else
#         (0.0, 1.0)
#     end

#     # --------------------------------------------------------------------------
#     # 4. Generate Subplots with Independent Colorbars
#     # --------------------------------------------------------------------------
#     fixed_title_str = join(["$k=$v" for (k, v) in fixed_variables], ", ")

#     p1 = heatmap(gammas, Es, gap_mat;
#         xlabel=L"\gamma", ylabel=L"E",
#         title="Localiser Gap - ($(fixed_title_str))",
#         colorbar_title=gap_title,
#         colorbar=true,
#         colormap=:plasma)

#     p2 = heatmap(gammas, Es, chern_mat;
#         xlabel=L"\gamma", ylabel=L"E",
#         title="Localiser Chern Number",
#         colorbar_title="Chern",
#         colorbar=true,
#         colormap=:RdBu)

#     # Panel 3: Dedicated panel for 2nd Max C(E) with its own colorbar
#     p3 = scatter(
#         xlabel=L"\gamma", ylabel=L"E",
#         title=L"2nd\ Max\ Cum.\ Chern\ C_{\mathrm{2nd}}(E)",
#         xlims=xlims, ylims=ylims
#     )

#     # Overlay C = 0 scatter points across panels
#     if !isempty(g_zero)
#         for p in (p1, p2, p3)
#             scatter!(p, g_zero, e_zero;
#                 label=L"C = %$(target_zero)",
#                 color=:white,
#                 markershape=:diamond,
#                 markersize=3.5,
#                 markerstrokewidth=0.6,
#                 markerstrokecolor=:black,
#                 alpha=0.9
#             )
#         end
#     end

#     # Overlay 2nd Max Cum Chern scatter points color-coded by c_2nd_max
#     if !isempty(g_2nd_max)
#         scatter!(p1, g_2nd_max, e_2nd_max;
#             marker_z=c_2nd_max,
#             label=L"C_{\mathrm{2nd\ max}}",
#             colormap=max_colormap,
#             clims=clims_2nd,
#             colorbar=false,
#             markershape=:circle,
#             markersize=4.0,
#             markerstrokewidth=0.4,
#             markerstrokecolor=:black,
#             alpha=0.9
#         )

#         scatter!(p2, g_2nd_max, e_2nd_max;
#             marker_z=c_2nd_max,
#             label=L"C_{\mathrm{2nd\ max}}",
#             colormap=max_colormap,
#             clims=clims_2nd,
#             colorbar=false,
#             markershape=:circle,
#             markersize=4.0,
#             markerstrokewidth=0.4,
#             markerstrokecolor=:black,
#             alpha=0.9
#         )

#         scatter!(p3, g_2nd_max, e_2nd_max;
#             marker_z=c_2nd_max,
#             label=L"C_{\mathrm{2nd\ max}}",
#             colorbar_title=L"C_{\mathrm{2nd\ max}}",
#             colormap=max_colormap,
#             clims=clims_2nd,
#             colorbar=true,
#             markershape=:circle,
#             markersize=4.0,
#             markerstrokewidth=0.4,
#             markerstrokecolor=:black,
#             alpha=0.9
#         )
#     end

#     plot!(p1, legend=:topright, xlims=xlims, ylims=ylims)
#     plot!(p2, legend=:topright, xlims=xlims, ylims=ylims)
#     plot!(p3, legend=:topright, xlims=xlims, ylims=ylims)

#     plt = plot(p1, p2, p3; layout=(1, 3), size=(2400, 600), margin=5Plots.mm)
    
#     isdir(dirname(filename)) || mkpath(dirname(filename))
#     savefig(plt, filename)

#     return plt
# end


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


for A in valid_Avals, B in valid_Bvals, m in valid_mvals
    
    # Extract (gamma, band_data) pairs for current slice
    sub_chern_pairs = [(item[4], item[5]) for item in crit_chern_valid_list if 
        isapprox(item[1], A; atol=1e-6) &&
        isapprox(item[2], B; atol=1e-6) &&
        isapprox(item[3], m; atol=1e-6)
    ]

    isempty(sub_chern_pairs) && continue

    fn = joinpath(output_path, "localiser_heatmap_A$(A)_B$(B)_m$(m).png")

    plt_specloc_gamma_vs_E_overlaid_crit_cum_chern(
        specloc_valid_df,
        sub_chern_pairs;
        logscale = true,
        target_p = 0.5,
        target_q = 0.0,
        tol = 0.005,
        xlims = (-3.0,3.0),
        ylims = (-6.0, 5.0),
        filename = fn,
        A = A, B = B, m = m
    )

    fn = joinpath(output_path, "localiser_rrrr_heatmap_A$(A)_B$(B)_m$(m).png")

    plt_specloc_gamma_vs_E_overlaid_crit_cum_chern_rrrr(
        specloc_valid_df,
        sub_chern_pairs;
        logscale = true,
        target_p = 0.5,
        target_q = 0.0,
        target_r=1.0,
        tol = 0.005,
        xlims = (-3.0,3.0),
        ylims = (-4.0, 4.0),
        filename = fn,
        A = A, B = B, m = m
    )

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
end
