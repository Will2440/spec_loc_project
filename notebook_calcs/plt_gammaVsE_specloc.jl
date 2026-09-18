using Plots
using LaTeXStrings
using DataFrames
using JLD2
using Glob
using ProgressMeter


data_path = "/Users/Will/Documents/spec_loc_project/notebook_calcs/data/qwz_disorder_lowspecloc/low_specloc_spectrum_nEvals4_asym_sin_sum_disord_none_nAvgReals1_PARAMS_As1.0-1-1.0_Bs1.0-1.0_ms-1.0--1.0_gammas-3.0-3.0_IndTypasym_sin_sum_Ws0.0-0.0_Lxs20-20_Lys20-20_keptSquaretrue_x0y0sweepfalse_Es-3.0-3.0_kappas0.2-0.2_embedPhis0.0-0.0_embedDs0.0-5.0"

output_path = "/Users/Will/Documents/spec_loc_project/notebook_calcs/results/qwz_disorder_lowspecloc/low_specloc_spectrum_nEvals4_asym_sin_sum_disord_none_nAvgReals1_PARAMS_As1.0-1-1.0_Bs1.0-1.0_ms-1.0--1.0_gammas-3.0-3.0_IndTypasym_sin_sum_Ws0.0-0.0_Lxs20-20_Lys20-20_keptSquaretrue_x0y0sweepfalse_Es-3.0-3.0_kappas0.2-0.2_embedPhis0.0-0.0_embedDs0.0-5.0"
isdir(output_path) || mkpath(output_path)

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

df = process_jld2_files_fluid(data_path)

Avals = sort(unique(df.A))
Bvals = sort(unique(df.B))
mvals = sort(unique(df.m))
gammavals = sort(unique(df.gamma))
Wvals = sort(unique(df.W))
kappavals = sort(unique(df.kappa))
Evals = sort(unique(df.E))
embedPhiVals = sort(unique(df.phi))
embedDVals = sort(unique(df.d))

println("Unique parameter values: | Type | min-length-max")
println("A: ", Avals, " | ", typeof(Avals[1]), " | ", minimum(Avals), " - ", length(Avals), " - ", maximum(Avals))
println("B: ", Bvals, " | ", typeof(Bvals[1]), " | ", minimum(Bvals), " - ", length(Bvals), " - ", maximum(Bvals))
println("m: ", mvals, " | ", typeof(mvals[1]), " | ", minimum(mvals), " - ", length(mvals), " - ", maximum(mvals))
println("gamma: ", gammavals, " | ", typeof(gammavals[1]), " | ", minimum(gammavals), " - ", length(gammavals), " - ", maximum(gammavals))
println("W: ", Wvals, " | ", typeof(Wvals[1]), " | ", minimum(Wvals), " - ", length(Wvals), " - ", maximum(Wvals))
println("kappa: ", kappavals, " | ", typeof(kappavals[1]), " | ", minimum(kappavals), " - ", length(kappavals), " - ", maximum(kappavals))
println("E: ", Evals, " | ", typeof(Evals[1]), " | ", minimum(Evals), " - ", length(Evals), " - ", maximum(Evals))
println("embedPhi: ", embedPhiVals, " | ", typeof(embedPhiVals[1]), " | ", minimum(embedPhiVals), " - ", length(embedPhiVals), " - ", maximum(embedPhiVals))
println("embedD: ", embedDVals, " | ", typeof(embedDVals[1]), " | ", minimum(embedDVals), " - ", length(embedDVals), " - ", maximum(embedDVals))

# # compatible with btoh specloc calcs (slow and fast) where specloc gap needs to be inferred form low-lying spectrum
function plt_heatmap_localiser_gamma_vs_E(
    df::DataFrame,
    gamma_data_pairs::Union{Vector{<:Tuple{Real, Any}}, Nothing}=nothing; 
    atol::Real=1e-8,
    logscale::Bool=false,
    target_p::Float64=0.5,
    target_q::Float64=0.0,
    tol::Float64=0.01,
    filename::String="plots/localiser_gamma_E_heatmap.png",
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
            chern_mat[yi, xi] = if hasproperty(row, :chern_marker) && !ismissing(row.chern_marker)
                Float64(row.chern_marker)
            elseif hasproperty(row, :signature) && !ismissing(row.signature)
                Float64(row.signature) / 2.0  # Convert signature to chern marker
            elseif hasproperty(row, :chern_number) && !ismissing(row.chern_number)
                Float64(row.chern_number)
            else
                NaN
            end
        end
    end

    # Log-scale transformation logic
    if logscale
        gap_mat = log10.(max.(gap_mat, 1e-12))
        gap_title = L"\log_{10}(\mathrm{min}|\lambda|)"
    else
        gap_title = L"\mathrm{min}|\lambda|"
    end

    # Explicit clims safety bounds
    valid_gaps = filter(!isnan, vec(gap_mat))
    c_min, c_max = isempty(valid_gaps) ? (0.0, 1.0) : (minimum(valid_gaps), maximum(valid_gaps))
    if c_min == c_max
        c_min -= 1e-5
        c_max += 1e-5
    end

    fixed_title_str = join(["$k=$v" for (k, v) in fixed_variables], ", ")

    p1 = heatmap(gammas, Es, gap_mat;
        clims = (c_min, c_max),
        xlabel = L"\gamma", 
        ylabel = L"E",
        title = "Localiser Gap - ($fixed_title_str)",
        colorbar_title = gap_title,
        colormap = :plasma,
        aspect_ratio = :auto,
        colorbar = true,
        legend = false)
    
    p2 = heatmap(gammas, Es, chern_mat;
        xlabel = L"\gamma", 
        ylabel = L"E",
        title = "Localiser Chern Number",
        colorbar_title = L"\mathrm{Chern}",
        colormap = :plasma,
        aspect_ratio = :auto,
        colorbar = true,
        legend = :outerright)

    # --------------------------------------------------------------------------
    # 3. Extract & Overlay Scatter Contour Points (If provided)
    # --------------------------------------------------------------------------
    if !isnothing(gamma_data_pairs)
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
            cum_1 = band_data.cum_chern_per_band[1, :, :][cite: 1]
            E_1   = band_data.plaquette_energies[1, :, :][cite: 1]
            cum_2 = band_data.cum_chern_per_band[2, :, :][cite: 1]
            E_2   = band_data.plaquette_energies[2, :, :][cite: 1]

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

        for (g_pts, e_pts, lbl, clr, shp) in series_configs
            if !isempty(g_pts)
                scatter!(p1, g_pts, e_pts; label=false, color=clr, markershape=shp, markersize=3.5, markerstrokewidth=0.4, markerstrokecolor=:black, alpha=0.85)
                scatter!(p2, g_pts, e_pts; label=lbl, color=clr, markershape=shp, markersize=3.5, markerstrokewidth=0.4, markerstrokecolor=:black, alpha=0.85)
            end
        end
    end

    plt = plot(p1, p2; layout=(1, 2), size=(1500, 480))

    if !isempty(filename)
        isdir(dirname(filename)) || mkpath(dirname(filename))
        savefig(plt, filename)
    end

    return plt
end

function plt_localiser_spectrum(
    df::DataFrame,
    x_var::Symbol;
    atol::Real=1e-8,
    filename::String="plots/localiser_spectrum.png",
    log_dos::Bool=false,
    fixed_variables...
)
    
    # Filter by parameters
    col_names_str = string.(names(df))
    string(x_var) in col_names_str || error("x_var $(x_var) is not a DataFrame column.")
    
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

    # Sort by x_var for plotting
    subdf = sort(subdf, x_var)
    
    # Extract spectra and flatten into plottable arrays
    spectrum_data = []
    x_plot = []
    
    for row in eachrow(subdf)
        spec = row.spectrum  # Vector of Complex eigenvalues
        x_val = row[x_var]
        
        # Convert to real if Complex
        if eltype(spec) <: Complex
            spec = real.(spec)
        end
        
        # Apply log scaling if requested
        if log_dos
            spec = log10.(abs.(spec))
        end
        
        # Add each eigenvalue with corresponding x value
        for eig in spec
            push!(spectrum_data, eig)
            push!(x_plot, x_val)
        end
    end
    
    ylabel_str = log_dos ? "log10(|Localiser eigenvalue|)" : "Localiser eigenvalue"
    
    p = scatter(x_plot, spectrum_data;
        xlabel=string(x_var),
        ylabel=ylabel_str,
        title="Localiser Spectrum vs $(x_var) - $(fixed_variables)",
        markersize=3,
        alpha=0.6,
        legend=false)
    
    isdir(dirname(filename)) || mkpath(dirname(filename))
    savefig(p, filename)
    
end


for (kappa, W, embedPhi, embedD, m) in Iterators.product(kappavals, Wvals, embedPhiVals, embedDVals, mvals)
    plt_heatmap_localiser_gamma_vs_E(df; 
        filename=joinpath(output_path, "localiser_gamma_E_heatmap_kappa$(kappa)_W$(W)_phi$(embedPhi)_d$(embedD)_m$(m).png"),
        logscale=true,
        kappa=kappa,
        W=W,
        phi=embedPhi,
        d=embedD,
        m=m,
    )
end
