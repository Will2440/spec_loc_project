using Plots
using LaTeXStrings
using DataFrames
using JLD2
using Glob
using ProgressMeter


data_path = "/Users/Will/Documents/spec_loc_project/notebook_calcs/data/qwZ_fast_low_spec/clean_low_specloc_spectrum_nEvals10_symmetric_disord_none_PARAMS_As1.0-1-1.0_Bs1.0-1.0_ms-1.0--1.0_gammas3.0-1.0_Ws0.0-0.0_Lxs50-50_Lys50-50_keptSquaretrue_x0y0sweepfalse_Es1.5-0.5_kappas0.001-1.0"

output_path = "/Users/Will/Documents/spec_loc_project/notebook_calcs/results/qwZ_fast_low_spec/clean_low_specloc_spectrum_nEvals10_symmetric_disord_none_PARAMS_As1.0-1-1.0_Bs1.0-1.0_ms-1.0--1.0_gammas3.0-1.0_Ws0.0-0.0_Lxs50-50_Lys50-50_keptSquaretrue_x0y0sweepfalse_Es1.5-0.5_kappas0.001-1.0"
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
Lxvals = sort(unique(df.Lx))
Lyvals = sort(unique(df.Ly))


function plt_specloc_spectrum_vs_E(
    df::DataFrame, 
    A::Real, B::Real, m::Real, gamma::Real, W::Real, 
    Lx::Int, Ly::Int, kappa::Real;
    save_path::String="",
    log_spec::Bool=false
)
    # Filter the DataFrame for specified parameters
    filtered_df = filter(row -> row.A == A && row.B == B && row.m == m && row.gamma == gamma &&
                                row.W == W && row.Lx == Lx && row.Ly == Ly &&
                                row.kappa == kappa,
                         df)

    if nrow(filtered_df) == 0
        @warn "No data found for specified parameters: A=$A, B=$B, m=$m, γ=$gamma, W=$W, Lx=$Lx, Ly=$Ly, κ=$kappa"
        return nothing
    end

    # Sort strictly by energy E
    sort!(filtered_df, :E)

    # Extract eigenvalue vectors and sort each vector by value (spectrum order)
    eval_vecs = [sort(real.(vec)) for vec in filtered_df.low_lying_evals]

    n_E = nrow(filtered_df)
    max_k = maximum(length, eval_vecs)
    evals_mat = fill(NaN, n_E, max_k)

    for i in 1:n_E
        k_i = length(eval_vecs[i])
        evals_mat[i, 1:k_i] .= eval_vecs[i]
    end

    if log_spec
        evals_mat = log10.(abs.(evals_mat) .+ eps())
    end

    label = log_spec ? "log10(|Spec[L]|)" : "Spec[L]"

    # plt = plot(filtered_df.E, evals_mat,
    #            xlabel="Energy (E)", 
    #            ylabel=label,
    #            title="Spec[L] vs E (Lx=$Lx, Ly=$Ly, γ=$gamma)",
    #            lc=:black,
    #            label="",
    #            legend=false
    # )

    plt = scatter(
            filtered_df.E, evals_mat,
            xlabel="Energy (E)", 
            ylabel=label,
            title="Spec[L] vs E (Lx=$Lx, Ly=$Ly, γ=$gamma, W=$W)",
            mc=:black,            # Marker color
            shape=:circle,        # Circle shape
            ms=1.5,               # Small marker size
            msw=0,                # No marker stroke/outline width
            msa=0,                # No marker stroke alpha
            label="",
            legend=false
    )

    hline!(plt, [0.0], color=:red, linestyle=:dash, label="E=0")
        
    if save_path != ""
        isdir(dirname(save_path)) || mkpath(dirname(save_path))
        savefig(plt, save_path)
        println("Plot saved to $save_path")
    else
        display(plt)
    end

    return plt
end

function plt_localiser_gap_xy_heatmaps_vs_E(
    df::DataFrame;
    save_dir::String="",
    logscale::Bool=false,
    atol::Real=1e-8,
    x_col::Symbol = hasproperty(df, :x0) ? :x0 : :x,
    y_col::Symbol = hasproperty(df, :y0) ? :y0 : :y,
    fixed_variables...
)
    # 1. Filter DataFrame by fixed parameter values
    subdf = df
    for (k, v) in fixed_variables
        string(k) in string.(names(df)) || error("Filter key $(k) is not a DataFrame column.")
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
        return Plots.Plot[]
    end

    xs = sort(unique(subdf[!, x_col]))
    ys = sort(unique(subdf[!, y_col]))
    Es = sort(unique(subdf.E))

    nx, ny = length(xs), length(ys)
    plots_vector = Plots.Plot[]
    fixed_title_str = join(["$k=$v" for (k, v) in fixed_variables], ", ")

    # 2. Iterate over energy levels E
    for E_val in Es
        e_mask = abs.(Float64.(subdf.E) .- Float64(E_val)) .<= Float64(atol)
        e_subdf = subdf[e_mask, :]

        # Warn if non-spatial parameters remain unfixed (causing overwrites)
        if nrow(e_subdf) > nx * ny
            @warn "e_subdf contains redundant rows ($(nrow(e_subdf)) rows for $nx×$ny grid). Ensure all parameter columns are passed in fixed_variables."
        end

        gap_mat = fill(NaN, ny, nx)

        for row in eachrow(e_subdf)
            xi = findfirst(==(row[x_col]), xs)
            yi = findfirst(==(row[y_col]), ys)

            gap_val = if hasproperty(row, :localiser_gap) && !ismissing(row.localiser_gap)
                Float64(row.localiser_gap)
            elseif hasproperty(row, :low_lying_evals) && !ismissing(row.low_lying_evals)
                # Compute magnitude directly on complex evaluations
                minimum(abs.(vec(row.low_lying_evals)))
            else
                NaN
            end

            if !isnothing(xi) && !isnothing(yi)
                gap_mat[yi, xi] = gap_val
            end
        end

        # Handle log scale transformation with floating-point floor to avoid -Inf
        if logscale
            gap_mat = log10.(max.(gap_mat, 1e-12))
            cb_title = L"\log_{10}(\text{Gap})"
        else
            cb_title = L"\text{Gap}"
        end

        # Force colormap limits based on current matrix min/max
        valid_vals = filter(!isnan, vec(gap_mat))
        c_min, c_max = isempty(valid_vals) ? (0.0, 1.0) : (minimum(valid_vals), maximum(valid_vals))
        if c_min == c_max
            c_min -= 1e-5
            c_max += 1e-5
        end

        plt = heatmap(
            xs, ys, gap_mat;
            clims = (c_min, c_max),
            xlabel = L"x_0",
            ylabel = L"y_0",
            title = "Localiser Gap at E = $(round(E_val, digits=4))\n($fixed_title_str)",
            colorbar_title = cb_title,
            colorbar = true,
            colormap = :plasma,
            aspect_ratio = :equal
        )

        push!(plots_vector, plt)

        if !isempty(save_dir)
            isdir(save_dir) || mkpath(save_dir)
            filename = joinpath(save_dir, "localiser_gap_xy_E$(round(E_val, digits=4)).png")
            savefig(plt, filename)
            println("Saved: $filename")
        end
    end

    return plots_vector
end

function plt_localiser_gap_E_vs_kappa_heatmap(
    df::DataFrame;
    save_path::String="",
    logscale::Bool=false,
    atol::Real=1e-8,
    fixed_variables...
)
    # 1. Filter DataFrame by fixed parameter values
    subdf = df
    for (k, v) in fixed_variables
        string(k) in string.(names(df)) || error("Filter key $(k) is not a DataFrame column.")
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

    # 2. Extract unique grid coordinates
    kappas = sort(unique(subdf.kappa))
    Es = sort(unique(subdf.E))
    nk, nE = length(kappas), length(Es)

    gap_mat = fill(NaN, nE, nk)

    # 3. Populate matrix from subdf rows
    for row in eachrow(subdf)
        ki = findfirst(==(row.kappa), kappas)
        ei = findfirst(==(row.E), Es)

        gap_val = if hasproperty(row, :localiser_gap) && !ismissing(row.localiser_gap)
            Float64(row.localiser_gap)
        elseif hasproperty(row, :low_lying_evals) && !ismissing(row.low_lying_evals)
            # Compute gap directly from low-lying eigenvalues magnitude
            minimum(abs.(vec(row.low_lying_evals)))
        else
            NaN
        end

        if !isnothing(ki) && !isnothing(ei)
            gap_mat[ei, ki] = gap_val
        end
    end

    # 4. Log scale & colorbar title logic
    if logscale
        gap_mat = log10.(max.(gap_mat, 1e-12))
        cb_title = L"\log_{10}(\text{Gap})"
    else
        cb_title = L"\text{Gap}"
    end

    # 5. Robust colormap limits
    valid_vals = filter(!isnan, vec(gap_mat))
    c_min, c_max = isempty(valid_vals) ? (0.0, 1.0) : (minimum(valid_vals), maximum(valid_vals))
    if c_min == c_max
        c_min -= 1e-5
        c_max += 1e-5
    end

    # 6. Generate plot
    fixed_title_str = join(["$k=$v" for (k, v) in fixed_variables], ", ")
    plt = heatmap(
        kappas, Es, gap_mat;
        clims = (c_min, c_max),
        xlabel = L"\kappa",
        ylabel = L"E",
        title = "Localiser Gap vs E and κ\n($fixed_title_str)",
        colorbar_title = cb_title,
        colorbar = true,
        colormap = :plasma,
        aspect_ratio = :auto
    )

    if !isempty(save_path)
        isdir(dirname(save_path)) || mkpath(dirname(save_path))
        savefig(plt, save_path)
        println("Plot saved to $save_path")
    end

    return plt
end


for (A, B, m, gamma, W, Lx, Ly, kappa) in Iterators.product(Avals, Bvals, mvals, gammavals, Wvals, Lxvals, Lyvals, kappavals)
    logscale = true
    
    # 1. Plot spectrum vs Energy across spatial configurations
    save_path = joinpath(output_path, "specloc_spectrum_log$(logscale)_A$(A)_B$(B)_m$(m)_gamma$(gamma)_W$(W)_Lx$(Lx)_Ly$(Ly)_kappa$(kappa).png")
    plt_specloc_spectrum_vs_E(
        df, A, B, m, gamma, W, Lx, Ly, kappa; 
        save_path = save_path,
        log_spec = logscale
    )

    # # 2. Plot spatial heatmaps per energy slice
    # heatmaps_dir = joinpath(output_path, "spatial_gap_heatmaps_log$(logscale)_A$(A)_B$(B)_m$(m)_gamma$(gamma)_W$(W)_Lx$(Lx)_Ly$(Ly)_kappa$(kappa)")
    # plt_localiser_gap_xy_heatmaps_vs_E(
    #     df;
    #     A = A,
    #     B = B,
    #     m = m,
    #     gamma = gamma,
    #     W = W,
    #     Lx = Lx,
    #     Ly = Ly,
    #     kappa = kappa,
    #     logscale = logscale,
    #     save_dir = heatmaps_dir
    # )

    # 3. Plot heatmap of localiser gap vs Energy and kappa
    save_path_kappa = joinpath(output_path, "localiser_gap_E_vs_kappa_log$(logscale)_A$(A)_B$(B)_m$(m)_gamma$(gamma)_W$(W)_Lx$(Lx)_Ly$(Ly).png")
    plt_localiser_gap_E_vs_kappa_heatmap(
        df;
        A = A,
        B = B,
        m = m,
        gamma = gamma,
        W = W,
        Lx = Lx,
        Ly = Ly,
        logscale = logscale,
        save_path = save_path_kappa
    )
end