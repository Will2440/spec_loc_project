using Plots
using LaTeXStrings
using DataFrames
using JLD2
using Glob
using ProgressMeter


data_path = "/Users/Will/Documents/spec_loc_project/notebook_calcs/data/qwz_disorder_lowspecloc/low_specloc_spectrum_nEvals4_symmetric_disord_anderson_nAvgReals5_PARAMS_As1.0-1-1.0_Bs1.0-1.0_ms-1.0--1.0_gammas0.5-0.5_Ws0.0-1.0_Lxs50-50_Lys50-50_keptSquaretrue_x0y0sweepfalse_Es0.0-3.0_kappas0.02-0.02"

output_path = "/Users/Will/Documents/spec_loc_project/notebook_calcs/results/qwz_disorder_lowspecloc/low_specloc_spectrum_nEvals4_symmetric_disord_anderson_nAvgReals5_PARAMS_As1.0-1-1.0_Bs1.0-1.0_ms-1.0--1.0_gammas0.5-0.5_Ws0.0-1.0_Lxs50-50_Lys50-50_keptSquaretrue_x0y0sweepfalse_Es0.0-3.0_kappas0.02-0.02"
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


# # compatible with btoh specloc calcs (slow and fast) where specloc gap needs to be inferred form low-lying spectrum
function plt_heatmap_localiser_W_vs_E(
    df::DataFrame; 
    atol::Real=1e-8,
    logscale::Bool=false,
    filename::String="plots/localiser_W_E_heatmap.png",
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
    Ws = sort(unique(subdf.W))
    Es = sort(unique(subdf.E))
    nW, nE = length(Ws), length(Es)

    gap_mat = fill(NaN, nE, nW)

    for row in eachrow(subdf)
        xi = findfirst(==(row.W), Ws)
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
    
    # Return early if no localiser data found
    if isempty(valid_gaps)
        @warn "No localiser data found for parameters: $fixed_variables"
        return nothing
    end
    
    c_min, c_max = minimum(valid_gaps), maximum(valid_gaps)
    if c_min == c_max
        c_min -= 1e-5
        c_max += 1e-5
    end

    fixed_title_str = join(["$k=$v" for (k, v) in fixed_variables], ", ")

    p1 = heatmap(Ws, Es, gap_mat;
        clims = (c_min, c_max),
        xlabel = L"W", 
        ylabel = L"E",
        title = "Localiser Gap - ($fixed_title_str)",
        colorbar_title = gap_title,
        colormap = :plasma,
        aspect_ratio = :auto,
        colorbar = true)

    if !isempty(filename)
        isdir(dirname(filename)) || mkpath(dirname(filename))
        savefig(p1, filename)
    end

    return p1
end



# Plot W vs E for each kappa value
plt_gammas = [0.0, 0.5]
L_value = df.Lx[1]  # Assuming Lx is constant across the DataFrame
for kappa in kappavals
    for gamma in plt_gammas
        plt_heatmap_localiser_W_vs_E(df; 
            filename=joinpath(output_path, "localiser_W_E_heatmap_kappa$(kappa)_gamma$(gamma).png"),
            logscale=true,
            kappa=kappa,
            gamma=gamma,
            Lx=L_value
        )
    end
end
