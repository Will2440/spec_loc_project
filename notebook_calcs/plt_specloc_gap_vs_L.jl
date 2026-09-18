using Plots
using LaTeXStrings
using DataFrames
using JLD2
using Glob
using ProgressMeter
using Printf


data_path = "/Users/Will/Documents/spec_loc_project/notebook_calcs/data/qwz_disorder_lowspecloc/low_specloc_spectrum_nEvals4_symmetric_disord_none_nAvgReals1_PARAMS_As1.0-1-1.0_Bs1.0-1.0_ms-1.0--1.0_gammas0.5-0.5_Ws0.0-0.0_Lxs10-200_Lys10-200_keptSquaretrue_x0y0sweepfalse_Es1.3-1.3_kappas0.02-0.02_embedPhis0.0-0.0_embedDs0.0-0.0"

output_path = "/Users/Will/Documents/spec_loc_project/notebook_calcs/results/qwz_disorder_lowspecloc/low_specloc_spectrum_nEvals4_symmetric_disord_none_nAvgReals1_PARAMS_As1.0-1-1.0_Bs1.0-1.0_ms-1.0--1.0_gammas0.5-0.5_Ws0.0-0.0_Lxs10-200_Lys10-200_keptSquaretrue_x0y0sweepfalse_Es1.3-1.3_kappas0.02-0.02_embedPhis0.0-0.0_embedDs0.0-0.0"
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

# ============================================================================
# Function to plot spectral localization gap vs L
# ============================================================================
function plt_specloc_gap_vs_L(
    df::DataFrame;
    atol::Real=1e-8,
    logscale::Bool=true,
    filename::String="plots/specloc_gap_vs_L.png",
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

    # Sort by Lx for consistent plotting
    subdf = sort(subdf, :Lx)

    # --------------------------------------------------------------------------
    # 2. Extract Gap vs L Data
    # --------------------------------------------------------------------------
    Lvals = Float64[]
    gaps = Float64[]
    kappavals_plot = Float64[]

    for row in eachrow(subdf)
        # Extract gap from localiser_gap or infer from low_lying_evals
        gap_val = if hasproperty(row, :localiser_gap) && !ismissing(row.localiser_gap)
            Float64(row.localiser_gap)
        elseif hasproperty(row, :low_lying_evals) && !ismissing(row.low_lying_evals)
            minimum(abs.(vec(row.low_lying_evals)))
        else
            NaN
        end

        if !isnan(gap_val)
            push!(Lvals, Float64(row.Lx))
            push!(gaps, gap_val)
            # Extract kappa value for annotation
            kappa_val = hasproperty(row, :kappa) && !ismissing(row.kappa) ? Float64(row.kappa) : NaN
            push!(kappavals_plot, kappa_val)
        end
    end

    if isempty(gaps)
        @warn "No valid gap data found after filtering"
        return nothing
    end

    # --------------------------------------------------------------------------
    # 3. Apply log scaling if requested
    # --------------------------------------------------------------------------
    if logscale
        gaps_plot = log10.(max.(gaps, 1e-12))
        ylabel_str = L"\log_{10}(\text{min}|\lambda|)"
    else
        gaps_plot = gaps
        ylabel_str = L"\text{Spectral Localization Gap (min}|\lambda|)"
    end

    fixed_title_str = join(["$k=$v" for (k, v) in fixed_variables], ", ")

    # --------------------------------------------------------------------------
    # 4. Create Plot
    # --------------------------------------------------------------------------
    p = plot(Lvals, gaps_plot;
        xlabel=L"L_x",
        ylabel=ylabel_str,
        title="Spectral Localization Gap vs L - ($fixed_title_str)",
        linewidth=2.5,
        marker=:circle,
        markersize=6,
        markerstrokewidth=1,
        legend=false,
        grid=true,
        gridstyle=:dash,
        gridalpha=0.3,
        size=(800, 550))

    # --------------------------------------------------------------------------
    # 5. Annotate each point with kappa value below the marker
    # --------------------------------------------------------------------------
    for (i, (x_val, y_val, kappa_val)) in enumerate(zip(Lvals, gaps_plot, kappavals_plot))
        if !isnan(kappa_val)
            # Format kappa value for annotation (round to 4 significant figures)
            kappa_str = @sprintf("%.3g", kappa_val)
            annotate!(p, x_val, y_val - 0.15 * (maximum(gaps_plot) - minimum(gaps_plot)), 
                     text(kappa_str, 7, :center, :top))
        end
    end

    if !isempty(filename)
        isdir(dirname(filename)) || mkpath(dirname(filename))
        savefig(p, filename)
        println("Plot saved to: $filename")
    end

    return p
end

# ============================================================================
# Call the plotting function for each fixed parameter combination
# ============================================================================
for (A, B, m, W, embedPhi, embedD) in Iterators.product(
    Avals, Bvals, mvals, Wvals, embedPhiVals, embedDVals
)
    # Plot for each unique (E, gamma) pair in the fixed parameter set
    subdf_params = filter(row ->
        abs(row.A - A) < 1e-8 &&
        abs(row.B - B) < 1e-8 &&
        abs(row.m - m) < 1e-8 &&
        abs(row.W - W) < 1e-8 &&
        abs(row.phi - embedPhi) < 1e-8 &&
        abs(row.d - embedD) < 1e-8,
        df
    )

    # Get unique E and gamma values for this parameter set
    unique_E = sort(unique(subdf_params.E))
    unique_gamma = sort(unique(subdf_params.gamma))

    for E_val in unique_E, gamma_val in unique_gamma
        plt_specloc_gap_vs_L(df;
            filename=joinpath(output_path, "specloc_gap_vs_L_E$(E_val)_gamma$(gamma_val)_A$(A)_B$(B)_m$(m)_W$(W)_phi$(embedPhi)_d$(embedD).png"),
            logscale=true,
            A=A, B=B, m=m, E=E_val, gamma=gamma_val, W=W, phi=embedPhi, d=embedD
        )
    end
end


