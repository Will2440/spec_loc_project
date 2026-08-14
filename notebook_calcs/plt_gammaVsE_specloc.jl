using Plots
using LaTeXStrings
using DataFrames
using JLD2
using Glob
using ProgressMeter


data_path = "/Users/Will/Documents/spec_loc_project/notebook_calcs/data/test_run"

output_path = "/Users/Will/Documents/spec_loc_project/notebook_calcs/results/test_run"
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

function pltheatmap_localiser_gamma_vs_E(
    df::DataFrame; 
    atol::Real=1e-8,
    logscale::Bool=false,
    filename::String="plots/localiser_gamma_E_heatmap.png",
    fixed_variables...
)
    
    # Filter by parameters
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

    # Extract unique gamma, E positions (as they appear in the data)
    gammas = sort(unique(subdf.gamma))
    Es = sort(unique(subdf.E))
    ngamma, nE = length(gammas), length(Es)

    # plotting heatmap where gamma is x axis and E is y axis and the colour values are (1) localiser gap and (2) chern number
    gap_mat = fill(NaN, nE, ngamma)  # (E, gamma) for heatmap orientation
    chern_mat = fill(NaN, nE, ngamma)  # (E, gamma) for heatmap orientation

    for row in eachrow(subdf)
        xi = findfirst(==(row.gamma), gammas)
        yi = findfirst(==(row.E), Es)
        gap_mat[yi, xi] = row.localiser_gap
        chern_mat[yi, xi] = ismissing(row.chern_number) ? NaN : Float64(row.chern_number)
    end

    if logscale
        gap_mat = log10.(gap_mat)
        gap_title = "log10(min|λ|)"
    else
        gap_title = "min|λ|"
    end

    p1 = heatmap(gammas, Es, gap_mat;
        xlabel="gamma", ylabel="E",
        title="Localiser min|λ| - $(fixed_variables)",
        colorbar_title=gap_title,
        aspect_ratio=:auto)
    
    p2 = heatmap(gammas, Es, chern_mat;
        xlabel="gamma", ylabel="E",
        title="Localiser-derived Chern",
        colorbar_title="chern",
        aspect_ratio=:auto)

    plt = plot(p1, p2; layout=(1, 2), size=(1400, 420))
    isdir(dirname(filename)) || mkpath(dirname(filename))
    savefig(plt, filename)

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


pltheatmap_localiser_gamma_vs_E(df; 
    filename=joinpath(output_path, "localiser_gamma_E_heatmap.png"),
    logscale=true
)