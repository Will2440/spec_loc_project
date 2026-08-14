using Plots
using LaTeXStrings
using DataFrames
using JLD2
using Glob
using ProgressMeter


data_path = "/Users/Will/Documents/spec_loc_project/data/qwz_model"

output_path = "/Users/Will/Documents/spec_loc_project/results/qwz_model"
isdir(output_path) || mkpath(output_path)

function process_jld2_files_fluid(
    folder_path::String
)
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
kappavals = sort(unique(df.kappa))
Evals = sort(unique(df.E))


function plt_localiser_position_heatmaps_sparse(
    df::DataFrame; 
    atol::Real=1e-8,
    logscale::Bool=false,
    filename::String="plots/localiser_position_heatmaps_sparse.png",
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

    # Extract unique x, y positions (as they appear in the data)
    xs = sort(unique(subdf.x))
    ys = sort(unique(subdf.y))
    nx, ny = length(xs), length(ys)
    
    # Create lookup dicts
    x_idx = Dict(xs[i] => i for i in 1:nx)
    y_idx = Dict(ys[i] => i for i in 1:ny)

    # Initialize matrices for the actual data points only
    gap_mat = fill(NaN, ny, nx)  # (y, x) for heatmap orientation
    # sig_mat = fill(NaN, ny, nx)
    chern_mat = fill(NaN, ny, nx)

    for row in eachrow(subdf)
        xi = x_idx[row.x]
        yi = y_idx[row.y]
        gap_mat[yi, xi] = row.localiser_gap
        # sig_mat[yi, xi] = row.signature
        chern_mat[yi, xi] = ismissing(row.chern_number) ? NaN : Float64(row.chern_number)
    end

    if logscale
        gap_mat = log10.(gap_mat)
    end

    p1 = heatmap(xs, ys, gap_mat;
        xlabel="x", ylabel="y",
        title="Localiser min|λ| - $(fixed_variables)",
        colorbar_title="min|λ|",
        aspect_ratio=:auto)
    
    p3 = heatmap(xs, ys, chern_mat;
        xlabel="x", ylabel="y",
        title="Localiser-derived Chern",
        colorbar_title="chern",
        # clims=(-2, 2),
        aspect_ratio=:auto)

    plt = plot(p1, p3; layout=(1, 2), size=(1400, 420))
    isdir(dirname(filename)) || mkpath(dirname(filename))
    savefig(plt, filename)
end


## plot for
A_plt = Avals
B_plt = Bvals
m_plt = mvals
E_plt = Evals
kappa_plt = kappavals
gamma_plt = gammavals


for A_val in A_plt, B_val in B_plt, m_val in m_plt, E_val in E_plt, kappa_val in kappa_plt, gamma_val in gamma_plt

    foldername = joinpath("plots", "perturbed_spatial_scan_larger", "m$(m_val)_gamma$(gamma_val)")
    isdir(foldername) || mkpath(foldername)

    plt_localiser_position_heatmaps_sparse(
        df; 
        A=A_val, 
        B=B_val, 
        m=m_val, 
        E=E_val, 
        kappa=kappa_val, 
        gamma=gamma_val,
        filename=joinpath(foldername, "localiser_position_heatmaps_A$(A_val)_B$(B_val)_m$(m_val)_E$(E_val)_kappa$(kappa_val)_gamma$(gamma_val).png")
    )
end
