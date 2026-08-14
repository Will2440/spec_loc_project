using Plots
using LaTeXStrings
using DataFrames
using JLD2
using Glob
using ProgressMeter


data_path = "/Users/Will/Documents/spec_loc_project/notebook_calcs/data/qwz_model_pert_disordered"

# output_path = "/Users/Will/Documents/spec_loc_project/notebook_calcs/plots/qwz_model_pert_disordered"
# isdir(output_path) || mkpath(output_path)

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


function plt_localiser_gamma_W_heatmap(
    df::DataFrame; 
    atol::Real=1e-8,
    logscale::Bool=false,
    filename::String="plots/localiser_gamma_W_heatmap.png",
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

    # Extract unique gamma, W positions (as they appear in the data)
    gammas = sort(unique(subdf.gamma))
    Ws = sort(unique(subdf.W))
    ng, nW = length(gammas), length(Ws)

    # plotting heatmap where W us x axis and gamm ais y axis and the colouir values are (1) localiser gap and (2) chern number
    gap_mat = fill(NaN, ng, nW)  # (gamma, W) for heatmap orientation
    chern_mat = fill(NaN, ng, nW)  # (gamma, W) for heatmap orientation

    for row in eachrow(subdf)
        xi = findfirst(==(row.W), Ws)
        yi = findfirst(==(row.gamma), gammas)
        gap_mat[yi, xi] = row.localiser_gap
        # sig_mat[yi, xi] = row.signature
        chern_mat[yi, xi] = ismissing(row.chern_number) ? NaN : Float64(row.chern_number)
    end

    if logscale
        gap_mat = log10.(gap_mat)
        gap_title = "log10(min|λ|)"
    else
        gap_title = "min|λ|"
    end

    p1 = heatmap(Ws, gammas, gap_mat;
        xlabel="W", ylabel="gamma",
        title="Localiser min|λ| - $(fixed_variables)",
        colorbar_title=gap_title,
        aspect_ratio=:auto)
    
    p3 = heatmap(Ws, gammas, chern_mat;
        xlabel="W", ylabel="gamma",
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
W_plt = Wvals


for A_val in A_plt, B_val in B_plt, m_val in m_plt, E_val in E_plt, kappa_val in kappa_plt

    foldername = joinpath("plots", "perturbed_disordered_gammaVsW", "m$(m_val)")
    isdir(foldername) || mkpath(foldername)

    plt_localiser_gamma_W_heatmap(
        df; 
        A=A_val, 
        B=B_val, 
        m=m_val, 
        E=E_val, 
        kappa=kappa_val, 
        logscale=true,
        filename=joinpath(foldername, "localiser_gamma_W_heatmap_A$(A_val)_B$(B_val)_m$(m_val)_E$(E_val)_kappa$(kappa_val).png")
    )
end
