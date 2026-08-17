using Dates
using ProgressMeter
using Distributed

# Force headless plotting by default (safe for CLI + parallel workers).
if !haskey(ENV, "GKSwstype")
    ENV["GKSwstype"] = "png"
end
if !haskey(ENV, "GKS_WSTYPE")
    ENV["GKS_WSTYPE"] = ENV["GKSwstype"]
end
if !haskey(ENV, "GKS_NO_GUI")
    ENV["GKS_NO_GUI"] = "1"
end

project_root = @__DIR__
include(joinpath(project_root, "plotting.jl"))
using .SpecLocPlotting

function newest_results_dir(root::String)
    dirs = filter(d -> isdir(d), readdir(root; join=true))
    isempty(dirs) && error("No result directories found under $(root)")
    sort!(dirs, by=d -> stat(d).mtime)
    return dirs[end]
end

function setup_parallel_workers(project_root::String, nworkers::Int)
    nworkers <= 1 && return

    needed = nworkers - nprocs()
    if needed > 0
        addprocs(needed)
    end

    plotting_file = joinpath(project_root, "plotting.jl")
    for p in workers()
        # Ensure each worker is also headless before loading Plots/GR.
        remotecall_wait(Core.eval, p, Main, :(ENV["GKSwstype"] = get(ENV, "GKSwstype", "png")))
        remotecall_wait(Core.eval, p, Main, :(ENV["GKS_WSTYPE"] = get(ENV, "GKS_WSTYPE", ENV["GKSwstype"])))
        remotecall_wait(Core.eval, p, Main, :(ENV["GKS_NO_GUI"] = get(ENV, "GKS_NO_GUI", "1")))
        remotecall_wait(include, p, plotting_file)
        remotecall_wait(Core.eval, p, Main, :(using .SpecLocPlotting))
    end
end

function main(; verbose::Bool=true, nworkers::Int=1)
    input_root = length(ARGS) >= 1 ? ARGS[1] : newest_results_dir(joinpath(project_root, "..", "data_collection", "hpc", "results"))
    output_root = length(ARGS) >= 2 ? ARGS[2] : joinpath(project_root, "plots")

    setup_parallel_workers(project_root, nworkers)

    if verbose
        println("Input root: $(input_root)")
        println("Output root: $(output_root)")
        println("Julia processes: $(nworkers)")
    end

    run_id = Dates.format(now(), "yyyymmdd_HHMMSS")
    index_file = process_all_cases(input_root, output_root; run_id=run_id, verbose=verbose, nworkers=nworkers)

    if verbose
        println("Done. Plot index file:")
        println(index_file)
    else
        println("Processing complete. Plot index: $(index_file)")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    verbose = lowercase(strip(get(ENV, "SPECLOC_VERBOSE", "true"))) != "false"
    nworkers = max(1, parse(Int, get(ENV, "SPECLOC_NWORKERS", "1")))
    main(; verbose=verbose, nworkers=nworkers)
end
