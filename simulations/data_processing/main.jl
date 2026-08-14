using Dates

project_root = @__DIR__
include(joinpath(project_root, "plotting.jl"))
using .SpecLocPlotting

function newest_results_dir(root::String)
    dirs = filter(d -> isdir(d), readdir(root; join=true))
    isempty(dirs) && error("No result directories found under $(root)")
    sort!(dirs, by=d -> stat(d).mtime)
    return dirs[end]
end

function main()
    input_root = length(ARGS) >= 1 ? ARGS[1] : newest_results_dir(joinpath(project_root, "..", "data_collection", "hpc", "results"))
    output_root = length(ARGS) >= 2 ? ARGS[2] : joinpath(project_root, "plots")

    println("Input root: $(input_root)")
    println("Output root: $(output_root)")

    run_id = Dates.format(now(), "yyyymmdd_HHMMSS")
    index_file = process_all_cases(input_root, output_root; run_id=run_id)

    println("Done. Plot index file:")
    println(index_file)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
