using JLD2
using Dates

# Force headless rendering for CLI calls from viewer.
# Use a file-producing workstation to avoid zero-byte PNG outputs.
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

function parse_float(x::AbstractString)
    s = strip(x)
    if isempty(s)
        return NaN
    end
    if lowercase(s) == "nan"
        return NaN
    end
    return parse(Float64, s)
end

function read_tsv(path::String)
    lines = readlines(path)
    isempty(lines) && return Dict{String, String}[]
    header = split(lines[1], '\t')
    rows = Vector{Dict{String, String}}()
    for line in lines[2:end]
        isempty(strip(line)) && continue
        vals = split(line, '\t')
        d = Dict{String, String}()
        for (k, v) in zip(header, vals)
            d[k] = v
        end
        push!(rows, d)
    end
    return rows
end

function group_merge_key(g::Dict{String, String})
    keys = [
        "A", "B", "m", "B_y", "perturbation_type", "disorder_type", "winding",
        "Lx_ribbon", "Lx_obc", "Ly_obc", "specloc_x", "specloc_y"
    ]
    return join([get(g, k, "") for k in keys], "|")
end

function merge_groups_by_physics(groups::Vector{Dict{String, String}})
    merged_by_key = Dict{String, Dict{String, String}}()
    merged_for_gid = Dict{String, Dict{String, String}}()

    for g in groups
        key = group_merge_key(g)
        if !haskey(merged_by_key, key)
            merged = copy(g)
            merged["run_group"] = "merged"
            merged_by_key[key] = merged
        else
            merged = merged_by_key[key]
            existing = Set([x for x in split(get(merged, "case_files", ""), "||") if !isempty(strip(x))])
            incoming = [x for x in split(get(g, "case_files", ""), "||") if !isempty(strip(x))]
            for cf in incoming
                push!(existing, cf)
            end
            merged["case_files"] = join(sort(collect(existing)), "||")
        end
        merged_for_gid[g["group_id"]] = merged_by_key[key]
    end

    return merged_for_gid
end

function load_packet_dcdE_clims(records_tsv::String)
    stats_file = joinpath(dirname(records_tsv), "lazy_packet_stats.tsv")
    isfile(stats_file) || return nothing
    rows = read_tsv(stats_file)
    dmin = NaN
    dmax = NaN
    for r in rows
        metric = get(r, "metric", "")
        if metric == "dcdE_min"
            dmin = parse_float(get(r, "value", "NaN"))
        elseif metric == "dcdE_max"
            dmax = parse_float(get(r, "value", "NaN"))
        end
    end
    if isfinite(dmin) && isfinite(dmax)
        return (dmin, dmax)
    end
    return nothing
end

function write_manifest(index_rows::Vector{NamedTuple}, out_file::String)
    mkpath(dirname(out_file))
    open(out_file, "w") do io
        println(io, "plot_type\tgamma\tW\tkappa\tE\tplot_path")
        for r in index_rows
            println(io, string(r.plot_type, '\t', r.gamma, '\t', r.W, '\t', r.kappa, '\t', r.E, '\t', r.plot_path))
        end
    end
end

function read_manifest(path::String)
    rows = read_tsv(path)
    out = NamedTuple[]
    for r in rows
        push!(out, (
            plot_type=r["plot_type"],
            gamma=parse_float(r["gamma"]),
            W=parse_float(r["W"]),
            kappa=parse_float(r["kappa"]),
            E=parse_float(r["E"]),
            plot_path=r["plot_path"],
        ))
    end
    return out
end

function manifest_valid(path::String)
    isfile(path) || return false
    rows = read_manifest(path)
    isempty(rows) && return false
    for r in rows
        p = r.plot_path
        isfile(p) || return false
        filesize(p) > 0 || return false
    end
    return true
end

function row_distance(man_row, req_gamma, req_W, req_kappa, req_E)
    dist = 0.0
    used = 0
    vals = ((man_row.gamma, req_gamma), (man_row.W, req_W), (man_row.kappa, req_kappa), (man_row.E, req_E))
    for (a, b) in vals
        if isfinite(a) && isfinite(b)
            dist += (a - b)^2
            used += 1
        end
    end
    return used == 0 ? 0.0 : dist
end

function resolve_path_from_manifest(man_rows, plot_type::String, req_gamma, req_W, req_kappa, req_E)
    candidates = [r for r in man_rows if r.plot_type == plot_type]
    isempty(candidates) && return ""
    best = candidates[argmin([row_distance(r, req_gamma, req_W, req_kappa, req_E) for r in candidates])]
    return best.plot_path
end

function render_case_bundle(record::Dict{String, String}, cache_root::String; fixed_dcdE_clims=nothing, force_rebuild::Bool=false)
    case_file = record["case_file"]
    case_key = string(hash(abspath(case_file)))
    bundle_dir = joinpath(cache_root, "case_" * case_key)
    manifest_file = joinpath(bundle_dir, "manifest.tsv")

    if force_rebuild && isdir(bundle_dir)
        rm(bundle_dir; recursive=true, force=true)
    end

    if manifest_valid(manifest_file)
        return read_manifest(manifest_file)
    elseif isdir(bundle_dir)
        rm(bundle_dir; recursive=true, force=true)
    end

    data = JLD2.load(case_file)
    result = data["result"]

    meta = result["metadata"]
    axes = result["axes"]
    obc = result["obc"]

    index_rows = NamedTuple[]

    ribbon_dir = SpecLocPlotting.ensure_dir(joinpath(bundle_dir, "ribbon"))
    band3d_dir = SpecLocPlotting.ensure_dir(joinpath(bundle_dir, "band3d"))

    for gamma in axes["gammas"]
        ribbon = result["ribbon_by_gamma"][gamma]
        band3d = result["band3d_by_gamma"][gamma]
        SpecLocPlotting.plot_ribbon_triplet!(index_rows, meta, "lazy", case_file, ribbon_dir, gamma, ribbon; fixed_dcdE_clims=fixed_dcdE_clims)
        SpecLocPlotting.plot_band3d!(index_rows, meta, "lazy", case_file, band3d_dir, gamma, band3d)
    end

    SpecLocPlotting.plot_obc!(index_rows, meta, "lazy", case_file, SpecLocPlotting.ensure_dir(joinpath(bundle_dir, "obc")), axes, obc)

    write_manifest(index_rows, manifest_file)
    return index_rows
end

function parse_csv_floats(s::String)
    isempty(strip(s)) && return Float64[]
    return [parse(Float64, x) for x in split(s, ",") if !isempty(strip(x))]
end

function render_specloc_group_bundle(record::Dict{String, String}, groups_by_id::Dict{String, Dict{String, String}}, cache_root::String; force_rebuild::Bool=false)
    rid = record["record_id"]
    gid = record["group_id"]
    haskey(groups_by_id, gid) || error("Group id not found: " * gid)
    grp = groups_by_id[gid]
    case_files = [String(x) for x in split(grp["case_files"], "||") if !isempty(strip(x))]
    bundle_key = string(hash(join(sort(case_files), "||")))
    bundle_dir = joinpath(cache_root, "specloc_" * bundle_key)
    manifest_file = joinpath(bundle_dir, "manifest.tsv")

    if force_rebuild && isdir(bundle_dir)
        rm(bundle_dir; recursive=true, force=true)
    end

    if manifest_valid(manifest_file)
        return read_manifest(manifest_file)
    elseif isdir(bundle_dir)
        rm(bundle_dir; recursive=true, force=true)
    end

    run_group = grp["run_group"]

    specloc_groups = Dict{Any, Any}()
    for cf in case_files
        payload = SpecLocPlotting.load_case_specloc_payload(cf)
        SpecLocPlotting.add_specloc_case!(specloc_groups, cf, run_group, payload.meta, payload.axes, payload.specloc)
    end

    index_rows = NamedTuple[]
    for (_, g) in specloc_groups
        axes_agg, specloc_agg = SpecLocPlotting.build_aggregated_specloc_payload(g)
        out_dir = SpecLocPlotting.ensure_dir(bundle_dir)
        case_label = "lazy_group::" * string(length(g[:case_files])) * " files"
        SpecLocPlotting.plot_specloc_cuts!(index_rows, g[:meta], "lazy", case_label, out_dir, axes_agg, specloc_agg)
    end

    write_manifest(index_rows, manifest_file)
    return index_rows
end

function main()
    if length(ARGS) != 4
        println("Usage: julia lazy_render.jl <records_tsv> <groups_tsv> <cache_root> <requests_tsv>")
        exit(1)
    end

    records_tsv = ARGS[1]
    groups_tsv = ARGS[2]
    cache_root = ARGS[3]
    requests_tsv = ARGS[4]

    mkpath(cache_root)

    records = read_tsv(records_tsv)
    groups = read_tsv(groups_tsv)
    requests = read_tsv(requests_tsv)

    rec_by_id = Dict(r["record_id"] => r for r in records)
    groups_by_id = merge_groups_by_physics(groups)
    fixed_dcdE_clims = load_packet_dcdE_clims(records_tsv)
    force_rebuild = lowercase(get(ENV, "SPECLOC_LAZY_FORCE_REBUILD", "false")) in ("1", "true", "yes", "on")

    suppress_gr_warnings = lowercase(get(ENV, "SPECLOC_SUPPRESS_GR_WARNINGS", "true")) in ("1", "true", "yes", "on")

    function warm_records!()
        rendered_local = Dict{String, Vector{NamedTuple}}()
        req_ids = unique([r["record_id"] for r in requests])
        for rid in req_ids
            haskey(rec_by_id, rid) || continue
            rec = rec_by_id[rid]
            if rec["render_kind"] == "case"
                rendered_local[rid] = render_case_bundle(rec, cache_root; fixed_dcdE_clims=fixed_dcdE_clims, force_rebuild=force_rebuild)
            elseif rec["render_kind"] == "specloc_group"
                rendered_local[rid] = render_specloc_group_bundle(rec, groups_by_id, cache_root; force_rebuild=force_rebuild)
            end
        end
        return rendered_local
    end

    rendered = if suppress_gr_warnings
        redirect_stderr(devnull) do
            warm_records!()
        end
    else
        warm_records!()
    end

    # Output resolved paths.
    for req in requests
        rid = req["record_id"]
        ptype = req["plot_type"]
        rec = rec_by_id[rid]
        man_rows = get(rendered, rid, NamedTuple[])

        req_gamma = parse_float(get(rec, "gamma", "NaN"))
        req_W = parse_float(get(rec, "W", "NaN"))
        req_kappa = parse_float(get(rec, "kappa", "NaN"))
        req_E = parse_float(get(rec, "E", "NaN"))

        path = resolve_path_from_manifest(man_rows, ptype, req_gamma, req_W, req_kappa, req_E)
        println(string(rid, '\t', ptype, '\t', path))
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
