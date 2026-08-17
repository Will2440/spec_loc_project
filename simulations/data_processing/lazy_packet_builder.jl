using JLD2
using Dates
using Printf

function getv(x::AbstractDict, key::String)
    if haskey(x, key)
        return x[key]
    end
    skey = Symbol(key)
    if haskey(x, skey)
        return x[skey]
    end
    error("Key not found: $(key)")
end

function getv(x::NamedTuple, key::String)
    return getproperty(x, Symbol(key))
end

function value_to_str(x)
    if x isa Real
        return @sprintf("%.10g", Float64(x))
    end
    return string(x)
end

function round_key(x::Real; digits::Int=10)
    return round(Float64(x); digits=digits)
end

function ensure_dir(path::String)
    isdir(path) || mkpath(path)
    return path
end

function find_case_files(input_root::String)
    files = String[]
    for (root, _, fs) in walkdir(input_root)
        for f in fs
            endswith(f, ".jld2") || continue
            push!(files, joinpath(root, f))
        end
    end
    sort!(files)
    return files
end

function group_key(meta, axes, run_group::String)
    return (
        run_group=run_group,
        A=round_key(getv(meta, "A")),
        B=round_key(getv(meta, "B")),
        m=round_key(getv(meta, "m")),
        B_y=round_key(getv(meta, "B_y")),
        perturbation_type=String(getv(meta, "perturbation_type")),
        disorder_type=String(getv(meta, "disorder_type")),
        winding=Int(getv(meta, "winding_number")),
        Lx_ribbon=Int(getv(meta, "Lx_ribbon")),
        Lx_obc=Int(getv(meta, "Lx_obc")),
        Ly_obc=Int(getv(meta, "Ly_obc")),
        specloc_x=Int(round(Float64(getv(meta, "specloc_x")))),
        specloc_y=Int(round(Float64(getv(meta, "specloc_y")))),
    )
end

function add_record!(records, next_id, render_kind, group_id, case_file, plot_type, meta; gamma=NaN, W=NaN, kappa=NaN, E=NaN)
    winding_val = try
        Int(getv(meta, "winding_number"))
    catch
        Int(getv(meta, "winding"))
    end

    push!(records, (
        record_id=next_id,
        render_kind=render_kind,
        group_id=group_id,
        case_file=case_file,
        plot_type=plot_type,
        A=Float64(getv(meta, "A")),
        B=Float64(getv(meta, "B")),
        m=Float64(getv(meta, "m")),
        B_y=Float64(getv(meta, "B_y")),
        perturbation_type=String(getv(meta, "perturbation_type")),
        disorder_type=String(getv(meta, "disorder_type")),
        winding=winding_val,
        Lx_ribbon=Int(getv(meta, "Lx_ribbon")),
        Lx_obc=Int(getv(meta, "Lx_obc")),
        Ly_obc=Int(getv(meta, "Ly_obc")),
        specloc_x=Float64(getv(meta, "specloc_x")),
        specloc_y=Float64(getv(meta, "specloc_y")),
        gamma=Float64(gamma),
        W=Float64(W),
        kappa=Float64(kappa),
        E=Float64(E),
    ))
    return next_id + 1
end

function write_records_tsv(records, out_file)
    ensure_dir(dirname(out_file))
    header = [
        "record_id", "render_kind", "group_id", "case_file", "plot_type",
        "A", "B", "m", "B_y", "perturbation_type", "disorder_type", "winding",
        "Lx_ribbon", "Lx_obc", "Ly_obc", "specloc_x", "specloc_y", "gamma", "W", "kappa", "E"
    ]

    open(out_file, "w") do io
        println(io, join(header, '\t'))
        for r in records
            vals = [
                string(r.record_id), r.render_kind, r.group_id, r.case_file, r.plot_type,
                value_to_str(r.A), value_to_str(r.B), value_to_str(r.m), value_to_str(r.B_y),
                r.perturbation_type, r.disorder_type, string(r.winding),
                string(r.Lx_ribbon), string(r.Lx_obc), string(r.Ly_obc),
                value_to_str(r.specloc_x), value_to_str(r.specloc_y),
                value_to_str(r.gamma), value_to_str(r.W), value_to_str(r.kappa), value_to_str(r.E),
            ]
            println(io, join(vals, '\t'))
        end
    end
end

function write_groups_tsv(groups, out_file)
    ensure_dir(dirname(out_file))
    header = [
        "group_id", "run_group", "case_files", "gammas", "Ws", "kappas", "Es",
        "A", "B", "m", "B_y", "perturbation_type", "disorder_type", "winding",
        "Lx_ribbon", "Lx_obc", "Ly_obc", "specloc_x", "specloc_y"
    ]

    open(out_file, "w") do io
        println(io, join(header, '\t'))
        for g in groups
            vals = [
                g[:group_id],
                g[:run_group],
                join(g[:case_files], "||"),
                join(value_to_str.(g[:gammas]), ","),
                join(value_to_str.(g[:Ws]), ","),
                join(value_to_str.(g[:kappas]), ","),
                join(value_to_str.(g[:Es]), ","),
                value_to_str(g[:A]), value_to_str(g[:B]), value_to_str(g[:m]), value_to_str(g[:B_y]),
                g[:perturbation_type], g[:disorder_type], string(g[:winding]),
                string(g[:Lx_ribbon]), string(g[:Lx_obc]), string(g[:Ly_obc]),
                string(g[:specloc_x]), string(g[:specloc_y]),
            ]
            println(io, join(vals, '\t'))
        end
    end
end

function build_lazy_packet(input_root::String, output_root::String; run_id::String=Dates.format(now(), "yyyymmdd_HHMMSS"), verbose::Bool=true)
    case_files = find_case_files(input_root)
    isempty(case_files) && error("No .jld2 files found under $(input_root)")

    if verbose
        println("Found $(length(case_files)) files. Building lazy packet...")
    end

    packet_dir = ensure_dir(joinpath(output_root, run_id, "lazy_packet"))

    records = NamedTuple[]
    groups_by_key = Dict{Any, Dict{Symbol, Any}}()
    next_id = 1
    global_dcdE_min = Inf
    global_dcdE_max = -Inf
    split_by_run_group = lowercase(get(ENV, "SPECLOC_LAZY_GROUP_BY_RUN", "false")) in ("1", "true", "yes", "on")

    for (i, cf) in enumerate(case_files)
        verbose && (i % 100 == 0) && println("Scanned $(i)/$(length(case_files)) files")
        data = JLD2.load(cf)
        result = data["result"]

        meta = result["metadata"]
        axes = result["axes"]

        gammas = Float64.(collect(getv(axes, "gammas")))
        Ws = Float64.(collect(getv(axes, "Ws")))
        kappas = Float64.(collect(getv(axes, "kappas")))
        Es = Float64.(collect(getv(axes, "Es")))

        rel = relpath(cf, input_root)
        parts = splitpath(rel)
        run_group = split_by_run_group ? (isempty(parts) ? "root" : parts[1]) : "all"

        gkey = group_key(meta, axes, run_group)
        if !haskey(groups_by_key, gkey)
            gid = "group_" * string(length(groups_by_key) + 1)
            groups_by_key[gkey] = Dict{Symbol, Any}(
                :group_id => gid,
                :run_group => run_group,
                :case_files => String[],
                :gammas => Float64[],
                :Ws => Float64[],
                :kappas => Float64[],
                :Es => Es,
                :A => Float64(getv(meta, "A")),
                :B => Float64(getv(meta, "B")),
                :m => Float64(getv(meta, "m")),
                :B_y => Float64(getv(meta, "B_y")),
                :perturbation_type => String(getv(meta, "perturbation_type")),
                :disorder_type => String(getv(meta, "disorder_type")),
                :winding => Int(getv(meta, "winding_number")),
                :Lx_ribbon => Int(getv(meta, "Lx_ribbon")),
                :Lx_obc => Int(getv(meta, "Lx_obc")),
                :Ly_obc => Int(getv(meta, "Ly_obc")),
                :specloc_x => Int(round(Float64(getv(meta, "specloc_x")))),
                :specloc_y => Int(round(Float64(getv(meta, "specloc_y")))),
            )
        end

        grp = groups_by_key[gkey]
        push!(grp[:case_files], cf)
        append!(grp[:gammas], gammas)
        append!(grp[:Ws], Ws)
        append!(grp[:kappas], kappas)

        for gamma in gammas
            ribbon = result["ribbon_by_gamma"][gamma]
            chd = getv(ribbon, "chern_contribution_density_on_ribbon")
            local_min = minimum(chd)
            local_max = maximum(chd)
            if local_min < global_dcdE_min
                global_dcdE_min = local_min
            end
            if local_max > global_dcdE_max
                global_dcdE_max = local_max
            end

            next_id = add_record!(records, next_id, "case", "", cf, "ribbon_ipr_fixed", meta; gamma=gamma)
            next_id = add_record!(records, next_id, "case", "", cf, "ribbon_dcdE_fixed", meta; gamma=gamma)
            next_id = add_record!(records, next_id, "case", "", cf, "ribbon_chern_acc_fixed", meta; gamma=gamma)
            next_id = add_record!(records, next_id, "case", "", cf, "ribbon_ipr_auto", meta; gamma=gamma)
            next_id = add_record!(records, next_id, "case", "", cf, "ribbon_dcdE_auto", meta; gamma=gamma)
            next_id = add_record!(records, next_id, "case", "", cf, "ribbon_chern_acc_auto", meta; gamma=gamma)
            next_id = add_record!(records, next_id, "case", "", cf, "band3d", meta; gamma=gamma)
        end

        obc = result["obc"]
        has_dos_by_gamma = haskey(obc, "dos_by_gamma")
        has_ldos_target_by_gamma = haskey(obc, "ldos_target_by_gamma")
        has_ldos_lowest_by_gamma = haskey(obc, "ldos_lowest_by_gamma")

        if has_dos_by_gamma
            for gamma in gammas
                next_id = add_record!(records, next_id, "case", "", cf, "dos", meta; gamma=gamma)
            end
        else
            next_id = add_record!(records, next_id, "case", "", cf, "dos", meta)
        end

        if has_ldos_target_by_gamma
            for gamma in gammas
                next_id = add_record!(records, next_id, "case", "", cf, "ldos_target", meta; gamma=gamma)
            end
        else
            next_id = add_record!(records, next_id, "case", "", cf, "ldos_target", meta)
        end

        if has_ldos_lowest_by_gamma
            for gamma in gammas
                next_id = add_record!(records, next_id, "case", "", cf, "ldos_lowest", meta; gamma=gamma)
            end
        else
            next_id = add_record!(records, next_id, "case", "", cf, "ldos_lowest", meta)
        end
    end

    specloc_plot_types = [
        "specloc_signature_vs_E",
        "specloc_loggap_vs_E",
        "specloc_gap_gamma_E",
        "specloc_sig_gamma_E",
        "specloc_spectrum_vs_gamma",
    ]

    groups = collect(values(groups_by_key))
    for grp in groups
        grp[:gammas] = sort(unique(Float64.(grp[:gammas])))
        grp[:Ws] = sort(unique(Float64.(grp[:Ws])))
        grp[:kappas] = sort(unique(Float64.(grp[:kappas])))

        g0 = grp[:gammas][cld(length(grp[:gammas]), 2)]
        W0 = grp[:Ws][cld(length(grp[:Ws]), 2)]
        k0 = grp[:kappas][cld(length(grp[:kappas]), 2)]
        E0 = grp[:Es][cld(length(grp[:Es]), 2)]

        for pt in specloc_plot_types
            next_id = add_record!(records, next_id, "specloc_group", grp[:group_id], "", pt, grp;
                gamma=g0, W=W0, kappa=k0, E=E0)
        end
    end

    records_file = joinpath(packet_dir, "lazy_packet_records.tsv")
    groups_file = joinpath(packet_dir, "lazy_packet_groups.tsv")
    stats_file = joinpath(packet_dir, "lazy_packet_stats.tsv")
    write_records_tsv(records, records_file)
    write_groups_tsv(groups, groups_file)
    open(stats_file, "w") do io
        println(io, "metric\tvalue")
        if isfinite(global_dcdE_min)
            println(io, "dcdE_min\t" * value_to_str(global_dcdE_min))
        end
        if isfinite(global_dcdE_max)
            println(io, "dcdE_max\t" * value_to_str(global_dcdE_max))
        end
    end

    if verbose
        println("Lazy packet created:")
        println("  records: $(records_file)")
        println("  groups:  $(groups_file)")
        println("  stats:   $(stats_file)")
        println("  records count: $(length(records))")
        println("  groups count:  $(length(groups))")
    end

    return (records_file=records_file, groups_file=groups_file)
end

function main()
    if length(ARGS) < 2
        println("Usage: julia lazy_packet_builder.jl <input_root> <output_root> [run_id]")
        exit(1)
    end

    input_root = ARGS[1]
    output_root = ARGS[2]
    run_id = length(ARGS) >= 3 ? ARGS[3] : Dates.format(now(), "yyyymmdd_HHMMSS")

    build_lazy_packet(input_root, output_root; run_id=run_id, verbose=true)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
