module SpecLocPlotting

using JLD2
using Plots
using Printf
using Dates
using ProgressMeter
using Distributed

export process_all_cases

getv(x::AbstractDict, key::String) = x[key]
getv(x::AbstractDict, key::Symbol) = x[key]
getv(x::NamedTuple, key::String) = getproperty(x, Symbol(key))
getv(x::NamedTuple, key::Symbol) = getproperty(x, key)
haskv(x::AbstractDict, key::String) = haskey(x, key)
haskv(x::AbstractDict, key::Symbol) = haskey(x, key)
haskv(x::NamedTuple, key::String) = hasproperty(x, Symbol(key))
haskv(x::NamedTuple, key::Symbol) = hasproperty(x, key)

function nearest_index(vals::AbstractVector{<:Real}, x::Real)
    return argmin(abs.(vals .- x))
end

function ensure_dir(path::String)
    isdir(path) || mkpath(path)
    return path
end

function value_to_str(x)
    if x isa Real
        return @sprintf("%.6g", Float64(x))
    end
    return string(x)
end

function approx_vec_equal(a::AbstractVector{<:Real}, b::AbstractVector{<:Real}; atol::Float64=1e-10)
    length(a) == length(b) || return false
    @inbounds for i in eachindex(a, b)
        if !isapprox(Float64(a[i]), Float64(b[i]); atol=atol, rtol=0.0)
            return false
        end
    end
    return true
end

function infer_gap_closure_from_ipr(ky_vals, energies, iprs; ipr_threshold::Float64=0.05, gap_tol::Float64=0.05)
    N_ky = length(ky_vals)
    nbands = size(energies, 1)
    size(iprs) == size(energies) || error("iprs and energies must have the same shape.")

    edge_gap_margin = fill(-Inf, N_ky)

    for i in 1:N_ky
        edge_idx = findall(iprs[:, i] .>= ipr_threshold)
        bulk_idx = findall(iprs[:, i] .< ipr_threshold)

        isempty(edge_idx) && continue
        isempty(bulk_idx) && continue

        bulk_energies = @view energies[bulk_idx, i]
        local_best = -Inf

        for idx in edge_idx
            E = energies[idx, i]
            lower_candidates = bulk_energies[bulk_energies .< E]
            upper_candidates = bulk_energies[bulk_energies .> E]
            (isempty(lower_candidates) || isempty(upper_candidates)) && continue

            lower_bulk = maximum(lower_candidates)
            upper_bulk = minimum(upper_candidates)
            margin = min(E - lower_bulk, upper_bulk - E)
            local_best = max(local_best, margin)
        end

        edge_gap_margin[i] = local_best
    end

    edge_in_gap_mask = edge_gap_margin .> gap_tol

    # Remove short numerical glitches that can split one physical interval into many.
    function smooth_mask(mask::BitVector; min_true_run::Int=3, max_false_gap::Int=2)
        out = copy(mask)
        n = length(out)

        i = 1
        while i <= n
            if out[i]
                j = i
                while j < n && out[j + 1]
                    j += 1
                end
                if (j - i + 1) < min_true_run
                    out[i:j] .= false
                end
                i = j + 1
            else
                i += 1
            end
        end

        i = 1
        while i <= n
            if !out[i]
                j = i
                while j < n && !out[j + 1]
                    j += 1
                end
                left_true = i > 1 && out[i - 1]
                right_true = j < n && out[j + 1]
                if left_true && right_true && (j - i + 1) <= max_false_gap
                    out[i:j] .= true
                end
                i = j + 1
            else
                i += 1
            end
        end

        return out
    end

    edge_in_gap_mask = smooth_mask(BitVector(edge_in_gap_mask))
    if !any(edge_in_gap_mask)
        return (
            ky_closure_bounds=Float64[],
            extent_in_pi=0.0,
            edge_in_gap_mask=edge_in_gap_mask,
            edge_gap_margin=edge_gap_margin,
        )
    end

    runs = Vector{Vector{Int}}()
    current = Int[]
    for i in 1:N_ky
        if edge_in_gap_mask[i]
            push!(current, i)
        elseif !isempty(current)
            push!(runs, current)
            current = Int[]
        end
    end
    !isempty(current) && push!(runs, current)

    if length(runs) > 1 && edge_in_gap_mask[1] && edge_in_gap_mask[end]
        merged = vcat(runs[end], runs[1])
        middle = length(runs) > 2 ? runs[2:end-1] : Vector{Vector{Int}}()
        runs = vcat([merged], middle)
    end

    run_lengths = [length(r) for r in runs]
    dom = runs[argmax(run_lengths)]
    dom_start = first(dom)
    dom_end = last(dom)

    period = 2 * pi
    ky_min, ky_max = extrema(ky_vals)

    function wrap_to_window(k::Float64)
        while k < ky_min
            k += period
        end
        while k > ky_max
            k -= period
        end
        return k
    end

    function interpolate_threshold(k1::Float64, g1::Float64, k2::Float64, g2::Float64)
        k2_adj = k2
        if k2_adj < k1
            k2_adj += period
        end
        if g1 == g2
            return wrap_to_window(0.5 * (k1 + k2_adj))
        end
        k_interp = k1 + (gap_tol - g1) * (k2_adj - k1) / (g2 - g1)
        return wrap_to_window(k_interp)
    end

    prev_idx = dom_start == 1 ? N_ky : (dom_start - 1)
    next_idx = dom_end == N_ky ? 1 : (dom_end + 1)

    left_bound = interpolate_threshold(
        Float64(ky_vals[prev_idx]),
        edge_gap_margin[prev_idx],
        Float64(ky_vals[dom_start]),
        edge_gap_margin[dom_start],
    )

    right_bound = interpolate_threshold(
        Float64(ky_vals[dom_end]),
        edge_gap_margin[dom_end],
        Float64(ky_vals[next_idx]),
        edge_gap_margin[next_idx],
    )

    return (
        ky_closure_bounds=[left_bound, right_bound],
        extent_in_pi=mod(right_bound - left_bound, period) / pi,
        edge_in_gap_mask=edge_in_gap_mask,
        edge_gap_margin=edge_gap_margin,
    )
end

function round_key(x::Real; digits::Int=10)
    return round(Float64(x); digits=digits)
end

function make_specloc_group_key(meta, axes, run_group::String)
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

function add_specloc_case!(groups::Dict{Any, Any}, case_file::String, run_group::String, meta, axes, specloc)
    key = make_specloc_group_key(meta, axes, run_group)
    gammas = Float64.(collect(getv(axes, "gammas")))
    Ws = Float64.(collect(getv(axes, "Ws")))
    kappas = Float64.(collect(getv(axes, "kappas")))
    Es = Float64.(collect(getv(axes, "Es")))

    gap = getv(specloc, "gap")
    sig = getv(specloc, "signature")

    if !haskey(groups, key)
        groups[key] = Dict{Symbol, Any}(
            :case_files => String[],
            :gammas => sort(unique(gammas)),
            :Ws => sort(unique(Ws)),
            :kappas => sort(unique(kappas)),
            :Es => Es,
            :gap_map => Dict{NTuple{4, Float64}, Float64}(),
            :sig_map => Dict{NTuple{4, Float64}, Float64}(),
            :spec_g_map => Dict{Float64, Vector{Float64}}(),
            :spec_W_map => Dict{Float64, Vector{Float64}}(),
            :spec_k_map => Dict{Float64, Vector{Float64}}(),
            :meta => meta,
        )
    end

    grp = groups[key]
    push!(grp[:case_files], case_file)

    grp[:gammas] = sort(unique(vcat(grp[:gammas], gammas)))
    grp[:Ws] = sort(unique(vcat(grp[:Ws], Ws)))
    grp[:kappas] = sort(unique(vcat(grp[:kappas], kappas)))
    grp[:Es] = sort(unique(vcat(Float64.(grp[:Es]), Es)))

    for (gi, g) in enumerate(gammas), (wi, w) in enumerate(Ws), (ki, k) in enumerate(kappas), (ei, e) in enumerate(Es)
        tuple_key = (round_key(g), round_key(w), round_key(k), round_key(e))
        grp[:gap_map][tuple_key] = Float64(gap[gi, wi, ki, ei])
        grp[:sig_map][tuple_key] = Float64(sig[gi, wi, ki, ei])
    end

    svg = getv(specloc, "spectrum_vs_gamma")
    svW = getv(specloc, "spectrum_vs_W")
    svk = getv(specloc, "spectrum_vs_kappa")

    for (i, vals) in enumerate(svg)
        grp[:spec_g_map][round_key(gammas[i])] = Float64.(collect(vals))
    end
    for (i, vals) in enumerate(svW)
        grp[:spec_W_map][round_key(Ws[i])] = Float64.(collect(vals))
    end
    for (i, vals) in enumerate(svk)
        grp[:spec_k_map][round_key(kappas[i])] = Float64.(collect(vals))
    end
end

function group_slug(key)
    return @sprintf(
        "run-%s_A%s_B%s_m%s_By%s_pt%s_dt%s_w%d_Lxr%d_Lxo%d_Lyo%d_sx%d_sy%d",
        key.run_group,
        value_to_str(key.A),
        value_to_str(key.B),
        value_to_str(key.m),
        value_to_str(key.B_y),
        key.perturbation_type,
        key.disorder_type,
        key.winding,
        key.Lx_ribbon,
        key.Lx_obc,
        key.Ly_obc,
        key.specloc_x,
        key.specloc_y,
    )
end

function build_aggregated_specloc_payload(group)
    gammas = Float64.(group[:gammas])
    Ws = Float64.(group[:Ws])
    kappas = Float64.(group[:kappas])
    Es = Float64.(group[:Es])

    gap_arr = fill(NaN, length(gammas), length(Ws), length(kappas), length(Es))
    sig_arr = fill(NaN, length(gammas), length(Ws), length(kappas), length(Es))

    gmap = Dict(round_key(g) => i for (i, g) in enumerate(gammas))
    Wmap = Dict(round_key(w) => i for (i, w) in enumerate(Ws))
    kmap = Dict(round_key(k) => i for (i, k) in enumerate(kappas))
    Emap = Dict(round_key(e) => i for (i, e) in enumerate(Es))

    for (k, v) in group[:gap_map]
        gi = gmap[k[1]]
        wi = Wmap[k[2]]
        ki = kmap[k[3]]
        ei = Emap[k[4]]
        gap_arr[gi, wi, ki, ei] = v
    end
    for (k, v) in group[:sig_map]
        gi = gmap[k[1]]
        wi = Wmap[k[2]]
        ki = kmap[k[3]]
        ei = Emap[k[4]]
        sig_arr[gi, wi, ki, ei] = v
    end

    svg = [get(group[:spec_g_map], round_key(g), Float64[]) for g in gammas]
    svW = [get(group[:spec_W_map], round_key(w), Float64[]) for w in Ws]
    svk = [get(group[:spec_k_map], round_key(k), Float64[]) for k in kappas]

    axes = Dict("gammas" => gammas, "Ws" => Ws, "kappas" => kappas, "Es" => Es)
    specloc = Dict(
        "gap" => gap_arr,
        "signature" => sig_arr,
        "spectrum_vs_gamma" => svg,
        "spectrum_vs_W" => svW,
        "spectrum_vs_kappa" => svk,
    )
    return axes, specloc
end

function save_plot_index_tsv(rows::Vector{NamedTuple}, out_file::String)
    ensure_dir(dirname(out_file))
    header = [
        "run_id", "case_file", "plot_type", "plot_path",
        "A", "B", "m", "B_y", "perturbation_type", "disorder_type", "winding",
        "Lx_ribbon", "Lx_obc", "Ly_obc", "specloc_x", "specloc_y", "gamma", "W", "kappa", "E"
    ]

    open(out_file, "w") do io
        println(io, join(header, '\t'))
        for r in rows
            vals = [
                r.run_id, r.case_file, r.plot_type, r.plot_path,
                value_to_str(r.A), value_to_str(r.B), value_to_str(r.m), value_to_str(r.B_y),
                string(r.perturbation_type), string(r.disorder_type), value_to_str(r.winding),
                value_to_str(r.Lx_ribbon), value_to_str(r.Lx_obc), value_to_str(r.Ly_obc),
                value_to_str(r.specloc_x), value_to_str(r.specloc_y),
                value_to_str(r.gamma), value_to_str(r.W), value_to_str(r.kappa), value_to_str(r.E),
            ]
            println(io, join(vals, '\t'))
        end
    end
end

function add_index_row!(rows, meta, run_id, case_file, plot_type, plot_path; gamma=NaN, W=NaN, kappa=NaN, E=NaN)
    x0 = getv(meta, "specloc_x")
    y0 = getv(meta, "specloc_y")
    push!(rows, (
        run_id=run_id,
        case_file=case_file,
        plot_type=plot_type,
        plot_path=plot_path,
        A=Float64(getv(meta, "A")),
        B=Float64(getv(meta, "B")),
        m=Float64(getv(meta, "m")),
        B_y=Float64(getv(meta, "B_y")),
        perturbation_type=String(getv(meta, "perturbation_type")),
        disorder_type=String(getv(meta, "disorder_type")),
        winding=Int(getv(meta, "winding_number")),
        Lx_ribbon=Int(getv(meta, "Lx_ribbon")),
        Lx_obc=Int(getv(meta, "Lx_obc")),
        Ly_obc=Int(getv(meta, "Ly_obc")),
        specloc_x=Float64(x0),
        specloc_y=Float64(y0),
        gamma=Float64(gamma),
        W=Float64(W),
        kappa=Float64(kappa),
        E=Float64(E),
    ))
end

function plot_ribbon_triplet!(index_rows, meta, run_id, case_file, out_dir, gamma, ribbon; fixed_dcdE_clims=nothing)
    ensure_dir(out_dir)

    ky = getv(ribbon, "ky_vals")
    energies = getv(ribbon, "energies")
    ipr = getv(ribbon, "iprs")
    chd = getv(ribbon, "chern_contribution_density_on_ribbon")
    accum_e = getv(ribbon, "accumulation_energies")
    accum_c = getv(ribbon, "cumulative_chern")
    max_ipr_val = maximum(ipr)

    gap_data = nothing
    if haskv(ribbon, "gap_closure_data")
        gap_data = getv(ribbon, "gap_closure_data")
    else
        gap_data = infer_gap_closure_from_ipr(ky, energies, ipr)
    end

    chlim = maximum(abs.(chd))

    # Function to add gap closure annotations
    function add_gap_annotations!(plt, ky, energies, gap_data, max_ipr_val)
        if !isnothing(gap_data)
            if haskv(gap_data, "edge_in_gap_mask")
                mask = getv(gap_data, "edge_in_gap_mask")
                if length(mask) == length(ky)
                    N = length(mask)
                    i = 1
                    while i <= N
                        if !mask[i]
                            j = i
                            while j < N && !mask[j + 1]
                                j += 1
                            end

                            k_left = i == 1 ? ky[1] : 0.5 * (ky[i - 1] + ky[i])
                            k_right = j == N ? ky[end] : 0.5 * (ky[j] + ky[j + 1])
                            vspan!(plt, [k_left, k_right]; color=:black, alpha=0.05, label=false)
                            i = j + 1
                        else
                            i += 1
                        end
                    end
                end
            end

            if haskv(gap_data, "ky_closure_bounds")
                closure_bounds = getv(gap_data, "ky_closure_bounds")
                if length(closure_bounds) == 2
                    for ky_closure in closure_bounds
                        vline!(plt, [ky_closure]; linestyle=:dot, color=:red, lw=1.5, alpha=0.7, label=false)
                    end
                end
            end

            if haskv(gap_data, "extent_in_pi")
                extent = getv(gap_data, "extent_in_pi")
                ky_range = extrema(ky)
                e_range = extrema(energies)
                ky_pos = ky_range[1] + 0.02 * (ky_range[2] - ky_range[1])
                e_pos = e_range[2] - 0.05 * (e_range[2] - e_range[1])
                annotation_text = @sprintf("Edge extent (gap): %.2fpi\\nMax IPR: %.3f", extent, max_ipr_val)
                annotate!(plt, ky_pos, e_pos, text(annotation_text, 8, :left, :top, :black))
            end
        end
    end

    # === FIXED LIMITS VERSION ===
    p1_fixed = plot(title="Ribbon IPR gamma=$(gamma)", xlabel="ky", ylabel="Energy", legend=false, colorbar=true, colorbar_title="IPR", clims=(0.0, 1.0))
    for b in 1:size(energies, 1)
        plot!(p1_fixed, ky, @view(energies[b, :]); line_z=@view(ipr[b, :]), lw=1.0, c=:viridis, label=false)
    end
    add_gap_annotations!(p1_fixed, ky, energies, gap_data, max_ipr_val)

    dcdE_clims = isnothing(fixed_dcdE_clims) ? (-chlim, chlim) : fixed_dcdE_clims
    p2_fixed = plot(title="Ribbon dC/dE gamma=$(gamma)", xlabel="ky", ylabel="Energy", legend=false, colorbar=true, colorbar_title="dC/dE", clims=dcdE_clims)
    for b in 1:size(energies, 1)
        plot!(p2_fixed, ky, @view(energies[b, :]); line_z=@view(chd[b, :]), lw=1.0, c=:RdBu, label=false)
    end

    p3_fixed = plot(accum_c, accum_e; xlabel="Accumulated Chern", ylabel="Energy", title="Chern accumulation gamma=$(gamma)", lw=2, color=:darkred, legend=false, xlims=(-1.0, 1.0))
    vline!(p3_fixed, [0.0]; linestyle=:dash, color=:black, lw=1, label=false)
    hline!(p3_fixed, [0.0]; linestyle=:dash, color=:black, lw=1, label=false)

    fp1_fixed = joinpath(out_dir, "ribbon_ipr_gamma$(gamma)_fixed.png")
    fp2_fixed = joinpath(out_dir, "ribbon_dcdE_gamma$(gamma)_fixed.png")
    fp3_fixed = joinpath(out_dir, "ribbon_chern_acc_gamma$(gamma)_fixed.png")

    savefig(p1_fixed, fp1_fixed)
    savefig(p2_fixed, fp2_fixed)
    savefig(p3_fixed, fp3_fixed)

    add_index_row!(index_rows, meta, run_id, case_file, "ribbon_ipr_fixed", fp1_fixed; gamma=gamma)
    add_index_row!(index_rows, meta, run_id, case_file, "ribbon_dcdE_fixed", fp2_fixed; gamma=gamma)
    add_index_row!(index_rows, meta, run_id, case_file, "ribbon_chern_acc_fixed", fp3_fixed; gamma=gamma)

    # === AUTO LIMITS VERSION ===
    p1_auto = plot(title="Ribbon IPR gamma=$(gamma) [auto]", xlabel="ky", ylabel="Energy", legend=false, colorbar=true, colorbar_title="IPR")
    for b in 1:size(energies, 1)
        plot!(p1_auto, ky, @view(energies[b, :]); line_z=@view(ipr[b, :]), lw=1.0, c=:viridis, label=false)
    end
    add_gap_annotations!(p1_auto, ky, energies, gap_data, max_ipr_val)

    p2_auto = plot(title="Ribbon dC/dE gamma=$(gamma) [auto]", xlabel="ky", ylabel="Energy", legend=false, colorbar=true, colorbar_title="dC/dE")
    for b in 1:size(energies, 1)
        plot!(p2_auto, ky, @view(energies[b, :]); line_z=@view(chd[b, :]), lw=1.0, c=:RdBu, label=false)
    end

    p3_auto = plot(accum_c, accum_e; xlabel="Accumulated Chern", ylabel="Energy", title="Chern accumulation gamma=$(gamma) [auto]", lw=2, color=:darkred, legend=false)
    vline!(p3_auto, [0.0]; linestyle=:dash, color=:black, lw=1, label=false)
    hline!(p3_auto, [0.0]; linestyle=:dash, color=:black, lw=1, label=false)

    fp1_auto = joinpath(out_dir, "ribbon_ipr_gamma$(gamma)_auto.png")
    fp2_auto = joinpath(out_dir, "ribbon_dcdE_gamma$(gamma)_auto.png")
    fp3_auto = joinpath(out_dir, "ribbon_chern_acc_gamma$(gamma)_auto.png")

    savefig(p1_auto, fp1_auto)
    savefig(p2_auto, fp2_auto)
    savefig(p3_auto, fp3_auto)

    add_index_row!(index_rows, meta, run_id, case_file, "ribbon_ipr_auto", fp1_auto; gamma=gamma)
    add_index_row!(index_rows, meta, run_id, case_file, "ribbon_dcdE_auto", fp2_auto; gamma=gamma)
    add_index_row!(index_rows, meta, run_id, case_file, "ribbon_chern_acc_auto", fp3_auto; gamma=gamma)
end

function plot_band3d!(index_rows, meta, run_id, case_file, out_dir, gamma, band3d)
    ensure_dir(out_dir)
    kx = getv(band3d, "kx_vals")
    ky = getv(band3d, "ky_vals")
    em = getv(band3d, "e_minus")
    ep = getv(band3d, "e_plus")

    plt = surface(kx, ky, em'; xlabel="kx", ylabel="ky", zlabel="E", color=:viridis, alpha=1.0, legend=false, title="3D Bandstructure gamma=$(gamma)")
    surface!(plt, kx, ky, ep'; color=:plasma, alpha=0.5, legend=false)

    fp = joinpath(out_dir, "band3d_gamma$(gamma).png")
    savefig(plt, fp)
    add_index_row!(index_rows, meta, run_id, case_file, "band3d", fp; gamma=gamma)
end

function plot_obc!(index_rows, meta, run_id, case_file, out_dir, axes, obc)
    ensure_dir(out_dir)

    E = getv(axes, "Es")
    Etarget = getv(obc, "ldos_target_E")
    meta_tag = @sprintf(
        "A=%s B=%s m=%s By=%s pt=%s dt=%s w=%d sx=%d sy=%d",
        value_to_str(getv(meta, "A")),
        value_to_str(getv(meta, "B")),
        value_to_str(getv(meta, "m")),
        value_to_str(getv(meta, "B_y")),
        String(getv(meta, "perturbation_type")),
        String(getv(meta, "disorder_type")),
        Int(getv(meta, "winding_number")),
        Int(round(Float64(getv(meta, "specloc_x")))),
        Int(round(Float64(getv(meta, "specloc_y")))),
    )

    has_dos_by_gamma = haskv(obc, "dos_by_gamma")
    has_ldos_target_by_gamma = haskv(obc, "ldos_target_by_gamma")
    has_ldos_lowest_by_gamma = haskv(obc, "ldos_lowest_by_gamma")

    dos_by_gamma = has_dos_by_gamma ? getv(obc, "dos_by_gamma") : nothing
    ldos_target_by_gamma = has_ldos_target_by_gamma ? getv(obc, "ldos_target_by_gamma") : nothing
    ldos_lowest_by_gamma = has_ldos_lowest_by_gamma ? getv(obc, "ldos_lowest_by_gamma") : nothing
    gammas = Float64.(collect(getv(axes, "gammas")))

    if has_dos_by_gamma
        for gamma in gammas
            gkey = gamma
            haskv(dos_by_gamma, gkey) || continue
            dos = getv(dos_by_gamma, gkey)
            p_dos = plot(E, dos; xlabel="Energy", ylabel="DOS", title="Total DOS gamma=$(gamma)\n$(meta_tag)", legend=false, color=:steelblue)
            gslug = replace(value_to_str(gamma), "." => "p", "-" => "m")
            fp_dos = joinpath(out_dir, "dos_vs_energy_gamma_$(gslug).png")
            savefig(p_dos, fp_dos)
            add_index_row!(index_rows, meta, run_id, case_file, "dos", fp_dos; gamma=gamma)
        end
    else
        dos = getv(obc, "dos")
        p_dos = plot(E, dos; xlabel="Energy", ylabel="DOS", title="Total DOS (legacy OBC ref slice)\n$(meta_tag)", legend=false, color=:steelblue)
        fp_dos = joinpath(out_dir, "dos_vs_energy.png")
        savefig(p_dos, fp_dos)
        add_index_row!(index_rows, meta, run_id, case_file, "dos", fp_dos)
    end

    if has_ldos_target_by_gamma
        for gamma in gammas
            gkey = gamma
            haskv(ldos_target_by_gamma, gkey) || continue
            ldos_target = getv(ldos_target_by_gamma, gkey)
            p_ldos_t = heatmap(ldos_target'; xlabel="x", ylabel="y", title="LDOS at E=$(Etarget), gamma=$(gamma)\n$(meta_tag)", aspect_ratio=1, color=:viridis, colorbar_title="rho")
            gslug = replace(value_to_str(gamma), "." => "p", "-" => "m")
            fp_lt = joinpath(out_dir, "ldos_target_gamma_$(gslug).png")
            savefig(p_ldos_t, fp_lt)
            add_index_row!(index_rows, meta, run_id, case_file, "ldos_target", fp_lt; gamma=gamma, E=Etarget)
        end
    else
        ldos_target = getv(obc, "ldos_target")
        p_ldos_t = heatmap(ldos_target'; xlabel="x", ylabel="y", title="LDOS at E=$(Etarget) (legacy OBC ref slice)\n$(meta_tag)", aspect_ratio=1, color=:viridis, colorbar_title="rho")
        fp_lt = joinpath(out_dir, "ldos_target.png")
        savefig(p_ldos_t, fp_lt)
        add_index_row!(index_rows, meta, run_id, case_file, "ldos_target", fp_lt; E=Etarget)
    end

    if has_ldos_lowest_by_gamma
        for gamma in gammas
            gkey = gamma
            haskv(ldos_lowest_by_gamma, gkey) || continue
            ldos_lowest = getv(ldos_lowest_by_gamma, gkey)
            p_ldos_l = heatmap(ldos_lowest'; xlabel="x", ylabel="y", title="LDOS lowest-|E| states gamma=$(gamma)\n$(meta_tag)", aspect_ratio=1, color=:viridis, colorbar_title="|psi|^2")
            gslug = replace(value_to_str(gamma), "." => "p", "-" => "m")
            fp_ll = joinpath(out_dir, "ldos_lowest_gamma_$(gslug).png")
            savefig(p_ldos_l, fp_ll)
            add_index_row!(index_rows, meta, run_id, case_file, "ldos_lowest", fp_ll; gamma=gamma)
        end
    else
        ldos_lowest = getv(obc, "ldos_lowest")
        p_ldos_l = heatmap(ldos_lowest'; xlabel="x", ylabel="y", title="LDOS lowest-|E| states (legacy OBC ref slice)\n$(meta_tag)", aspect_ratio=1, color=:viridis, colorbar_title="|psi|^2")
        fp_ll = joinpath(out_dir, "ldos_lowest.png")
        savefig(p_ldos_l, fp_ll)
        add_index_row!(index_rows, meta, run_id, case_file, "ldos_lowest", fp_ll)
    end
end

function plot_specloc_cuts!(index_rows, meta, run_id, case_file, out_dir, axes, specloc)
    ensure_dir(out_dir)

    gammas = getv(axes, "gammas")
    Ws = getv(axes, "Ws")
    kappas = getv(axes, "kappas")
    Es = getv(axes, "Es")

    gap = getv(specloc, "gap")
    sig = getv(specloc, "signature")

    ig0 = cld(length(gammas), 2)
    iW0 = cld(length(Ws), 2)
    ik0 = cld(length(kappas), 2)
    iE0 = cld(length(Es), 2)

    g0 = gammas[ig0]
    W0 = Ws[iW0]
    k0 = kappas[ik0]
    E0 = Es[iE0]

    # Prefer the central reference slice for stable lookup semantics, and only
    # fall back to the densest finite slice when the reference is empty.
    ref_score = count(isfinite, @view(sig[ig0, iW0, ik0, :])) + count(isfinite, @view(gap[ig0, iW0, ik0, :]))
    best_score = -1
    best_idx = (ig0, iW0, ik0)
    for gi in eachindex(gammas), wi in eachindex(Ws), ki in eachindex(kappas)
        sig_slice = @view sig[gi, wi, ki, :]
        gap_slice = @view gap[gi, wi, ki, :]
        score = count(isfinite, sig_slice) + count(isfinite, gap_slice)
        if score > best_score
            best_score = score
            best_idx = (gi, wi, ki)
        end
    end

    ig_ref, iW_ref, ik_ref = ref_score > 0 ? (ig0, iW0, ik0) : best_idx
    g_ref = gammas[ig_ref]
    W_ref = Ws[iW_ref]
    k_ref = kappas[ik_ref]

    # signature and log(gap) vs E
    sig_slice = vec(sig[ig_ref, iW_ref, ik_ref, :])
    gap_slice = vec(gap[ig_ref, iW_ref, ik_ref, :])
    sig_mask = isfinite.(sig_slice)
    gap_mask = isfinite.(gap_slice)

    # Energy reference for vertical marker line in the rendered plot.
    E_ref = E0
    E_ref_env = get(ENV, "SPECLOC_SIGGAP_E_REF", "")
    if !isempty(strip(E_ref_env))
        E_try = try
            parse(Float64, E_ref_env)
        catch
            NaN
        end
        if isfinite(E_try)
            E_ref = E_try
        end
    end

    p_sigE = plot(Es[sig_mask], sig_slice[sig_mask]; xlabel="E", ylabel="signature", title="Signature vs E (g=$(g_ref), W=$(W_ref), k=$(k_ref))", legend=false, lw=2)
    p_gapE = plot(Es[gap_mask], log10.(gap_slice[gap_mask] .+ 1e-14); xlabel="E", ylabel="log10(gap)", title="log(gap) vs E (g=$(g_ref), W=$(W_ref), k=$(k_ref))", legend=false, lw=2)
    vline!(p_sigE, [E_ref]; linestyle=:dash, color=:red, lw=1.5, alpha=0.8, label=false)
    vline!(p_gapE, [E_ref]; linestyle=:dash, color=:red, lw=1.5, alpha=0.8, label=false)

    fp_sigE = joinpath(out_dir, "specloc_signature_vs_E.png")
    fp_gapE = joinpath(out_dir, "specloc_loggap_vs_E.png")
    savefig(p_sigE, fp_sigE)
    savefig(p_gapE, fp_gapE)

    add_index_row!(index_rows, meta, run_id, case_file, "specloc_signature_vs_E", fp_sigE; gamma=g_ref, W=W_ref, kappa=k_ref)
    add_index_row!(index_rows, meta, run_id, case_file, "specloc_loggap_vs_E", fp_gapE; gamma=g_ref, W=W_ref, kappa=k_ref)

    # Full gamma-E heatmaps at fixed (W, k), preferring central reference.
    wk_ref_score = count(isfinite, @view(sig[:, iW0, ik0, :])) + count(isfinite, @view(gap[:, iW0, ik0, :]))
    wk_best_score = -1
    wk_best = (iW0, ik0)
    for wi in eachindex(Ws), ki in eachindex(kappas)
        score = count(isfinite, @view(sig[:, wi, ki, :])) + count(isfinite, @view(gap[:, wi, ki, :]))
        if score > wk_best_score
            wk_best_score = score
            wk_best = (wi, ki)
        end
    end
    iW_heat, ik_heat = wk_ref_score > 0 ? (iW0, ik0) : wk_best
    W_heat = Ws[iW_heat]
    k_heat = kappas[ik_heat]

    h_gap_gE = log10.(gap[:, iW_heat, ik_heat, :] .+ 1e-14)
    h_sig_gE = sig[:, iW_heat, ik_heat, :]
    p_gE_gap = heatmap(Es, gammas, h_gap_gE'; xlabel="E", ylabel="gamma", title="log(gap): E vs gamma (W=$(W_heat), k=$(k_heat))", colorbar_title="log10(gap)")
    p_gE_sig = heatmap(Es, gammas, h_sig_gE'; xlabel="E", ylabel="gamma", title="signature: E vs gamma (W=$(W_heat), k=$(k_heat))", colorbar_title="signature")

    fp_gE_gap = joinpath(out_dir, "specloc_gap_gamma_E.png")
    fp_gE_sig = joinpath(out_dir, "specloc_sig_gamma_E.png")
    savefig(p_gE_gap, fp_gE_gap)
    savefig(p_gE_sig, fp_gE_sig)

    add_index_row!(index_rows, meta, run_id, case_file, "specloc_gap_gamma_E", fp_gE_gap; W=W_heat, kappa=k_heat)
    add_index_row!(index_rows, meta, run_id, case_file, "specloc_sig_gamma_E", fp_gE_sig; W=W_heat, kappa=k_heat)

    # Heatmaps at fixed E0
    h_gap_gW = log10.(gap[:, :, ik0, iE0] .+ 1e-14)
    h_sig_gW = sig[:, :, ik0, iE0]
    p_gW_gap = heatmap(Ws, gammas, h_gap_gW; xlabel="W", ylabel="gamma", title="log(gap): gamma vs W (k=$(k0), E=$(E0))", colorbar_title="log10(gap)")
    p_gW_sig = heatmap(Ws, gammas, h_sig_gW; xlabel="W", ylabel="gamma", title="signature: gamma vs W (k=$(k0), E=$(E0))", colorbar_title="signature")

    fp_gW_gap = joinpath(out_dir, "specloc_gap_gamma_W.png")
    fp_gW_sig = joinpath(out_dir, "specloc_sig_gamma_W.png")
    savefig(p_gW_gap, fp_gW_gap)
    savefig(p_gW_sig, fp_gW_sig)

    add_index_row!(index_rows, meta, run_id, case_file, "specloc_gap_gamma_W", fp_gW_gap; kappa=k0, E=E0)
    add_index_row!(index_rows, meta, run_id, case_file, "specloc_sig_gamma_W", fp_gW_sig; kappa=k0, E=E0)

    h_gap_gk = log10.(gap[:, iW0, :, iE0] .+ 1e-14)
    h_sig_gk = sig[:, iW0, :, iE0]
    p_gk_gap = heatmap(kappas, gammas, h_gap_gk; xlabel="kappa", ylabel="gamma", title="log(gap): gamma vs kappa (W=$(W0), E=$(E0))", colorbar_title="log10(gap)")
    p_gk_sig = heatmap(kappas, gammas, h_sig_gk; xlabel="kappa", ylabel="gamma", title="signature: gamma vs kappa (W=$(W0), E=$(E0))", colorbar_title="signature")

    fp_gk_gap = joinpath(out_dir, "specloc_gap_gamma_kappa.png")
    fp_gk_sig = joinpath(out_dir, "specloc_sig_gamma_kappa.png")
    savefig(p_gk_gap, fp_gk_gap)
    savefig(p_gk_sig, fp_gk_sig)

    add_index_row!(index_rows, meta, run_id, case_file, "specloc_gap_gamma_kappa", fp_gk_gap; W=W0, E=E0)
    add_index_row!(index_rows, meta, run_id, case_file, "specloc_sig_gamma_kappa", fp_gk_sig; W=W0, E=E0)

    h_gap_Wk = log10.(gap[ig0, :, :, iE0] .+ 1e-14)
    h_sig_Wk = sig[ig0, :, :, iE0]
    p_Wk_gap = heatmap(kappas, Ws, h_gap_Wk; xlabel="kappa", ylabel="W", title="log(gap): W vs kappa (g=$(g0), E=$(E0))", colorbar_title="log10(gap)")
    p_Wk_sig = heatmap(kappas, Ws, h_sig_Wk; xlabel="kappa", ylabel="W", title="signature: W vs kappa (g=$(g0), E=$(E0))", colorbar_title="signature")

    fp_Wk_gap = joinpath(out_dir, "specloc_gap_W_kappa.png")
    fp_Wk_sig = joinpath(out_dir, "specloc_sig_W_kappa.png")
    savefig(p_Wk_gap, fp_Wk_gap)
    savefig(p_Wk_sig, fp_Wk_sig)

    add_index_row!(index_rows, meta, run_id, case_file, "specloc_gap_W_kappa", fp_Wk_gap; gamma=g0, E=E0)
    add_index_row!(index_rows, meta, run_id, case_file, "specloc_sig_W_kappa", fp_Wk_sig; gamma=g0, E=E0)

    # Spectrum line plots against gamma / W / kappa
    svg = getv(specloc, "spectrum_vs_gamma")
    svW = getv(specloc, "spectrum_vs_W")
    svk = getv(specloc, "spectrum_vs_kappa")

    p_svg = plot(; xlabel="gamma", ylabel="Localiser eigenvalue", title="Specloc spectrum vs gamma", legend=false)
    nonempty_pairs = [(Float64(gammas[i]), Float64.(collect(svg[i]))) for i in eachindex(svg) if !isempty(svg[i])]
    if !isempty(nonempty_pairs)
        nlev = minimum(length(last(p)) for p in nonempty_pairs)
        for (_, vals) in nonempty_pairs
            sort!(vals)
        end
        for j in 1:nlev
            xs = [p[1] for p in nonempty_pairs]
            ys = [p[2][j] for p in nonempty_pairs]
            plot!(p_svg, xs, ys; lw=0.7, alpha=0.45, color=:steelblue, label=false)
        end
    end

    ymin_env = get(ENV, "SPECLOC_SPECTRUM_YMIN", "")
    ymax_env = get(ENV, "SPECLOC_SPECTRUM_YMAX", "")
    if !isempty(strip(ymin_env)) && !isempty(strip(ymax_env))
        ymin = try parse(Float64, ymin_env) catch; NaN end
        ymax = try parse(Float64, ymax_env) catch; NaN end
        if isfinite(ymin) && isfinite(ymax) && ymin < ymax
            plot!(p_svg; ylims=(ymin, ymax))
        end
    end

    xW = Float64[]; yW = Float64[]
    for (i, vals) in enumerate(svW)
        append!(xW, fill(Ws[i], length(vals)))
        append!(yW, vals)
    end
    p_svW = scatter(xW, yW; xlabel="W", ylabel="Localiser eigenvalue", markersize=2, alpha=0.5, title="Specloc spectrum vs W", legend=false)

    xk = Float64[]; yk = Float64[]
    for (i, vals) in enumerate(svk)
        append!(xk, fill(kappas[i], length(vals)))
        append!(yk, vals)
    end
    p_svk = scatter(xk, yk; xlabel="kappa", ylabel="Localiser eigenvalue", markersize=2, alpha=0.5, title="Specloc spectrum vs kappa", legend=false)

    fp_svg = joinpath(out_dir, "specloc_spectrum_vs_gamma.png")
    fp_svW = joinpath(out_dir, "specloc_spectrum_vs_W.png")
    fp_svk = joinpath(out_dir, "specloc_spectrum_vs_kappa.png")
    savefig(p_svg, fp_svg)
    savefig(p_svW, fp_svW)
    savefig(p_svk, fp_svk)

    add_index_row!(index_rows, meta, run_id, case_file, "specloc_spectrum_vs_gamma", fp_svg; W=W0, kappa=k0, E=E0)
    add_index_row!(index_rows, meta, run_id, case_file, "specloc_spectrum_vs_W", fp_svW; gamma=g0, kappa=k0, E=E0)
    add_index_row!(index_rows, meta, run_id, case_file, "specloc_spectrum_vs_kappa", fp_svk; gamma=g0, W=W0, E=E0)
end

function process_case(case_file::String, output_root::String, run_id::String, index_rows::Vector{NamedTuple}; plot_specloc::Bool=true, return_payload::Bool=true)
    data = JLD2.load(case_file)
    result = data["result"]

    meta = result["metadata"]
    axes = result["axes"]
    specloc = result["specloc"]
    obc = result["obc"]

    case_name = splitext(basename(case_file))[1]
    out_dir = ensure_dir(joinpath(output_root, run_id, case_name))

    # Ribbon + band3d plots per gamma
    ribbon_dir = ensure_dir(joinpath(out_dir, "ribbon"))
    band3d_dir = ensure_dir(joinpath(out_dir, "band3d"))
    for gamma in axes["gammas"]
        gkey = gamma
        ribbon = result["ribbon_by_gamma"][gkey]
        band3d = result["band3d_by_gamma"][gkey]
        plot_ribbon_triplet!(index_rows, meta, run_id, case_file, ribbon_dir, gamma, ribbon)
        plot_band3d!(index_rows, meta, run_id, case_file, band3d_dir, gamma, band3d)
    end

    plot_obc!(index_rows, meta, run_id, case_file, ensure_dir(joinpath(out_dir, "obc")), axes, obc)
    if plot_specloc
        plot_specloc_cuts!(index_rows, meta, run_id, case_file, ensure_dir(joinpath(out_dir, "specloc")), axes, specloc)
    end

    if return_payload
        return (meta=meta, axes=axes, specloc=specloc)
    end

    return nothing
end

function load_case_specloc_payload(case_file::String)
    data = JLD2.load(case_file)
    result = data["result"]
    return (meta=result["metadata"], axes=result["axes"], specloc=result["specloc"])
end

function process_case_plots_only(case_file::String, output_root::String, run_id::String)
    local_rows = NamedTuple[]
    process_case(case_file, output_root, run_id, local_rows; plot_specloc=false, return_payload=false)
    return local_rows
end

function process_all_cases(input_root::String, output_root::String; run_id::String=Dates.format(now(), "yyyymmdd_HHMMSS"), verbose::Bool=true, nworkers::Int=1)
    case_files = String[]
    for (root, _, files) in walkdir(input_root)
        for f in files
            endswith(f, ".jld2") || continue
            push!(case_files, joinpath(root, f))
        end
    end
    sort!(case_files)

    isempty(case_files) && error("No .jld2 files found under $(input_root)")

    if verbose
        println("Found $(length(case_files)) case files to process.")
    end

    index_rows = NamedTuple[]
    specloc_groups = Dict{Any, Any}()

    parallel_enabled = nworkers > 1 && Distributed.nprocs() >= nworkers

    if parallel_enabled
        proc_ids = Distributed.procs()[1:nworkers]
        pool = Distributed.WorkerPool(proc_ids)
        if verbose
            println("Generating per-case plots in parallel using $(length(proc_ids)) Julia processes...")
        end

        row_batches = pmap(pool, case_files) do cf
            process_case_plots_only(cf, output_root, run_id)
        end

        for rows in row_batches
            append!(index_rows, rows)
        end
    else
        if nworkers > 1 && verbose
            println("Parallel mode requested (nworkers=$(nworkers)) but not enough Julia processes available; falling back to serial.")
        end
        prog = verbose ? nothing : Progress(length(case_files); desc="Processing", dt=0.5)
        for cf in case_files
            if verbose
                println("Processing: $(cf)")
            else
                next!(prog)
            end
            process_case(cf, output_root, run_id, index_rows; plot_specloc=false, return_payload=false)
        end
    end

    if verbose
        println("Collecting and aggregating specloc data across all files...")
    end
    prog_collect = verbose ? nothing : Progress(length(case_files); desc="Collecting specloc", dt=0.5)
    for cf in case_files
        payload = load_case_specloc_payload(cf)
        rel = relpath(cf, input_root)
        parts = splitpath(rel)
        run_group = isempty(parts) ? "root" : parts[1]
        add_specloc_case!(specloc_groups, cf, run_group, payload.meta, payload.axes, payload.specloc)
        if !verbose
            next!(prog_collect)
        end
    end

    if verbose
        println("Building aggregated specloc plots across chunked files...")
    end
    
    prog_agg = verbose ? nothing : Progress(length(specloc_groups); desc="Aggregating specloc", dt=0.5)
    
    for (key, grp) in specloc_groups
        out_dir = ensure_dir(joinpath(output_root, run_id, "specloc_aggregated", group_slug(key)))
        axes_agg, specloc_agg = build_aggregated_specloc_payload(grp)
        case_label = "aggregated::" * string(length(grp[:case_files])) * " files"
        plot_specloc_cuts!(index_rows, grp[:meta], run_id, case_label, out_dir, axes_agg, specloc_agg)
        if !verbose
            next!(prog_agg)
        end
    end

    index_file = joinpath(output_root, run_id, "plot_index.tsv")
    save_plot_index_tsv(index_rows, index_file)
    
    if verbose
        println("Saved plot index: $(index_file)")
    end
    
    return index_file
end

end # module
