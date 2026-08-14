module SpecLocPlotting

using JLD2
using Plots
using Printf
using Dates

export process_all_cases

getv(x::AbstractDict, key::String) = x[key]
getv(x::AbstractDict, key::Symbol) = x[key]
getv(x::NamedTuple, key::String) = getproperty(x, Symbol(key))
getv(x::NamedTuple, key::Symbol) = getproperty(x, key)

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

function plot_ribbon_triplet!(index_rows, meta, run_id, case_file, out_dir, gamma, ribbon)
    ensure_dir(out_dir)

    ky = getv(ribbon, "ky_vals")
    energies = getv(ribbon, "energies")
    ipr = getv(ribbon, "iprs")
    chd = getv(ribbon, "chern_contribution_density_on_ribbon")
    accum_e = getv(ribbon, "accumulation_energies")
    accum_c = getv(ribbon, "cumulative_chern")

    p1 = plot(title="Ribbon IPR gamma=$(gamma)", xlabel="k_y", ylabel="Energy", legend=false, colorbar=true, colorbar_title="IPR")
    for b in 1:size(energies, 1)
        plot!(p1, ky, @view(energies[b, :]); line_z=@view(ipr[b, :]), lw=1.0, c=:viridis, label=false)
    end

    chlim = maximum(abs.(chd))
    p2 = plot(title="Ribbon dC/dE gamma=$(gamma)", xlabel="k_y", ylabel="Energy", legend=false, colorbar=true, colorbar_title="dC/dE", clims=(-chlim, chlim))
    for b in 1:size(energies, 1)
        plot!(p2, ky, @view(energies[b, :]); line_z=@view(chd[b, :]), lw=1.0, c=:RdBu, label=false)
    end

    p3 = plot(accum_c, accum_e; xlabel="Accumulated Chern", ylabel="Energy", title="Chern accumulation gamma=$(gamma)", lw=2, color=:darkred, legend=false)
    vline!(p3, [0.0]; linestyle=:dash, color=:black, lw=1, label=false)
    hline!(p3, [0.0]; linestyle=:dash, color=:black, lw=1, label=false)

    fp1 = joinpath(out_dir, "ribbon_ipr_gamma$(gamma).png")
    fp2 = joinpath(out_dir, "ribbon_dcdE_gamma$(gamma).png")
    fp3 = joinpath(out_dir, "ribbon_chern_acc_gamma$(gamma).png")

    savefig(p1, fp1)
    savefig(p2, fp2)
    savefig(p3, fp3)

    add_index_row!(index_rows, meta, run_id, case_file, "ribbon_ipr", fp1; gamma=gamma)
    add_index_row!(index_rows, meta, run_id, case_file, "ribbon_dcdE", fp2; gamma=gamma)
    add_index_row!(index_rows, meta, run_id, case_file, "ribbon_chern_acc", fp3; gamma=gamma)
end

function plot_band3d!(index_rows, meta, run_id, case_file, out_dir, gamma, band3d)
    ensure_dir(out_dir)
    kx = getv(band3d, "kx_vals")
    ky = getv(band3d, "ky_vals")
    em = getv(band3d, "e_minus")
    ep = getv(band3d, "e_plus")

    plt = surface(kx, ky, em'; xlabel="k_x", ylabel="k_y", zlabel="E", color=:viridis, alpha=0.8, legend=false, title="3D Bandstructure gamma=$(gamma)")
    surface!(plt, kx, ky, ep'; color=:plasma, alpha=0.8, legend=false)

    fp = joinpath(out_dir, "band3d_gamma$(gamma).png")
    savefig(plt, fp)
    add_index_row!(index_rows, meta, run_id, case_file, "band3d", fp; gamma=gamma)
end

function plot_obc!(index_rows, meta, run_id, case_file, out_dir, axes, obc)
    ensure_dir(out_dir)

    E = getv(axes, "Es")
    dos = getv(obc, "dos")
    ldos_target = getv(obc, "ldos_target")
    ldos_lowest = getv(obc, "ldos_lowest")
    Etarget = getv(obc, "ldos_target_E")

    p_dos = plot(E, dos; xlabel="Energy", ylabel="DOS", title="Total DOS", legend=false, color=:steelblue)
    p_ldos_t = heatmap(ldos_target'; xlabel="x", ylabel="y", title="LDOS at E=$(Etarget)", aspect_ratio=1, color=:viridis, colorbar_title="rho")
    p_ldos_l = heatmap(ldos_lowest'; xlabel="x", ylabel="y", title="LDOS lowest-|E| states", aspect_ratio=1, color=:viridis, colorbar_title="|psi|^2")

    fp_dos = joinpath(out_dir, "dos_vs_energy.png")
    fp_lt = joinpath(out_dir, "ldos_target.png")
    fp_ll = joinpath(out_dir, "ldos_lowest.png")

    savefig(p_dos, fp_dos)
    savefig(p_ldos_t, fp_lt)
    savefig(p_ldos_l, fp_ll)

    add_index_row!(index_rows, meta, run_id, case_file, "dos", fp_dos)
    add_index_row!(index_rows, meta, run_id, case_file, "ldos_target", fp_lt; E=Etarget)
    add_index_row!(index_rows, meta, run_id, case_file, "ldos_lowest", fp_ll)
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

    # signature and log(gap) vs E
    p_sigE = plot(Es, vec(sig[ig0, iW0, ik0, :]); xlabel="E", ylabel="signature", title="Signature vs E (g=$(g0), W=$(W0), k=$(k0))", legend=false, lw=2)
    p_gapE = plot(Es, log10.(vec(gap[ig0, iW0, ik0, :]) .+ 1e-14); xlabel="E", ylabel="log10(gap)", title="log(gap) vs E (g=$(g0), W=$(W0), k=$(k0))", legend=false, lw=2)

    fp_sigE = joinpath(out_dir, "specloc_signature_vs_E.png")
    fp_gapE = joinpath(out_dir, "specloc_loggap_vs_E.png")
    savefig(p_sigE, fp_sigE)
    savefig(p_gapE, fp_gapE)

    add_index_row!(index_rows, meta, run_id, case_file, "specloc_signature_vs_E", fp_sigE; gamma=g0, W=W0, kappa=k0)
    add_index_row!(index_rows, meta, run_id, case_file, "specloc_loggap_vs_E", fp_gapE; gamma=g0, W=W0, kappa=k0)

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

    xg = Float64[]; yg = Float64[]
    for (i, vals) in enumerate(svg)
        append!(xg, fill(gammas[i], length(vals)))
        append!(yg, vals)
    end
    p_svg = scatter(xg, yg; xlabel="gamma", ylabel="Localiser eigenvalue", markersize=2, alpha=0.5, title="Specloc spectrum vs gamma", legend=false)

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

function process_case(case_file::String, output_root::String, run_id::String, index_rows::Vector{NamedTuple})
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
    plot_specloc_cuts!(index_rows, meta, run_id, case_file, ensure_dir(joinpath(out_dir, "specloc")), axes, specloc)
end

function process_all_cases(input_root::String, output_root::String; run_id::String=Dates.format(now(), "yyyymmdd_HHMMSS"))
    case_files = String[]
    for (root, _, files) in walkdir(input_root)
        for f in files
            endswith(f, ".jld2") || continue
            push!(case_files, joinpath(root, f))
        end
    end

    isempty(case_files) && error("No .jld2 files found under $(input_root)")

    println("Found $(length(case_files)) case files to process.")

    index_rows = NamedTuple[]
    for cf in case_files
        println("Processing: $(cf)")
        process_case(cf, output_root, run_id, index_rows)
    end

    index_file = joinpath(output_root, run_id, "plot_index.tsv")
    save_plot_index_tsv(index_rows, index_file)
    println("Saved plot index: $(index_file)")
    return index_file
end

end # module
