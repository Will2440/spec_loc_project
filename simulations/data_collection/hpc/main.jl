using Dates
using DelimitedFiles
using Printf
using JLD2: @save

project_root = @__DIR__
include(joinpath(project_root, "solver.jl"))
using .SpecLocSolver

parse_lit(s::AbstractString) = Base.invokelatest(eval, Meta.parse(s))

# =====================================================================
# Dynamic energy range estimation
# =====================================================================

function scalar_shift_energy(kx, ky, gamma, pt::Symbol)
    pt == :sym_cos_sum  || pt == :symmetric  && return gamma*(cos(kx)+cos(ky))
    pt == :sym_cos_diff  && return gamma*(cos(kx)-cos(ky))
    pt == :sym_cos_add   && return gamma*cos(kx+ky)
    pt == :sym_cos_sub   && return gamma*cos(kx-ky)
    pt == :asym_sin_sum  && return gamma*(sin(kx)+sin(ky))
    pt == :asym_sin_diff && return gamma*(sin(kx)-sin(ky))
    pt == :asym_sin_add  && return gamma*sin(kx+ky)
    pt == :asym_sin_sub  && return gamma*sin(kx-ky)
    pt == :tilt          && return -gamma*sin(kx)
    pt == :none          && return 0.0
    error("Unknown perturbation type: $pt")
end

function estimate_band_extrema(A, B, m, gamma, pt::Symbol; nk=81)
    ks = collect(range(-π, π; length=nk+1))[1:end-1]
    emin, emax = Inf, -Inf
    @inbounds for kx in ks, ky in ks
        dz = m + 2B - B*cos(kx) - B*cos(ky)
        dn = sqrt((A*sin(kx))^2 + (A*sin(ky))^2 + dz^2)
        sc = scalar_shift_energy(kx, ky, gamma, pt)
        e1 = sc-dn; e2 = sc+dn
        e1 < emin && (emin=e1);  e2 > emax && (emax=e2)
    end
    return emin, emax
end

function compute_energy_vals(
    A, B, m, gamma, pt::Symbol,
    mode::Symbol, npts::Int, margin::Float64, nk::Int,
    fixed::Vector{Float64},
)
    mode == :fixed && return fixed
    mode == :dynamic_band || error("Unknown energy_range_mode: $mode")
    lo, hi = estimate_band_extrema(Float64(A), Float64(B), Float64(m), Float64(gamma), pt; nk=nk)
    w = hi - lo;  w > 0 || (w = max(abs(hi), 1.0))
    pad = margin * w
    return collect(range(lo-pad, hi+pad; length=npts))
end

# =====================================================================
# Parameter file reading (23 columns)
# =====================================================================

function read_parameter_rows(path::String)
    raw   = DelimitedFiles.readdlm(path, '\t', String; header=true)
    table = raw[1]
    rows  = NamedTuple[]
    for i in 1:size(table, 1)
        r = table[i, :]
        push!(rows, (
            As=Vector{Float64}(parse_lit(r[1])),
            Bs=Vector{Float64}(parse_lit(r[2])),
            ms=Vector{Float64}(parse_lit(r[3])),
            perturbation_type=Symbol(parse_lit(r[4])),
            disorder_type=Symbol(parse_lit(r[5])),
            Lx_obc=Int(parse_lit(r[6])),
            Ly_obc=Int(parse_lit(r[7])),
            gamma_vals=Vector{Float64}(parse_lit(r[8])),
            W_vals=Vector{Float64}(parse_lit(r[9])),
            kappa_vals=Vector{Float64}(parse_lit(r[10])),
            energy_range_mode=Symbol(parse_lit(r[11])),
            energy_points=Int(parse_lit(r[12])),
            energy_margin_fraction=Float64(parse_lit(r[13])),
            energy_scan_nk=Int(parse_lit(r[14])),
            fixed_E_vals=Vector{Float64}(parse_lit(r[15])),
            specloc_x=Int(parse_lit(r[16])),
            specloc_y=Int(parse_lit(r[17])),
            seed=Int(parse_lit(r[18])),
            orbital_displacement=Float64(parse_lit(r[19])),
            phi=Float64(parse_lit(r[20])),
            n_disorder_realisations=Int(parse_lit(r[21])),
            scale_kappa_to_L=Bool(parse_lit(r[22])),
            kappa_scales=Vector{Float64}(parse_lit(r[23])),
        ))
    end
    return rows
end

function newest_params_file()
    dir   = joinpath(project_root, "param_sets")
    files = filter(f -> endswith(f, ".dat"), readdir(dir; join=true))
    isempty(files) && error("No .dat files found in $dir. Run param_prep.jl first.")
    sort!(files, by=f -> stat(f).mtime)
    return files[end]
end

resolve_run_id(ri) = let e = strip(get(ENV,"SPECLOC_RUN_ID",""))
    if !isempty(e)
        return e
    else
        # Include SLURM job ID if running on HPC
        job_id = strip(get(ENV,"SLURM_ARRAY_JOB_ID",""))
        array_task = strip(get(ENV,"SLURM_ARRAY_TASK_ID",""))
        ts = Dates.format(now(),"yyyymmdd_HHMMSS")
        if !isempty(job_id)
            return "job_$(job_id)_task_$(array_task)_$(ts)_r$ri"
        else
            return ts*"_r$ri"
        end
    end
end

format_dur(s) = let t=max(0,round(Int,s)); @sprintf("%02d:%02d:%02d",t÷3600,(t%3600)÷60,t%60) end

stdout_tty() = try Base.isatty(stdout) catch; false end

# =====================================================================
# Row runner
# =====================================================================

function run_row(row_index::Int, params_path::String)
    rows = read_parameter_rows(params_path)
    1 <= row_index <= length(rows) || error("row_index=$row_index out of range 1–$(length(rows))")
    row = rows[row_index]

    run_id  = resolve_run_id(row_index)
    out_dir = joinpath(project_root, "results", run_id)
    isdir(out_dir) || mkpath(out_dir)

    total = length(row.As)*length(row.Bs)*length(row.ms)*length(row.gamma_vals)
    println("Row $row_index: chunks=$total  pt=$(row.perturbation_type)  dt=$(row.disorder_type)  n_real=$(row.n_disorder_realisations)  scale_κ=$(row.scale_kappa_to_L)  d=$(row.orbital_displacement)  φ=$(row.phi)")
    flush(stdout)

    ci = 0;  t0 = time();  tlast = t0;  tcum = 0.0

    for A in row.As, B in row.Bs, m in row.ms, gamma in row.gamma_vals
        ci += 1

        E_vals = compute_energy_vals(
            A, B, m, gamma, row.perturbation_type,
            row.energy_range_mode, row.energy_points,
            row.energy_margin_fraction, row.energy_scan_nk, row.fixed_E_vals)

        cfg = SolverCaseConfig(
            A=A, B=B, m=m,
            perturbation_type=row.perturbation_type,
            disorder_type=row.disorder_type,
            Lx_obc=row.Lx_obc, Ly_obc=row.Ly_obc,
            gamma_vals=[Float64(gamma)],
            W_vals=row.W_vals, kappa_vals=row.kappa_vals,
            E_vals=E_vals,
            specloc_x=row.specloc_x, specloc_y=row.specloc_y,
            seed=row.seed,
            orbital_displacement=row.orbital_displacement,
            phi=row.phi,
            n_disorder_realisations=row.n_disorder_realisations,
            scale_kappa_to_L=row.scale_kappa_to_L,
            kappa_scales=row.kappa_scales,
        )

        result = run_case(cfg)
        fname  = joinpath(out_dir,
            "row$(row_index)_chunk$(ci)_A$(A)_B$(B)_m$(m)_g$(gamma)" *
            "_pt$(row.perturbation_type)_dt$(row.disorder_type)" *
            "_d$(row.orbital_displacement)_phi$(row.phi)" *
            "_sx$(row.specloc_x)_sy$(row.specloc_y).jld2")
        @save fname result

        dt = time()-tlast;  tcum += dt;  eta = (total-ci)*(tcum/ci)
        line = @sprintf("[Row %d] %d/%d | dt=%.1fs | avg=%.1fs | ETA=%s",
                        row_index, ci, total, dt, tcum/ci, format_dur(eta))
        stdout_tty() ? print("\r"*line) : println(line)
        flush(stdout);  tlast = time()
    end

    stdout_tty() && println()
    ttot = time()-t0
    println(@sprintf("[Row %d] done %d/%d in %s | avg=%.1fs", row_index, ci, total, format_dur(ttot), ttot/max(ci,1)))
end

# =====================================================================
# Entry point
# =====================================================================

function main()
    length(ARGS) >= 1 || (println("Usage: julia main.jl <row_index> [params_file]"); return)
    params_path = length(ARGS) >= 2 ? ARGS[2] : newest_params_file()
    println("Params: $params_path")
    run_row(parse(Int, ARGS[1]), params_path)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end