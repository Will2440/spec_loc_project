using Dates
using DelimitedFiles
using Printf
using JLD2: @save

project_root = @__DIR__
include(joinpath(project_root, "solver.jl"))
using .SpecLocSolver

parse_lit(s::AbstractString) = Base.invokelatest(eval, Meta.parse(s))

# =====================================================================
# Per-case energy computation (mirrors param_prep.jl helpers)
# =====================================================================

function scalar_shift_energy(kx, ky, gamma, perturbation_type::Symbol)
    if perturbation_type == :symmetric
        return gamma * (cos(kx) + cos(ky))
    elseif perturbation_type == :tilt
        return gamma * sin(kx)
    elseif perturbation_type == :none
        return 0.0
    else
        error("Unknown perturbation type: $perturbation_type")
    end
end

function estimate_band_extrema_case(
    A::Float64,
    B::Float64,
    m::Float64,
    gamma::Float64,
    B_y::Float64,
    perturbation_type::Symbol;
    nk::Int=81,
)
    kx_vals = collect(range(-π, π; length=nk + 1))[1:end-1]
    ky_vals = collect(range(-π, π; length=nk + 1))[1:end-1]

    e_min = Inf
    e_max = -Inf

    @inbounds for kx in kx_vals, ky in ky_vals
        dx = A * sin(kx)
        dy = A * sin(ky)
        dz = m + 2B - B * cos(kx) - B * B_y * cos(ky)
        scalar = scalar_shift_energy(kx, ky, gamma, perturbation_type)
        dnorm = sqrt(dx * dx + dy * dy + dz * dz)
        e1 = scalar - dnorm
        e2 = scalar + dnorm
        if e1 < e_min
            e_min = e1
        end
        if e2 > e_max
            e_max = e2
        end
    end

    return e_min, e_max
end

function compute_energy_vals_for_case(
    A::Float64,
    B::Float64,
    m::Float64,
    gamma::Float64,
    B_y::Float64,
    perturbation_type::Symbol,
    energy_range_mode::Symbol,
    energy_points::Int,
    energy_margin_fraction::Float64,
    energy_scan_nk::Int,
    fixed_E_vals::Vector{Float64},
)
    if energy_range_mode == :fixed
        return fixed_E_vals
    elseif energy_range_mode == :dynamic_band
        lo, hi = estimate_band_extrema_case(A, B, m, gamma, B_y, perturbation_type; nk=energy_scan_nk)
        width = hi - lo
        if !(width > 0)
            width = max(abs(hi), 1.0)
        end
        pad = energy_margin_fraction * width
        return collect(range(lo - pad, hi + pad; length=energy_points))
    else
        error("Unknown energy_range_mode: $energy_range_mode")
    end
end

function read_parameter_rows(params_path::String)
    raw = DelimitedFiles.readdlm(params_path, '\t', String; header=true)
    table = raw[1]
    rows = NamedTuple[]

    for i in 1:size(table, 1)
        r = table[i, :]
        push!(rows, (
            As=Vector{Float64}(parse_lit(r[1])),
            Bs=Vector{Float64}(parse_lit(r[2])),
            ms=Vector{Float64}(parse_lit(r[3])),
            B_y=Float64(parse_lit(r[4])),
            perturbation_type=Symbol(parse_lit(r[5])),
            disorder_type=Symbol(parse_lit(r[6])),
            winding_number=Int(parse_lit(r[7])),
            Lx_ribbon=Int(parse_lit(r[8])),
            Lx_obc=Int(parse_lit(r[9])),
            Ly_obc=Int(parse_lit(r[10])),
            gamma_vals=Vector{Float64}(parse_lit(r[11])),
            W_vals=Vector{Float64}(parse_lit(r[12])),
            kappa_vals=Vector{Float64}(parse_lit(r[13])),
            energy_range_mode=Symbol(parse_lit(r[14])),
            energy_points=Int(parse_lit(r[15])),
            energy_margin_fraction=Float64(parse_lit(r[16])),
            energy_scan_nk=Int(parse_lit(r[17])),
            fixed_E_vals=Vector{Float64}(parse_lit(r[18])),
            N_ky=Int(parse_lit(r[19])),
            Nkx_bulk=Int(parse_lit(r[20])),
            Nky_bulk=Int(parse_lit(r[21])),
            Nkx_3d=Int(parse_lit(r[22])),
            Nky_3d=Int(parse_lit(r[23])),
            ldos_target_E=Float64(parse_lit(r[24])),
            ldos_eta=Float64(parse_lit(r[25])),
            dos_eta=Float64(parse_lit(r[26])),
            lowest_energy_count=Int(parse_lit(r[27])),
            specloc_x=Int(parse_lit(r[28])),
            specloc_y=Int(parse_lit(r[29])),
            seed=Int(parse_lit(r[30])),
        ))
    end
    return rows
end

function newest_params_file()
    params_dir = joinpath(project_root, "param_sets")
    files = filter(f -> endswith(f, ".dat"), readdir(params_dir; join=true))
    isempty(files) && error("No parameter .dat files found in $(params_dir). Run parameter.jl first.")
    sort!(files, by=f -> stat(f).mtime)
    return files[end]
end

function resolve_run_id(row_index::Int)
    run_id_env = strip(get(ENV, "SPECLOC_RUN_ID", ""))
    if !isempty(run_id_env)
        return run_id_env
    end
    return Dates.format(now(), "yyyymmdd_HHMMSS") * "_r$(row_index)"
end

function format_duration(seconds::Real)
    total = max(0, round(Int, seconds))
    h = total ÷ 3600
    m = (total % 3600) ÷ 60
    s = total % 60
    return @sprintf("%02d:%02d:%02d", h, m, s)
end

function stdout_is_tty()
    if isdefined(Base, :isatty)
        return Base.isatty(stdout)
    elseif isdefined(Base, :Libc) && isdefined(Base.Libc, :isatty)
        return Base.Libc.isatty(1) != 0
    else
        return false
    end
end

function run_row(row_index::Int, params_path::String)
    rows = read_parameter_rows(params_path)
    row_index < 1 && error("row_index must be >= 1")
    row_index > length(rows) && error("row_index=$(row_index) exceeds available rows=$(length(rows))")

    row = rows[row_index]

    run_id = resolve_run_id(row_index)
    out_dir = joinpath(project_root, "results", run_id)
    isdir(out_dir) || mkpath(out_dir)

    total_chunks = length(row.As) * length(row.Bs) * length(row.ms) * length(row.gamma_vals)
    println("Row $(row_index): total chunks=$(total_chunks)")

    chunk_index = 0
    row_start_time = time()
    last_chunk_time = row_start_time
    cumulative_chunk_time = 0.0

    for A in row.As, B in row.Bs, m in row.ms, gamma in row.gamma_vals
        chunk_index += 1
        # Compute energy range per individual (A,B,m,gamma) case.
        E_vals = compute_energy_vals_for_case(
            Float64(A),
            Float64(B),
            Float64(m),
            Float64(gamma),
            Float64(row.B_y),
            row.perturbation_type,
            row.energy_range_mode,
            row.energy_points,
            row.energy_margin_fraction,
            row.energy_scan_nk,
            row.fixed_E_vals,
        )

        cfg = SolverCaseConfig(
            A=A,
            B=B,
            m=m,
            B_y=row.B_y,
            perturbation_type=row.perturbation_type,
            disorder_type=row.disorder_type,
            winding_number=row.winding_number,
            Lx_ribbon=row.Lx_ribbon,
            Lx_obc=row.Lx_obc,
            Ly_obc=row.Ly_obc,
            gamma_vals=[Float64(gamma)],
            W_vals=row.W_vals,
            kappa_vals=row.kappa_vals,
            E_vals=E_vals,
            N_ky=row.N_ky,
            Nkx_bulk=row.Nkx_bulk,
            Nky_bulk=row.Nky_bulk,
            Nkx_3d=row.Nkx_3d,
            Nky_3d=row.Nky_3d,
            ldos_target_E=row.ldos_target_E,
            ldos_eta=row.ldos_eta,
            dos_eta=row.dos_eta,
            lowest_energy_count=row.lowest_energy_count,
            specloc_x=row.specloc_x,
            specloc_y=row.specloc_y,
            seed=row.seed,
        )

        result = run_case(cfg)
        filename = joinpath(
            out_dir,
            "row$(row_index)_chunk$(chunk_index)_A$(A)_B$(B)_m$(m)_g$(gamma)_pt$(row.perturbation_type)_dt$(row.disorder_type)_sx$(row.specloc_x)_sy$(row.specloc_y).jld2",
        )
        @save filename result

        now_t = time()
        chunk_time = now_t - last_chunk_time
        cumulative_chunk_time += chunk_time
        avg_chunk_time = cumulative_chunk_time / chunk_index
        remaining_chunks = total_chunks - chunk_index
        eta_seconds = remaining_chunks * avg_chunk_time

        progress_line = @sprintf(
            "[Row %d] chunk %d/%d | dt=%.1fs | avg=%.1fs | ETA=%s",
            row_index,
            chunk_index,
            total_chunks,
            chunk_time,
            avg_chunk_time,
            format_duration(eta_seconds),
        )

        if stdout_is_tty()
            print("\r" * progress_line)
        else
            println(progress_line)
        end
        flush(stdout)

        last_chunk_time = now_t
    end

    if stdout_is_tty()
        println()
    end

    total_time = time() - row_start_time
    println(@sprintf("[Row %d] completed %d/%d chunks in %s | avg=%.1fs", row_index, chunk_index, total_chunks, format_duration(total_time), total_time / max(chunk_index, 1)))
end

function main()
    if length(ARGS) < 1
        println("Usage:")
        println("  julia main.jl <row_index> [params_file]")
        return
    end

    row_index = parse(Int, ARGS[1])
    params_path = length(ARGS) >= 2 ? ARGS[2] : newest_params_file()
    println("Using params file: $(params_path)")
    run_row(row_index, params_path)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
