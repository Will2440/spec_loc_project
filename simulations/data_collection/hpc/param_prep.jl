using DelimitedFiles
using Dates
using Printf
using Statistics

project_root = @__DIR__
outdir = joinpath(project_root, "param_sets")
isdir(outdir) || mkpath(outdir)

# =====================================================================
# User config (EDIT THESE)
# =====================================================================

As = collect(0.5:0.5:3.0)
Bs = collect(0.5:0.5:3.0)
ms = collect(-5.0:0.5:1.0) #[-1.0, -2.0]
B_ys = [1.0]

perturbation_types = [:symmetric, :tilt]  # :none, :symmetric, :tilt
disorder_types = [:anderson] # :none, :anderson, :mass
boundary_cond_twists = [0]

Lx_ribbons = [50]
Lx_obcs = [20]
Ly_obcs = [20]

gamma_vals = collect(-3.0:0.5:3.0)
W_vals = collect(0.0:0.5:4.0)
kappa_vals = [2e-1]

# Energy range mode
energy_range_mode = :dynamic_band  # :fixed, :dynamic_band
fixed_E_vals = collect(range(-5.0, 5.0; length=101))
energy_points = 51
energy_margin_fraction = 0.10  # 10% padding above/below estimated extrema
energy_scan_nk = 11            # lightweight band extrema scan grid per k-axis
# When true, dynamic-band jobs are written at per-case granularity for
# (A,B,m,gamma), improving HPC parallel scaling and avoiding bundled cases.
dynamic_energy_split_cases = true

# Grid resolutions for downstream calculations
N_ky = 101 ## for ribbon geometry
Nkx_bulk = 81 ## for 2D bulk geometry
Nky_bulk = 81 ## for 2D bulk geometry
Nkx_3d = 101 ## for 3D bandstructure plot
Nky_3d = 101 ## for 3D bandstructure plot

# DOS/LDOS settings
ldos_target_E = 0.0
ldos_eta = 0.1
dos_eta = 0.1
lowest_energy_count = 4

# Spectral localiser placement
# Supported modes:
#   :manual (calculate at a single user-specified point)
#   :all
#   :all_plus_outside
#   :centre_point
#   :centre_region
#   :mid_edge_point_x, :mid_edge_point_y
#   :mid_edge_region_x, :mid_edge_region_y
#   :corner_point_00, :corner_point_01, :corner_point_10, :corner_point_11
#   :corner_region_00, :corner_region_01, :corner_region_10, :corner_region_11
specloc_mode = :centre_point
specloc_x = 7
specloc_y = 7
specloc_resolution_x = 1
specloc_resolution_y = 1
specloc_outside_fraction = 0.10
specloc_centre_fraction = 0.30
specloc_mid_edge_fraction = 0.30
specloc_corner_fraction = 0.20

# Reproducibility
seed = 1234

# ---------------------------------------------------------------------
# Row allocation strategy
# ---------------------------------------------------------------------
# :explicit_vals_per_row  -> you set *_vals_per_row directly
# :target_rows            -> auto-select chunk counts to approach target_number_of_rows
allocation_mode = :target_rows

target_number_of_rows = 200

# In :target_rows mode we prioritize reaching the requested output row count.
# Per-case dynamic splitting can dominate row count and make targeting ineffective,
# so it is disabled automatically in that mode.
dynamic_energy_split_cases_active = dynamic_energy_split_cases
dynamic_split_disabled_for_target_rows = false
if allocation_mode == :target_rows && dynamic_energy_split_cases_active
    dynamic_energy_split_cases_active = false
    dynamic_split_disabled_for_target_rows = true
end

# Used when allocation_mode == :explicit_vals_per_row
A_vals_per_row = max(1, length(As))
B_vals_per_row = max(1, length(Bs))
m_vals_per_row = max(1, length(ms))
gamma_vals_per_row = max(1, length(gamma_vals))
W_vals_per_row = max(1, length(W_vals))
kappa_vals_per_row = max(1, length(kappa_vals))

# Split policies
split_mode = :contiguous  # :contiguous or :roundrobin

# =====================================================================
# Helpers
# =====================================================================

function float3(x)
    return replace(string(x), "." => "p")
end

function range_info(vec)
    isempty(vec) && return "none"
    if eltype(vec) <: Integer
        return "$(vec[1])-$(vec[end])-$(length(vec))"
    end
    return string(float3(vec[1]), "-", float3(vec[end]), "-", length(vec))
end

function chunk_vector(v::Vector, vals_per_chunk::Int; mode::Symbol=:contiguous)
    isempty(v) && return [v]
    vals_per_chunk = max(1, vals_per_chunk)

    n_chunks = ceil(Int, length(v) / vals_per_chunk)
    if mode == :contiguous
        chunks = Vector{typeof(v)}()
        i = 1
        while i <= length(v)
            j = min(length(v), i + vals_per_chunk - 1)
            push!(chunks, v[i:j])
            i = j + 1
        end
        return chunks
    elseif mode == :roundrobin
        chunks = [typeof(v)() for _ in 1:n_chunks]
        for (i, x) in enumerate(v)
            push!(chunks[((i - 1) % n_chunks) + 1], x)
        end
        return chunks
    else
        error("Unknown split mode: $mode")
    end
end

function stride_axis(lo::Int, hi::Int, step::Int)
    step = max(1, step)
    if lo > hi
        return Int[]
    end
    vals = collect(lo:step:hi)
    if isempty(vals) || vals[end] != hi
        push!(vals, hi)
    end
    return vals
end

function bounded_fraction(x::Real)
    return clamp(Float64(x), 0.0, 1.0)
end

function center_index(L::Int)
    return cld(L, 2)
end

function center_interval(L::Int, frac::Real)
    f = bounded_fraction(frac)
    c = center_index(L)
    half = floor(Int, 0.5 * f * L)
    lo = max(1, c - half)
    hi = min(L, c + half)
    return lo, hi
end

function edge_span_interval(L::Int, frac::Real, side_low::Bool)
    f = bounded_fraction(frac)
    span = max(1, ceil(Int, f * L))
    if side_low
        return 1, min(L, span)
    end
    return max(1, L - span + 1), L
end

function parse_corner_suffix(mode::Symbol, prefix::String)
    s = String(mode)
    startswith(s, prefix) || return nothing
    suffix = s[length(prefix) + 1:end]
    length(suffix) == 2 || return nothing
    xbit, ybit = suffix[1], suffix[2]
    (xbit in ('0', '1') && ybit in ('0', '1')) || return nothing
    return (xbit == '0', ybit == '0')
end

function uniquify_points(points::Vector{Tuple{Int, Int}})
    return sort(unique(points), by=p -> (p[1], p[2]))
end

function build_specloc_points(
    mode::Symbol,
    Lx::Int,
    Ly::Int;
    base_x::Int,
    base_y::Int,
    res_x::Int,
    res_y::Int,
    outside_frac::Real,
    centre_frac::Real,
    mid_edge_frac::Real,
    corner_frac::Real,
)
    rx = max(1, res_x)
    ry = max(1, res_y)

    if mode == :manual
        return [(base_x, base_y)]
    elseif mode == :all
        xs = stride_axis(1, Lx, rx)
        ys = stride_axis(1, Ly, ry)
        return [(x, y) for x in xs for y in ys]
    elseif mode == :all_plus_outside
        ox = ceil(Int, bounded_fraction(outside_frac) * Lx)
        oy = ceil(Int, bounded_fraction(outside_frac) * Ly)
        xs = stride_axis(1 - ox, Lx + ox, rx)
        ys = stride_axis(1 - oy, Ly + oy, ry)
        return [(x, y) for x in xs for y in ys]
    elseif mode == :centre_point
        return [(center_index(Lx), center_index(Ly))]
    elseif mode == :centre_region
        xlo, xhi = center_interval(Lx, centre_frac)
        ylo, yhi = center_interval(Ly, centre_frac)
        xs = stride_axis(xlo, xhi, rx)
        ys = stride_axis(ylo, yhi, ry)
        return [(x, y) for x in xs for y in ys]
    elseif mode == :mid_edge_point_x
        return [(1, center_index(Ly)), (Lx, center_index(Ly))]
    elseif mode == :mid_edge_point_y
        return [(center_index(Lx), 1), (center_index(Lx), Ly)]
    elseif mode == :mid_edge_region_x
        ylo, yhi = center_interval(Ly, mid_edge_frac)
        ys = stride_axis(ylo, yhi, ry)
        return vcat([(1, y) for y in ys], [(Lx, y) for y in ys])
    elseif mode == :mid_edge_region_y
        xlo, xhi = center_interval(Lx, mid_edge_frac)
        xs = stride_axis(xlo, xhi, rx)
        return vcat([(x, 1) for x in xs], [(x, Ly) for x in xs])
    end

    corner_pt = parse_corner_suffix(mode, "corner_point_")
    if corner_pt !== nothing
        low_x, low_y = corner_pt
        x = low_x ? 1 : Lx
        y = low_y ? 1 : Ly
        return [(x, y)]
    end

    corner_reg = parse_corner_suffix(mode, "corner_region_")
    if corner_reg !== nothing
        low_x, low_y = corner_reg
        xlo, xhi = edge_span_interval(Lx, corner_frac, low_x)
        ylo, yhi = edge_span_interval(Ly, corner_frac, low_y)
        xs = stride_axis(xlo, xhi, rx)
        ys = stride_axis(ylo, yhi, ry)
        return [(x, y) for x in xs for y in ys]
    end

    error("Unknown specloc_mode: $mode")
end

# Greedy balancing: allocate chunk counts across dimensions until product meets/exceeds target.
function choose_chunk_counts(lengths::Dict{Symbol, Int}, target_rows::Int)
    target_rows = max(1, target_rows)
    dims = collect(keys(lengths))
    nchunks = Dict(d => 1 for d in dims)

    function current_rows()
        prod(values(nchunks))
    end

    while current_rows() < target_rows
        # Expand the dimension with largest remaining reducible span.
        best_dim = nothing
        best_score = -Inf
        for d in dims
            L = lengths[d]
            c = nchunks[d]
            if c < L
                score = L / c
                if score > best_score
                    best_score = score
                    best_dim = d
                end
            end
        end
        best_dim === nothing && break
        nchunks[best_dim] += 1
    end

    return nchunks
end

function effective_chunk_count(len::Int, requested_chunks::Int)
    requested_chunks = clamp(requested_chunks, 1, max(1, len))
    vals_per_chunk = ceil(Int, len / requested_chunks)
    return ceil(Int, len / vals_per_chunk)
end

function emitted_rows_for_chunk_target(
    lengths::Dict{Symbol, Int},
    chunk_target::Int,
    rows_per_chunk_product::Int,
)
    n_chunks = choose_chunk_counts(lengths, max(1, chunk_target))
    effective_product = prod(
        effective_chunk_count(lengths[d], n_chunks[d])
        for d in keys(lengths)
    )
    return effective_product * rows_per_chunk_product
end

function nearest_allocator_rows(
    lengths::Dict{Symbol, Int},
    rows_per_chunk_product::Int,
    target_rows::Int,
)
    target_rows = max(1, target_rows)
    rows_per_chunk_product = max(1, rows_per_chunk_product)
    base_chunk_target = ceil(Int, target_rows / rows_per_chunk_product)
    max_chunk_target = prod(values(lengths))

    lower = nothing
    upper = nothing
    seen_rows = Set{Int}()

    delta = 0
    while true
        candidates = delta == 0 ? (base_chunk_target,) : (base_chunk_target - delta, base_chunk_target + delta)
        for cand in candidates
            if cand < 1 || cand > max_chunk_target
                continue
            end
            rows = emitted_rows_for_chunk_target(lengths, cand, rows_per_chunk_product)
            if rows in seen_rows
                continue
            end
            push!(seen_rows, rows)
            if rows <= target_rows && (lower === nothing || rows > lower)
                lower = rows
            end
            if rows >= target_rows && (upper === nothing || rows < upper)
                upper = rows
            end
        end

        if lower !== nothing && upper !== nothing
            break
        end
        if base_chunk_target - delta <= 1 && base_chunk_target + delta >= max_chunk_target
            break
        end
        delta += 1
    end

    return lower, upper
end

# =====================================================================
# Lightweight analytic band-edge estimator (for dynamic E range)
# =====================================================================

function scalar_shift(kx, ky, gamma, perturbation_type::Symbol)
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

function estimate_band_extrema(
    A::Float64,
    B::Float64,
    m::Float64,
    gamma::Float64,
    B_y::Float64,
    perturbation_type::Symbol;
    nk::Int=81,
)
    kx_vals = collect(range(-pi, pi; length=nk + 1))[1:end-1]
    ky_vals = collect(range(-pi, pi; length=nk + 1))[1:end-1]

    e_min = Inf
    e_max = -Inf

    @inbounds for kx in kx_vals, ky in ky_vals
        dx = A * sin(kx)
        dy = A * sin(ky)
        dz = m + 2B - B * cos(kx) - B * B_y * cos(ky)
        scalar = scalar_shift(kx, ky, gamma, perturbation_type)
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

function dynamic_energy_range_for_row(
    As_row,
    Bs_row,
    ms_row,
    gammas_row,
    B_y,
    perturbation_type::Symbol;
    margin_frac::Float64=0.10,
    nE::Int=101,
    nk::Int=81,
)
    emin = Inf
    emax = -Inf

    for A in As_row, B in Bs_row, m in ms_row, gamma in gammas_row
        lo, hi = estimate_band_extrema(Float64(A), Float64(B), Float64(m), Float64(gamma), Float64(B_y), perturbation_type; nk=nk)
        if lo < emin
            emin = lo
        end
        if hi > emax
            emax = hi
        end
    end

    width = emax - emin
    if !(width > 0)
        width = max(abs(emax), 1.0)
    end

    pad = margin_frac * width
    return collect(range(emin - pad, emax + pad; length=nE)), (emin, emax)
end

# =====================================================================
# Build chunking plan
# =====================================================================

chunk_lengths = Dict(
    :A => length(As),
    :B => length(Bs),
    :m => length(ms),
    :gamma => length(gamma_vals),
    :W => length(W_vals),
    :kappa => length(kappa_vals),
)

vals_per_row = Dict{Symbol, Int}()

# Rows are emitted per (chunk tuple) multiplied by fixed outer combinations.
# Compute the fixed multiplier so target_rows can map to a chunk-product target.
specloc_points_per_outer_combo = sum(
    begin
        pts = build_specloc_points(
            specloc_mode,
            Lx_obc,
            Ly_obc;
            base_x=specloc_x,
            base_y=specloc_y,
            res_x=specloc_resolution_x,
            res_y=specloc_resolution_y,
            outside_frac=specloc_outside_fraction,
            centre_frac=specloc_centre_fraction,
            mid_edge_frac=specloc_mid_edge_fraction,
            corner_frac=specloc_corner_fraction,
        )
        pts = uniquify_points(pts)
        isempty(pts) && error("specloc_mode=$(specloc_mode) generated zero points for Lx=$(Lx_obc), Ly=$(Ly_obc)")
        length(pts)
    end
    for Lx_obc in Lx_obcs, Ly_obc in Ly_obcs
)

fixed_outer_combos =
    length(B_ys) *
    length(perturbation_types) *
    length(disorder_types) *
    length(boundary_cond_twists) *
    length(Lx_ribbons)

rows_per_chunk_product = fixed_outer_combos * specloc_points_per_outer_combo

if allocation_mode == :explicit_vals_per_row
    vals_per_row[:A] = max(1, A_vals_per_row)
    vals_per_row[:B] = max(1, B_vals_per_row)
    vals_per_row[:m] = max(1, m_vals_per_row)
    vals_per_row[:gamma] = max(1, gamma_vals_per_row)
    vals_per_row[:W] = max(1, W_vals_per_row)
    vals_per_row[:kappa] = max(1, kappa_vals_per_row)
elseif allocation_mode == :target_rows
    chunk_product_target = ceil(Int, target_number_of_rows / max(1, rows_per_chunk_product))
    n_chunks = choose_chunk_counts(chunk_lengths, max(1, chunk_product_target))
    vals_per_row[:A] = ceil(Int, length(As) / n_chunks[:A])
    vals_per_row[:B] = ceil(Int, length(Bs) / n_chunks[:B])
    vals_per_row[:m] = ceil(Int, length(ms) / n_chunks[:m])
    vals_per_row[:gamma] = ceil(Int, length(gamma_vals) / n_chunks[:gamma])
    vals_per_row[:W] = ceil(Int, length(W_vals) / n_chunks[:W])
    vals_per_row[:kappa] = ceil(Int, length(kappa_vals) / n_chunks[:kappa])
else
    error("Unknown allocation_mode: $allocation_mode")
end

A_chunks = chunk_vector(As, vals_per_row[:A]; mode=split_mode)
B_chunks = chunk_vector(Bs, vals_per_row[:B]; mode=split_mode)
m_chunks = chunk_vector(ms, vals_per_row[:m]; mode=split_mode)
gamma_chunks = chunk_vector(gamma_vals, vals_per_row[:gamma]; mode=split_mode)
W_chunks = chunk_vector(W_vals, vals_per_row[:W]; mode=split_mode)
kappa_chunks = chunk_vector(kappa_vals, vals_per_row[:kappa]; mode=split_mode)

# =====================================================================
# Build rows
# =====================================================================

rows = NamedTuple[]
row_energy_bounds = Tuple{Float64, Float64}[]
specloc_counts_per_geometry = Int[]

nearest_target_rows_lower = nothing
nearest_target_rows_upper = nothing
if allocation_mode == :target_rows
    nearest_target_rows_lower, nearest_target_rows_upper = nearest_allocator_rows(
        chunk_lengths,
        rows_per_chunk_product,
        target_number_of_rows,
    )
end

for A_chunk in A_chunks,
    B_chunk in B_chunks,
    m_chunk in m_chunks,
    gamma_chunk in gamma_chunks,
    W_chunk in W_chunks,
    kappa_chunk in kappa_chunks,
    B_y in B_ys,
    perturbation_type in perturbation_types,
    disorder_type in disorder_types,
    winding_number in boundary_cond_twists,
    Lx_ribbon in Lx_ribbons,
    Lx_obc in Lx_obcs,
    Ly_obc in Ly_obcs

    # Store energy metadata for per-case computation in main.jl
    # Don't pre-compute E_vals; let each case compute its own range

    specloc_points = build_specloc_points(
        specloc_mode,
        Lx_obc,
        Ly_obc;
        base_x=specloc_x,
        base_y=specloc_y,
        res_x=specloc_resolution_x,
        res_y=specloc_resolution_y,
        outside_frac=specloc_outside_fraction,
        centre_frac=specloc_centre_fraction,
        mid_edge_frac=specloc_mid_edge_fraction,
        corner_frac=specloc_corner_fraction,
    )
    specloc_points = uniquify_points(specloc_points)
    isempty(specloc_points) && error("specloc_mode=$(specloc_mode) generated zero points for Lx=$(Lx_obc), Ly=$(Ly_obc)")
    push!(specloc_counts_per_geometry, length(specloc_points))

    # For dynamic E mode, optionally split (A,B,m,gamma) so each emitted row maps
    # to one physical case and can receive its own energy range in main.jl.
    A_groups = (energy_range_mode == :dynamic_band && dynamic_energy_split_cases_active) ? [[a] for a in A_chunk] : [A_chunk]
    B_groups = (energy_range_mode == :dynamic_band && dynamic_energy_split_cases_active) ? [[b] for b in B_chunk] : [B_chunk]
    m_groups = (energy_range_mode == :dynamic_band && dynamic_energy_split_cases_active) ? [[mm] for mm in m_chunk] : [m_chunk]
    gamma_groups = (energy_range_mode == :dynamic_band && dynamic_energy_split_cases_active) ? [[g] for g in gamma_chunk] : [gamma_chunk]

    for A_group in A_groups,
        B_group in B_groups,
        m_group in m_groups,
        gamma_group in gamma_groups,
        (sx, sy) in specloc_points
        push!(rows, (
            As=A_group,
            Bs=B_group,
            ms=m_group,
            B_y=B_y,
            perturbation_type=perturbation_type,
            disorder_type=disorder_type,
            winding_number=winding_number,
            Lx_ribbon=Lx_ribbon,
            Lx_obc=Lx_obc,
            Ly_obc=Ly_obc,
            gamma_vals=gamma_group,
            W_vals=W_chunk,
            kappa_vals=kappa_chunk,
            energy_range_mode=energy_range_mode,
            energy_points=energy_points,
            energy_margin_fraction=energy_margin_fraction,
            energy_scan_nk=energy_scan_nk,
            fixed_E_vals=fixed_E_vals,
            N_ky=N_ky,
            Nkx_bulk=Nkx_bulk,
            Nky_bulk=Nky_bulk,
            Nkx_3d=Nkx_3d,
            Nky_3d=Nky_3d,
            ldos_target_E=ldos_target_E,
            ldos_eta=ldos_eta,
            dos_eta=dos_eta,
            lowest_energy_count=lowest_energy_count,
            specloc_x=sx,
            specloc_y=sy,
            seed=seed,
        ))
        push!(row_energy_bounds, (Inf, -Inf))  # placeholder; actual bounds computed per case
    end
end

# =====================================================================
# Save .dat
# =====================================================================

function row_to_strings(row)
    return [
        repr(row.As),
        repr(row.Bs),
        repr(row.ms),
        repr(row.B_y),
        repr(row.perturbation_type),
        repr(row.disorder_type),
        repr(row.winding_number),
        repr(row.Lx_ribbon),
        repr(row.Lx_obc),
        repr(row.Ly_obc),
        repr(row.gamma_vals),
        repr(row.W_vals),
        repr(row.kappa_vals),
        repr(row.energy_range_mode),
        repr(row.energy_points),
        repr(row.energy_margin_fraction),
        repr(row.energy_scan_nk),
        repr(row.fixed_E_vals),
        repr(row.N_ky),
        repr(row.Nkx_bulk),
        repr(row.Nky_bulk),
        repr(row.Nkx_3d),
        repr(row.Nky_3d),
        repr(row.ldos_target_E),
        repr(row.ldos_eta),
        repr(row.dos_eta),
        repr(row.lowest_energy_count),
        repr(row.specloc_x),
        repr(row.specloc_y),
        repr(row.seed),
    ]
end

header = [
    "As", "Bs", "ms", "B_y", "perturbation_type", "disorder_type", "winding_number",
    "Lx_ribbon", "Lx_obc", "Ly_obc",
    "gamma_vals", "W_vals", "kappa_vals",
    "energy_range_mode", "energy_points", "energy_margin_fraction", "energy_scan_nk", "fixed_E_vals",
    "N_ky", "Nkx_bulk", "Nky_bulk", "Nkx_3d", "Nky_3d",
    "ldos_target_E", "ldos_eta", "dos_eta", "lowest_energy_count",
    "specloc_x", "specloc_y", "seed",
]

timestamp = Dates.format(now(), "yyyymmdd_HHMMSS")
alloc_tag = allocation_mode == :target_rows ? "target$(target_number_of_rows)" : "explicit"
energy_tag = energy_range_mode == :dynamic_band ? "Eauto" : "Efixed"
out_file = joinpath(outdir, "params_$(timestamp)_$(alloc_tag)_$(energy_tag)_rows$(length(rows)).dat")

open(out_file, "w") do io
    println(io, join(header, '\t'))
    for r in rows
        println(io, join(row_to_strings(r), '\t'))
    end
end

# =====================================================================
# Diagnostics summary
# =====================================================================

row_combo_counts = Int[]
for r in rows
    ncomb = length(r.As) * length(r.Bs) * length(r.ms) * length(r.gamma_vals) * length(r.W_vals) * length(r.kappa_vals) * r.energy_points
    push!(row_combo_counts, ncomb)
end

global_emin = Inf
global_emax = -Inf

println("\n================ Parameter Prep Summary ================")
println("Output file: $(out_file)")
println("Rows: $(length(rows))")
if allocation_mode == :target_rows
    println("Target rows requested: $(target_number_of_rows)")
    println("Nearest achievable rows: lower=$(nearest_target_rows_lower), upper=$(nearest_target_rows_upper)")
end
println("Allocation mode: $(allocation_mode)")
println("Split mode: $(split_mode)")
println("Energy mode: $(energy_range_mode)")
println("Energy margin fraction: $(energy_margin_fraction)")
println("Energy scan nk: $(energy_scan_nk)")
println("Dynamic energy split cases (configured): $(dynamic_energy_split_cases)")
println("Dynamic energy split cases (active): $(dynamic_energy_split_cases_active)")
if dynamic_split_disabled_for_target_rows
    println("  note: disabled in :target_rows mode so row targeting remains effective")
end
println("Specloc mode: $(specloc_mode)")
println("--------------------------------------------------------")
println("Ranges:    min - max - count")
println("  A:       $(range_info(As))")
println("  B:       $(range_info(Bs))")
println("  m:       $(range_info(ms))")
println("  gamma:   $(range_info(gamma_vals))")
println("  W:       $(range_info(W_vals))")
println("  kappa:   $(range_info(kappa_vals))")
println("  B_y:     $(range_info(B_ys))")
println("  Lx_rib:  $(range_info(Lx_ribbons))")
println("  Lx_obc:  $(range_info(Lx_obcs))")
println("  Ly_obc:  $(range_info(Ly_obcs))")
println("--------------------------------------------------------")
println("Chunk settings (vals per row):")
println("  A=$(vals_per_row[:A]), B=$(vals_per_row[:B]), m=$(vals_per_row[:m]), gamma=$(vals_per_row[:gamma]), W=$(vals_per_row[:W]), kappa=$(vals_per_row[:kappa])")
println("Chunk counts:")
println("  A=$(length(A_chunks)), B=$(length(B_chunks)), m=$(length(m_chunks)), gamma=$(length(gamma_chunks)), W=$(length(W_chunks)), kappa=$(length(kappa_chunks))")
println("--------------------------------------------------------")
println("Specloc points per geometry row:")
println("  min=$(minimum(specloc_counts_per_geometry)), mean=$(round(mean(specloc_counts_per_geometry), digits=1)), max=$(maximum(specloc_counts_per_geometry))")
println("  x-range=$(minimum(getindex.(rows, :specloc_x))) to $(maximum(getindex.(rows, :specloc_x)))")
println("  y-range=$(minimum(getindex.(rows, :specloc_y))) to $(maximum(getindex.(rows, :specloc_y)))")
println("--------------------------------------------------------")
println("Per-row combinations (A*B*m*gamma*W*kappa*E):")
println("  min=$(minimum(row_combo_counts)), mean=$(round(mean(row_combo_counts), digits=1)), max=$(maximum(row_combo_counts))")
println("Energy computation:")
if energy_range_mode == :fixed
    println("  mode: fixed E range [$(minimum(fixed_E_vals)), $(maximum(fixed_E_vals))], points=$(length(fixed_E_vals))")
    println("  all cases use same energy grid")
else
    println("  mode: dynamic_band (per-case extrema)")
    println("  each case will compute its own E range with:")
    println("    margin_frac=$(energy_margin_fraction), points=$(energy_points), scan_nk=$(energy_scan_nk)")
end
println("========================================================\n")
