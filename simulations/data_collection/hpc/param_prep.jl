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

As = [1.0] #collect(0.5:0.5:3.0)
Bs = [1.0] #collect(0.5:0.5:3.0)
ms = [-1.0] #collect(-5.0:0.5:1.0)

# Available perturbation types:
#   :none, :sym_cos_sum, :sym_cos_diff, :sym_cos_add, :sym_cos_sub,
#   :asym_sin_sum, :asym_sin_diff, :asym_sin_add, :asym_sin_sub, :tilt
perturbation_types = [:sym_cos_sum]

# Available disorder types: :none, :anderson, :mass
disorder_types = [:none]

Lx_obcs = [20]
Ly_obcs = [20]

gamma_vals = collect(-3.0:0.5:3.0)
W_vals     = [0.0]
kappa_vals = [2e-1]          # used when scale_kappa_to_L = false

n_disorder_realisations = 1  # >1 only meaningful when disorder_type != :none and W > 0

scale_kappa_to_L = false
kappa_scales     = [0.0004]  # effective κ = scale × Lx_obc when scale_kappa_to_L = true

orbital_displacements = [0.0]
phis                  = [0.0]

energy_range_mode       = :dynamic_band   # :fixed or :dynamic_band
fixed_E_vals            = collect(range(-5.0, 5.0; length=101))
energy_points           = 51
energy_margin_fraction  = 0.10
energy_scan_nk          = 11
dynamic_energy_split_cases = true

# Spectral localiser placement mode
# Modes: :manual, :all, :all_plus_outside, :centre_point, :centre_region,
#        :mid_edge_point_x/y, :mid_edge_region_x/y,
#        :corner_point_00/01/10/11, :corner_region_00/01/10/11
specloc_mode              = :centre_point
specloc_x                 = 7
specloc_y                 = 7
specloc_resolution_x      = 1
specloc_resolution_y      = 1
specloc_outside_fraction  = 0.10
specloc_centre_fraction   = 0.30
specloc_mid_edge_fraction = 0.30
specloc_corner_fraction   = 0.20

seed = 1234

allocation_mode       = :target_rows   # :explicit_vals_per_row or :target_rows
target_number_of_rows = 200

dynamic_energy_split_cases_active      = dynamic_energy_split_cases
dynamic_split_disabled_for_target_rows = false
if allocation_mode == :target_rows && dynamic_energy_split_cases_active
    dynamic_energy_split_cases_active      = false
    dynamic_split_disabled_for_target_rows = true
end

# Used when allocation_mode == :explicit_vals_per_row
A_vals_per_row     = max(1, length(As))
B_vals_per_row     = max(1, length(Bs))
m_vals_per_row     = max(1, length(ms))
gamma_vals_per_row = max(1, length(gamma_vals))
W_vals_per_row     = max(1, length(W_vals))
kappa_axis_per_row = max(1, length(scale_kappa_to_L ? kappa_scales : kappa_vals))
d_vals_per_row     = max(1, length(orbital_displacements))
phi_vals_per_row   = max(1, length(phis))

split_mode = :contiguous   # :contiguous or :roundrobin

# =====================================================================
# Helpers
# =====================================================================

function range_info(vec)
    isempty(vec) && return "none"
    eltype(vec) <: Integer && return "$(vec[1])-$(vec[end])-$(length(vec))"
    return "$(vec[1])-$(vec[end])-$(length(vec))"
end

function chunk_vector(v::Vector, n::Int; mode::Symbol=:contiguous)
    isempty(v) && return [v]
    n = max(1, n)
    if mode == :contiguous
        chunks = Vector{typeof(v)}()
        i = 1
        while i <= length(v)
            push!(chunks, v[i:min(length(v), i+n-1)]); i += n
        end
        return chunks
    elseif mode == :roundrobin
        nc = ceil(Int, length(v)/n)
        ch = [typeof(v)() for _ in 1:nc]
        for (i, x) in enumerate(v); push!(ch[((i-1)%nc)+1], x); end
        return ch
    else
        error("Unknown split mode: $mode")
    end
end

function stride_axis(lo::Int, hi::Int, step::Int)
    step = max(1, step);  lo > hi && return Int[]
    vals = collect(lo:step:hi)
    (isempty(vals) || vals[end] != hi) && push!(vals, hi)
    return vals
end

bounded_fraction(x) = clamp(Float64(x), 0.0, 1.0)
center_index(L)     = cld(L, 2)

function center_interval(L, frac)
    c = center_index(L);  half = floor(Int, 0.5*bounded_fraction(frac)*L)
    return max(1, c-half), min(L, c+half)
end

function edge_span_interval(L, frac, side_low::Bool)
    span = max(1, ceil(Int, bounded_fraction(frac)*L))
    side_low ? (1, min(L, span)) : (max(1, L-span+1), L)
end

function parse_corner_suffix(mode::Symbol, prefix::String)
    s = String(mode);  startswith(s, prefix) || return nothing
    sfx = s[length(prefix)+1:end];  length(sfx) == 2 || return nothing
    xb, yb = sfx[1], sfx[2]
    (xb in ('0','1') && yb in ('0','1')) || return nothing
    return (xb=='0', yb=='0')
end

uniquify_points(pts) = sort(unique(pts), by=p->(p[1],p[2]))

function build_specloc_points(
    mode::Symbol, Lx::Int, Ly::Int;
    base_x, base_y, res_x, res_y,
    outside_frac, centre_frac, mid_edge_frac, corner_frac,
)
    rx = max(1, res_x);  ry = max(1, res_y)
    if     mode == :manual;        return [(base_x, base_y)]
    elseif mode == :all;           return [(x,y) for x in stride_axis(1,Lx,rx) for y in stride_axis(1,Ly,ry)]
    elseif mode == :all_plus_outside
        ox = ceil(Int, bounded_fraction(outside_frac)*Lx)
        oy = ceil(Int, bounded_fraction(outside_frac)*Ly)
        return [(x,y) for x in stride_axis(1-ox,Lx+ox,rx) for y in stride_axis(1-oy,Ly+oy,ry)]
    elseif mode == :centre_point;  return [(center_index(Lx), center_index(Ly))]
    elseif mode == :centre_region
        xlo,xhi = center_interval(Lx, centre_frac);  ylo,yhi = center_interval(Ly, centre_frac)
        return [(x,y) for x in stride_axis(xlo,xhi,rx) for y in stride_axis(ylo,yhi,ry)]
    elseif mode == :mid_edge_point_x; return [(1,center_index(Ly)),(Lx,center_index(Ly))]
    elseif mode == :mid_edge_point_y; return [(center_index(Lx),1),(center_index(Lx),Ly)]
    elseif mode == :mid_edge_region_x
        ylo,yhi = center_interval(Ly,mid_edge_frac);  ys = stride_axis(ylo,yhi,ry)
        return vcat([(1,y) for y in ys], [(Lx,y) for y in ys])
    elseif mode == :mid_edge_region_y
        xlo,xhi = center_interval(Lx,mid_edge_frac);  xs = stride_axis(xlo,xhi,rx)
        return vcat([(x,1) for x in xs], [(x,Ly) for x in xs])
    end
    cp = parse_corner_suffix(mode, "corner_point_")
    if cp !== nothing; lx,ly = cp; return [(lx ? 1 : Lx, ly ? 1 : Ly)]; end
    cr = parse_corner_suffix(mode, "corner_region_")
    if cr !== nothing
        lx,ly = cr
        xlo,xhi = edge_span_interval(Lx,corner_frac,lx);  ylo,yhi = edge_span_interval(Ly,corner_frac,ly)
        return [(x,y) for x in stride_axis(xlo,xhi,rx) for y in stride_axis(ylo,yhi,ry)]
    end
    error("Unknown specloc_mode: $mode")
end

function choose_chunk_counts(lengths::Dict{Symbol,Int}, target::Int)
    target = max(1, target);  dims = collect(keys(lengths))
    nc = Dict(d=>1 for d in dims)
    while prod(values(nc)) < target
        best, bscore = nothing, -Inf
        for d in dims
            nc[d] < lengths[d] || continue
            s = lengths[d] / nc[d];  s > bscore && (bscore=s; best=d)
        end
        best === nothing && break;  nc[best] += 1
    end
    return nc
end

function emitted_rows(lengths, chunk_target, rpc)
    nc = choose_chunk_counts(lengths, max(1, chunk_target))
    prod(begin
        c = clamp(nc[d], 1, max(1, lengths[d]))
        ceil(Int, lengths[d] / ceil(Int, lengths[d]/c))
    end for d in keys(lengths)) * rpc
end

function nearest_neighbor_rows(lengths, rpc, actual)
    seen = Dict{Int,Int}()
    for ct in 1:prod(values(lengths))
        r = emitted_rows(lengths, ct, rpc);  haskey(seen,r) || (seen[r]=ct)
    end
    sorted = sort(collect(keys(seen)))
    lo = up = nothing
    for (i,r) in enumerate(sorted)
        r == actual || continue
        i>1              && (lo = sorted[i-1])
        i<length(sorted) && (up = sorted[i+1])
        break
    end
    return lo, up
end

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
        dz  = m + 2B - B*cos(kx) - B*cos(ky)
        dn  = sqrt((A*sin(kx))^2 + (A*sin(ky))^2 + dz^2)
        sc  = scalar_shift_energy(kx, ky, gamma, pt)
        e1 = sc-dn; e2 = sc+dn
        e1 < emin && (emin=e1);  e2 > emax && (emax=e2)
    end
    return emin, emax
end

# =====================================================================
# Active kappa axis and chunk lengths
# =====================================================================

kappa_axis = scale_kappa_to_L ? kappa_scales : kappa_vals

chunk_lengths = Dict(
    :A     => length(As),
    :B     => length(Bs),
    :m     => length(ms),
    :gamma => length(gamma_vals),
    :W     => length(W_vals),
    :kappa => length(kappa_axis),
    :d     => length(orbital_displacements),
    :phi   => length(phis),
)

vals_per_row = Dict{Symbol,Int}()

specloc_pts_per_geo = sum(
    length(uniquify_points(build_specloc_points(
        specloc_mode, Lx_obc, Ly_obc;
        base_x=specloc_x, base_y=specloc_y,
        res_x=specloc_resolution_x, res_y=specloc_resolution_y,
        outside_frac=specloc_outside_fraction, centre_frac=specloc_centre_fraction,
        mid_edge_frac=specloc_mid_edge_fraction, corner_frac=specloc_corner_fraction)))
    for Lx_obc in Lx_obcs, Ly_obc in Ly_obcs
)

fixed_outer_combos     = length(perturbation_types) * length(disorder_types)
rows_per_chunk_product = fixed_outer_combos * specloc_pts_per_geo

if allocation_mode == :explicit_vals_per_row
    vals_per_row[:A]     = A_vals_per_row
    vals_per_row[:B]     = B_vals_per_row
    vals_per_row[:m]     = m_vals_per_row
    vals_per_row[:gamma] = gamma_vals_per_row
    vals_per_row[:W]     = W_vals_per_row
    vals_per_row[:kappa] = kappa_axis_per_row
    vals_per_row[:d]     = d_vals_per_row
    vals_per_row[:phi]   = phi_vals_per_row
elseif allocation_mode == :target_rows
    n_chunks = choose_chunk_counts(chunk_lengths,
                   ceil(Int, target_number_of_rows / max(1, rows_per_chunk_product)))
    for k in keys(chunk_lengths)
        vals_per_row[k] = ceil(Int, chunk_lengths[k] / n_chunks[k])
    end
else
    error("Unknown allocation_mode: $allocation_mode")
end

A_chunks     = chunk_vector(As,                  vals_per_row[:A];     mode=split_mode)
B_chunks     = chunk_vector(Bs,                  vals_per_row[:B];     mode=split_mode)
m_chunks     = chunk_vector(ms,                  vals_per_row[:m];     mode=split_mode)
gamma_chunks = chunk_vector(gamma_vals,          vals_per_row[:gamma]; mode=split_mode)
W_chunks     = chunk_vector(W_vals,              vals_per_row[:W];     mode=split_mode)
kappa_chunks = chunk_vector(kappa_axis,          vals_per_row[:kappa]; mode=split_mode)
d_chunks     = chunk_vector(orbital_displacements, vals_per_row[:d];   mode=split_mode)
phi_chunks   = chunk_vector(phis,                vals_per_row[:phi];   mode=split_mode)

# =====================================================================
# Build rows
# =====================================================================

rows = NamedTuple[]
specloc_counts = Int[]

for A_chunk  in A_chunks,
    B_chunk  in B_chunks,
    m_chunk  in m_chunks,
    gamma_chunk in gamma_chunks,
    W_chunk  in W_chunks,
    kappa_chunk in kappa_chunks,
    d_chunk  in d_chunks,
    phi_chunk in phi_chunks,
    perturbation_type in perturbation_types,
    disorder_type     in disorder_types,
    Lx_obc in Lx_obcs,
    Ly_obc in Ly_obcs

    pts = uniquify_points(build_specloc_points(
        specloc_mode, Lx_obc, Ly_obc;
        base_x=specloc_x, base_y=specloc_y,
        res_x=specloc_resolution_x, res_y=specloc_resolution_y,
        outside_frac=specloc_outside_fraction, centre_frac=specloc_centre_fraction,
        mid_edge_frac=specloc_mid_edge_fraction, corner_frac=specloc_corner_fraction))
    isempty(pts) && error("specloc_mode=$specloc_mode generated zero points")
    push!(specloc_counts, length(pts))

    kappa_vals_row   = scale_kappa_to_L ? kappa_vals  : kappa_chunk
    kappa_scales_row = scale_kappa_to_L ? kappa_chunk : kappa_scales

    A_groups     = (energy_range_mode==:dynamic_band && dynamic_energy_split_cases_active) ? [[a] for a in A_chunk] : [A_chunk]
    B_groups     = (energy_range_mode==:dynamic_band && dynamic_energy_split_cases_active) ? [[b] for b in B_chunk] : [B_chunk]
    m_groups     = (energy_range_mode==:dynamic_band && dynamic_energy_split_cases_active) ? [[mm] for mm in m_chunk] : [m_chunk]
    gamma_groups = (energy_range_mode==:dynamic_band && dynamic_energy_split_cases_active) ? [[g] for g in gamma_chunk] : [gamma_chunk]

    for A_group in A_groups, B_group in B_groups, m_group in m_groups,
        gamma_group in gamma_groups, d in d_chunk, phi in phi_chunk,
        (sx, sy) in pts

        push!(rows, (
            As=A_group, Bs=B_group, ms=m_group,
            perturbation_type=perturbation_type, disorder_type=disorder_type,
            Lx_obc=Lx_obc, Ly_obc=Ly_obc,
            gamma_vals=gamma_group, W_vals=W_chunk, kappa_vals=kappa_vals_row,
            energy_range_mode=energy_range_mode, energy_points=energy_points,
            energy_margin_fraction=energy_margin_fraction, energy_scan_nk=energy_scan_nk,
            fixed_E_vals=fixed_E_vals,
            specloc_x=sx, specloc_y=sy, seed=seed,
            orbital_displacement=d, phi=phi,
            n_disorder_realisations=n_disorder_realisations,
            scale_kappa_to_L=scale_kappa_to_L, kappa_scales=kappa_scales_row,
        ))
    end
end

nearest_lower = nearest_upper = nothing
if allocation_mode == :target_rows
    nearest_lower, nearest_upper =
        nearest_neighbor_rows(chunk_lengths, rows_per_chunk_product, length(rows))
end

# =====================================================================
# Save .dat  (23 columns)
# =====================================================================

header = [
    "As","Bs","ms","perturbation_type","disorder_type",
    "Lx_obc","Ly_obc",
    "gamma_vals","W_vals","kappa_vals",
    "energy_range_mode","energy_points","energy_margin_fraction","energy_scan_nk","fixed_E_vals",
    "specloc_x","specloc_y","seed",
    "orbital_displacement","phi","n_disorder_realisations","scale_kappa_to_L","kappa_scales",
]

function row_to_strings(r)
    [repr(r.As), repr(r.Bs), repr(r.ms),
     repr(r.perturbation_type), repr(r.disorder_type),
     repr(r.Lx_obc), repr(r.Ly_obc),
     repr(r.gamma_vals), repr(r.W_vals), repr(r.kappa_vals),
     repr(r.energy_range_mode), repr(r.energy_points),
     repr(r.energy_margin_fraction), repr(r.energy_scan_nk), repr(r.fixed_E_vals),
     repr(r.specloc_x), repr(r.specloc_y), repr(r.seed),
     repr(r.orbital_displacement), repr(r.phi),
     repr(r.n_disorder_realisations), repr(r.scale_kappa_to_L), repr(r.kappa_scales)]
end

ts      = Dates.format(now(), "yyyymmdd_HHMMSS")
atag    = allocation_mode == :target_rows ? "target$(target_number_of_rows)" : "explicit"
etag    = energy_range_mode == :dynamic_band ? "Eauto" : "Efixed"
outfile = joinpath(outdir, "params_$(ts)_$(atag)_$(etag)_rows$(length(rows)).dat")

open(outfile, "w") do io
    println(io, join(header, '\t'))
    for r in rows; println(io, join(row_to_strings(r), '\t')); end
end

# =====================================================================
# Summary
# =====================================================================

println("\n================ Parameter Prep Summary ================")
println("Output:          $outfile")
println("Rows:            $(length(rows))")
allocation_mode == :target_rows &&
    println("Target / nbrs:   $(target_number_of_rows)  lower=$(nearest_lower)  upper=$(nearest_upper)")
println("Allocation:      $allocation_mode  |  Split: $split_mode  |  Energy: $energy_range_mode")
println("--------------------------------------------------------")
println("Parameter ranges (min–max–count):")
for (lbl, v) in [("A",As),("B",Bs),("m",ms),("gamma",gamma_vals),("W",W_vals),
                  ("kappa_axis (scale=$(scale_kappa_to_L))",kappa_axis),
                  ("d",orbital_displacements),("phi",phis),
                  ("Lx_obc",Lx_obcs),("Ly_obc",Ly_obcs)]
    println("  $(lpad(lbl,32)): $(range_info(v))")
end
println("Perturbation types:  $perturbation_types")
println("Disorder types:      $disorder_types  |  n_realisations=$n_disorder_realisations")
println("scale_kappa_to_L=$scale_kappa_to_L  kappa_scales=$kappa_scales")
println("orbital_displacements=$orbital_displacements  phis=$phis")
println("--------------------------------------------------------")
println("Chunk counts:")
for k in [:A,:B,:m,:gamma,:W,:kappa,:d,:phi]
    println("  $(k): $(length(eval(Symbol(string(k,"_chunks")))))  ($(vals_per_row[k]) vals/row)")
end
println("Specloc pts/geo:  min=$(minimum(specloc_counts))  max=$(maximum(specloc_counts))")
println("  x-range: $(minimum(getindex.(rows,:specloc_x)))–$(maximum(getindex.(rows,:specloc_x)))")
println("  y-range: $(minimum(getindex.(rows,:specloc_y)))–$(maximum(getindex.(rows,:specloc_y)))")
if energy_range_mode == :dynamic_band
    println("Energy: dynamic_band  margin=$(energy_margin_fraction)  pts=$(energy_points)  scan_nk=$(energy_scan_nk)  split_cases=$(dynamic_energy_split_cases_active)")
    dynamic_split_disabled_for_target_rows && println("  (split disabled for :target_rows mode)")
else
    println("Energy: fixed [$(minimum(fixed_E_vals)), $(maximum(fixed_E_vals))]  pts=$(length(fixed_E_vals))")
end
println("========================================================\n")
