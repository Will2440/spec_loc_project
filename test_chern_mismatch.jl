using Printf
include("notebook_calcs/qwz_pbc_module.jl")
using .QWZ_Model

println("Investigating: Per-Band vs Global Accumulated Chern")
println("=" ^ 70)

# Create test data with perturbation
data = compute_bulk_band_berry_data(
    Nkx=51,
    Nky=51,
    A=1.0,
    B=1.0,
    m=-1.0,
    gamma=3.0,  # User mentioned gamma=3.0
    phi=0.0,
    orbital_displacement=0.0,
    perturbation_type=:sym_cos_sum
)

println("\n1. GLOBAL ACCUMULATED CHERN (what we query for energy lookup)")
println("-" ^ 70)
println("Range: [" * @sprintf("%.4f", minimum(data.cumulative_chern)) * 
        ", " * @sprintf("%.4f", maximum(data.cumulative_chern)) * "]")

# Show some key values
println("\nSample points in global accumulation:")
n_points = length(data.cumulative_chern)
indices = [1, n_points÷4, n_points÷2, 3*n_points÷4, n_points]
for idx in indices
    if idx <= n_points
        @printf("  Index %5d: E = %8.4f, C_global = %8.4f\n", 
                idx, data.accumulation_energies[idx], data.cumulative_chern[idx])
    end
end

println("\n2. PER-BAND ACCUMULATED CHERN (what contours show on plot)")
println("-" ^ 70)

for band in 1:2
    cum_band = data.cum_chern_per_band[band, :, :]
    min_c = minimum(cum_band)
    max_c = maximum(cum_band)
    @printf("Band %d: Range [%.4f, %.4f]\n", band, min_c, max_c)
end

println("\n3. COMPARISON: Critical Chern Values")
println("-" ^ 70)

test_cherns = [0.0, 0.5, 1.0, -0.5]
for c_target in test_cherns
    E_found = get_energy_at_chern(data.accumulation_energies, data.cumulative_chern, c_target)
    
    in_global_range = (c_target >= minimum(data.cumulative_chern)) && 
                      (c_target <= maximum(data.cumulative_chern))
    
    print("C_target = " * @sprintf("%7.4f", c_target) * " : ")
    if isnan(E_found)
        print("OUTSIDE RANGE (N/A) - ")
    else
        print("Found at E = " * @sprintf("%8.4f", E_found) * " - ")
    end
    
    in_band_1 = (c_target >= 0) && (c_target <= 1)
    in_band_2 = (c_target >= -1) && (c_target <= 0)
    
    if in_band_1 || in_band_2
        band_info = in_band_1 ? "in Band 1 per-band range" : "in Band 2 per-band range"
        println(band_info)
    else
        println("NOT in any per-band range")
    end
end

println("\n" ^ 1)
println("EXPLANATION:")
println("=" ^ 70)
println("""
The N/A "outside range" happens because:

1. The CONTOURS on the C(E, kx, ky) plots show PER-BAND accumulated Chern
   - Band 1 accumulates from 0 to its total Chern (≈+1)
   - Band 2 accumulates from 0 to its total Chern (≈-1)
   - Each band is ordered independently by energy within that band

2. The ENERGY LOOKUP uses GLOBAL accumulated Chern
   - All k-points and bands are flattened and ordered by energy globally
   - Chern values accumulate as you sweep through all states by energy
   - The range depends on the band structure overlap in energy

3. These are DIFFERENT things!
   - A Chern value visible as a contour in per-band space might not exist
     in the global energy-ordered accumulation
   - For example: -0.5 is in Band 2's per-band range, but might not appear
     when states are globally ordered by energy (if bands don't overlap 
     energetically in the right way)
""")

println("\n4. VERIFICATION at gamma=3.0")
println("-" ^ 70)
@printf("At gamma=3.0:\n")
@printf("  Global Chern range: [%.4f, %.4f]\n", 
        minimum(data.cumulative_chern), maximum(data.cumulative_chern))
println("  This determines which contour values can have energies.")
