using Printf
include("notebook_calcs/qwz_pbc_module.jl")
using .QWZ_Model

println("Testing Energy of Critical Chern Computation...")
println("=" ^ 60)

# Create test data
data = compute_bulk_band_berry_data(
    Nkx=51,
    Nky=51,
    A=1.0,
    B=1.0,
    m=-1.0,
    gamma=0.0,
    phi=0.0,
    orbital_displacement=0.0,
    perturbation_type=:none
)

println("\nAccumulated Chern Range:")
println("  Min: " * @sprintf("%.4f", minimum(data.cumulative_chern)))
println("  Max: " * @sprintf("%.4f", maximum(data.cumulative_chern)))
println("  Energy range: [" * @sprintf("%.4f", minimum(data.accumulation_energies)) * ", " * @sprintf("%.4f", maximum(data.accumulation_energies)) * "]")

println("\n" ^ 1)
println("Testing get_energy_at_chern function:")
println("-" ^ 60)

# Test cases
test_cases = [
    (0.0, "Lower band Chern (at gap)"),
    (0.5, "Halfway to lower band Chern"),
    (1.0, "Full lower band Chern"),
    (0.25, "Quarter of lower band Chern"),
    (-1.0, "Outside range (too low)"),
    (2.0, "Outside range (too high)")
]

for (chern_val, description) in test_cases
    E = get_energy_at_chern(data.accumulation_energies, data.cumulative_chern, chern_val)
    if isnan(E)
        @printf("C = %7.4f : E = N/A (outside range) - %s\n", chern_val, description)
    else
        @printf("C = %7.4f : E = %8.4f - %s\n", chern_val, E, description)
    end
end

println("\n" ^ 1)
println("Verification: Re-compute accumulated Chern at found energies:")
println("-" ^ 60)

for (chern_val, description) in [(0.5, "50% Chern"), (1.0, "Full lower band")]
    E = get_energy_at_chern(data.accumulation_energies, data.cumulative_chern, chern_val)
    if !isnan(E)
        # Find the index closest to this energy
        idx = searchsortedlast(data.accumulation_energies, E)
        if idx > 0 && idx < length(data.cumulative_chern)
            C_at_E = data.cumulative_chern[idx]
            @printf("C_target = %.4f → E = %.4f → C_computed = %.4f (diff = %.6f)\n", 
                    chern_val, E, C_at_E, abs(chern_val - C_at_E))
        end
    end
end

println("\n" ^ 60)
println("Energy of Critical Chern computation test completed! ✓")
