#!/usr/bin/env julia
# setup_julia_packages.jl
# Run once per login to ensure required Julia packages are installed in the cluster
# Usage: julia setup_julia_packages.jl

import Pkg

# List of required packages
packages = [
    "JLD2",
    "KrylovKit",
]

println("Setting up Julia packages for SpecLoc HPC pipeline...")
println("Current depot path: $(first(DEPOT_PATH))")

for pkg in packages
    println("\nChecking/installing: $pkg")
    try
        Pkg.add(pkg)
        println("✓ $pkg installed/updated")
    catch e
        println("✗ Error installing $pkg: $e")
    end
end

println("\n" * repeat("=", 60))
println("Package setup complete!")
println(repeat("=", 60))

# Verify installations
println("\nVerifying installations:")
for pkg in packages
    try
        eval(:(using $(Symbol(pkg))))
        println("✓ $pkg loads successfully")
    catch
        println("✗ $pkg failed to load")
    end
end
