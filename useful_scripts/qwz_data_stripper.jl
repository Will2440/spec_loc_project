"""
    qwz_data_stripper.jl

Strip data-heavy eigensystems from JLD2 result files saved by the spectral localizer solver.
Removes Hamiltonian matrices and eigenvectors to reclaim disk space while preserving
all computed observables (specloc gap/signature, DOS/LDOS, ribbon spectra, band structure).

Usage:
    julia qwz_data_stripper.jl <job_folder_path> [--inplace]

    job_folder_path: Path to folder containing subfolders with .jld2 result files
    --inplace:       Replace original files instead of creating .stripped.jld2 copies (default: safe mode with copies)

Example:
    julia qwz_data_stripper.jl /path/to/job_results
    julia qwz_data_stripper.jl /path/to/job_results --inplace
"""

using JLD2
using Printf

function strip_result_file(filepath::String; inplace::Bool=false)
    """
    Load a result file, strip heavy data, and save it back (or to new file).
    Returns: (original_size_mb, stripped_size_mb, success::Bool, message::String)
    """
    try
        # Get original file size
        original_size = filesize(filepath) / (1024^2)  # Convert to MB
        
        # Load the data
        data = load(filepath)
        
        # Verify structure
        if !haskey(data, "hamiltonian_eigensystems") && !haskey(data, "localiser_eigensystems")
            return original_size, original_size, false, "File structure unrecognized (missing expected keys)"
        end
        
        # Strip Hamiltonian eigensystems (but keep gamma_w_pairs for reference)
        if haskey(data, "hamiltonian_eigensystems")
            ham_sys = data["hamiltonian_eigensystems"]
            if isa(ham_sys, Dict)
                # Keep gamma_w_pairs, delete matrices and eigendecompositions
                stripped_ham = Dict(
                    "gamma_w_pairs" => get(ham_sys, "gamma_w_pairs", nothing),
                    "matrices" => nothing,
                    "eigenvalues" => nothing,
                    "eigenvectors" => nothing,
                )
                data["hamiltonian_eigensystems"] = stripped_ham
            end
        end
        
        # Strip Localiser eigensystems (keep structure but nullify data)
        if haskey(data, "localiser_eigensystems")
            loc_sys = data["localiser_eigensystems"]
            if isa(loc_sys, Dict)
                data["localiser_eigensystems"] = Dict(
                    "full_saved" => get(loc_sys, "full_saved", false),
                    "eigenvalues" => nothing,
                    "eigenvectors" => nothing,
                )
            end
        end
        
        # Determine output path
        if inplace
            output_path = filepath
        else
            # Create .stripped.jld2 version
            output_path = replace(filepath, ".jld2" => ".stripped.jld2")
        end
        
        # Save stripped data
        save(output_path, data)
        
        # Get new file size
        stripped_size = filesize(output_path) / (1024^2)  # Convert to MB
        
        return original_size, stripped_size, true, "Success"
        
    catch e
        return 0.0, 0.0, false, "Error: $(string(e))"
    end
end

function find_jld2_files(root_path::String)::Vector{String}
    """Recursively find all .jld2 files in root_path and subdirectories."""
    jld2_files = String[]
    
    for (dirpath, dirs, files) in walkdir(root_path)
        for file in files
            if endswith(file, ".jld2") && !endswith(file, ".stripped.jld2")
                push!(jld2_files, joinpath(dirpath, file))
            end
        end
    end
    
    return sort(jld2_files)
end

function main()
    # Parse command line arguments
    if length(ARGS) < 1
        println("""
        Usage: julia qwz_data_stripper.jl <job_folder_path> [--inplace]
        
        job_folder_path: Path to folder containing subfolders with .jld2 result files
        --inplace:       Replace original files (default: create .stripped.jld2 copies)
        """)
        exit(1)
    end
    
    job_folder = ARGS[1]
    inplace = any(arg == "--inplace" for arg in ARGS[2:end])
    
    # Validate path
    if !isdir(job_folder)
        println("ERROR: Directory not found: $job_folder")
        exit(1)
    end
    
    # Find all .jld2 files
    println("Scanning for .jld2 files in: $job_folder")
    jld2_files = find_jld2_files(job_folder)
    
    if isempty(jld2_files)
        println("No .jld2 files found in the directory tree.")
        exit(0)
    end
    
    println("Found $(length(jld2_files)) .jld2 files to process")
    if !inplace
        println("Running in SAFE MODE: creating .stripped.jld2 copies (originals preserved)")
    else
        println("Running in INPLACE MODE: replacing original files")
    end
    println()
    
    # Process files
    total_original = 0.0
    total_stripped = 0.0
    successful = 0
    failed = 0
    failures = Tuple{String, String}[]
    
    for (i, filepath) in enumerate(jld2_files)
        rel_path = relpath(filepath, job_folder)
        print("[$i/$(length(jld2_files))] Processing: $rel_path ... ")
        
        orig_size, strip_size, success, msg = strip_result_file(filepath; inplace=inplace)
        
        if success
            total_original += orig_size
            total_stripped += strip_size
            savings = orig_size - strip_size
            percent = (savings / orig_size) * 100
            @printf("✓ %.1f MB → %.1f MB (saved %.1f MB, %.1f%%)\n", 
                    orig_size, strip_size, savings, percent)
            successful += 1
        else
            println("✗ $msg")
            failed += 1
            push!(failures, (rel_path, msg))
        end
    end
    
    # Summary
    println()
    println("=" ^ 70)
    println("SUMMARY")
    println("=" ^ 70)
    println("Files processed: $(successful + failed)")
    println("  Successful: $successful")
    println("  Failed: $failed")
    
    if total_original > 0
        total_savings = total_original - total_stripped
        percent = (total_savings / total_original) * 100
        @printf("Total disk space: %.1f MB → %.1f MB (saved %.1f MB, %.1f%%)\n",
                total_original, total_stripped, total_savings, percent)
    end
    
    if !isempty(failures)
        println()
        println("Failed files:")
        for (filepath, msg) in failures
            println("  • $filepath: $msg")
        end
    end
    
    if !inplace && successful > 0
        println()
        println("Original files preserved. Stripped versions saved as .stripped.jld2")
        println("Verify the stripped files before deleting originals.")
    end
    
    exit(failed > 0 ? 1 : 0)
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main()
end
