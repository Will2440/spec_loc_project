"""
    qwz_data_stripper.jl

Strip data-heavy eigensystems from JLD2 result files saved by the spectral localizer solver.
Removes Hamiltonian matrices and eigenvectors to reclaim disk space while preserving
all computed observables (specloc gap/signature, DOS/LDOS, ribbon spectra, band structure).

Creates a parallel folder structure with .stripped appended to the folder name, allowing
stripped and original datasets to be processed and managed independently.

Usage:
    julia qwz_data_stripper.jl <job_folder_path>

    job_folder_path: Path to folder containing subfolders with .jld2 result files

Example:
    julia qwz_data_stripper.jl /path/to/job_results
    
    This creates: /path/to/job_results.stripped/
    with identical folder structure and .stripped.jld2 files
"""

using JLD2
using Printf

function find_present_key(d::AbstractDict, key::String)
    if haskey(d, key)
        return key
    end
    sym = Symbol(key)
    if haskey(d, sym)
        return sym
    end
    return nothing
end

function strip_payload!(payload::AbstractDict)
    changed = false

    ham_key = find_present_key(payload, "hamiltonian_eigensystems")
    if ham_key !== nothing
        ham_sys = payload[ham_key]
        if ham_sys isa AbstractDict
            gamma_pairs_key = find_present_key(ham_sys, "gamma_w_pairs")
            gamma_pairs = gamma_pairs_key === nothing ? nothing : ham_sys[gamma_pairs_key]

            payload[ham_key] = Dict(
                "gamma_w_pairs" => gamma_pairs,
                "matrices" => nothing,
                "eigenvalues" => nothing,
                "eigenvectors" => nothing,
            )
            changed = true
        end
    end

    loc_key = find_present_key(payload, "localiser_eigensystems")
    if loc_key !== nothing
        loc_sys = payload[loc_key]
        if loc_sys isa AbstractDict
            full_saved_key = find_present_key(loc_sys, "full_saved")
            full_saved = full_saved_key === nothing ? false : loc_sys[full_saved_key]

            payload[loc_key] = Dict(
                "full_saved" => full_saved,
                "eigenvalues" => nothing,
                "eigenvectors" => nothing,
            )
            changed = true
        end
    end

    return changed
end

function strip_result_file(input_filepath::String, output_filepath::String)
    """
    Load a result file, strip heavy data, and save to output path.
    Returns: (original_size_mb, stripped_size_mb, success::Bool, message::String)
    """
    try
        # Get original file size
        original_size = filesize(input_filepath) / (1024^2)  # Convert to MB
        
        # Load the data
        data = load(input_filepath)
        
        changed = false

        # Two supported layouts:
        # 1) payload keys at top level
        # 2) @save filename result   => top-level key "result" containing payload
        if data isa AbstractDict
            changed |= strip_payload!(data)

            result_key = find_present_key(data, "result")
            if result_key !== nothing && data[result_key] isa AbstractDict
                changed |= strip_payload!(data[result_key])
            end
        end

        if !changed
            return original_size, original_size, false, "File structure unrecognized (no strip-target keys found in top-level or result payload)"
        end
        
        # Ensure output directory exists
        mkpath(dirname(output_filepath))
        
        # Save stripped data
        if data isa AbstractDict
            # Preserve original top-level variable naming where possible.
            # If source file used @save filename result, keep that shape.
            result_key = find_present_key(data, "result")
            if result_key !== nothing
                save(output_filepath, String(result_key), data[result_key])
            else
                save(output_filepath, data)
            end
        else
            save(output_filepath, "result", data)
        end
        
        # Get new file size
        stripped_size = filesize(output_filepath) / (1024^2)  # Convert to MB
        
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
        Usage: julia qwz_data_stripper.jl <job_folder_path>
        
        job_folder_path: Path to folder containing subfolders with .jld2 result files
        
        Creates: <job_folder_path>.stripped/ with identical structure and .stripped.jld2 files
        """)
        exit(1)
    end
    
    job_folder = ARGS[1]
    
    # Validate path
    if !isdir(job_folder)
        println("ERROR: Directory not found: $job_folder")
        exit(1)
    end
    
    # Create output folder at same level
    output_root = job_folder * ".stripped"
    
    # Find all .jld2 files
    println("Scanning for .jld2 files in: $job_folder")
    jld2_files = find_jld2_files(job_folder)
    
    if isempty(jld2_files)
        println("No .jld2 files found in the directory tree.")
        exit(0)
    end
    
    println("Found $(length(jld2_files)) .jld2 files to process")
    println("Output folder: $output_root")
    println()
    
    # Process files
    total_original = 0.0
    total_stripped = 0.0
    successful = 0
    failed = 0
    failures = Tuple{String, String}[]
    
    for (i, input_filepath) in enumerate(jld2_files)
        # Get relative path from job_folder
        rel_path = relpath(input_filepath, job_folder)
        
        # Create output path with .stripped filename
        output_filepath = joinpath(output_root, replace(rel_path, ".jld2" => ".stripped.jld2"))
        
        print("[$i/$(length(jld2_files))] Processing: $rel_path ... ")
        
        orig_size, strip_size, success, msg = strip_result_file(input_filepath, output_filepath)
        
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
    
    if successful > 0
        println()
        println("Original folder preserved: $job_folder")
        println("Stripped folder created: $output_root")
        println("Verify the stripped version before deleting the original folder.")
    end
    
    exit(failed > 0 ? 1 : 0)
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main()
end
