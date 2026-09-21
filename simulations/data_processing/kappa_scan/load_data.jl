using JLD2
using Statistics
using Printf

"""
    load_specloc_result(filepath::String)

Load a spectral localiser JLD2 result file and extract key data.

Returns a NamedTuple with:
  - metadata: Dict of simulation parameters
  - axes: Dict with vectors (gammas, Ws, kappas, Es)
  - gap: Array of shape (ng, nW, nk, nE)
  - signature: Array of shape (ng, nW, nk, nE)
  - index: Array of shape (ng, nW, nk, nE)  [= 0.5 * signature]
"""
function load_specloc_result(filepath::String)
    data = JLD2.load(filepath)
    result = data["result"]
    
    return (
        metadata=result["metadata"],
        axes=result["axes"],
        gap=result["specloc"]["gap"],
        signature=result["specloc"]["signature"],
        index=result["specloc"]["index"],
    )
end

"""
    average_over_energy(data::AbstractArray, axis_vals::Vector)

Average an array over its last dimension (energy) and return
(averaged_data, axis_vals_name) tuple.
"""
function average_over_energy(data::AbstractArray, axis_name::String="Energy")
    averaged = mean(data; dims=ndims(data))
    squeezed = dropdims(averaged; dims=ndims(data))
    return squeezed, axis_name
end

"""
    extract_gap_vs_kappa(result::NamedTuple)

Extract gap vs kappa for the given result.
Assumes single gamma, single W (ng=1, nW=1).

Returns (kappas::Vector, gap_vs_kappa::Vector)
"""
function extract_gap_vs_kappa(result::NamedTuple)
    kappas = result.axes["kappas"]
    gap = result.gap
    
    # Expect shape (1, 1, nk, nE)
    @assert size(gap, 1) == 1 "Expected single gamma dimension, got $(size(gap, 1))"
    @assert size(gap, 2) == 1 "Expected single W dimension, got $(size(gap, 2))"
    
    # Extract and squeeze to shape (nk, nE)
    gap_data = dropdims(gap; dims=(1, 2))
    
    # Average over energy
    gap_vs_kappa = vec(mean(gap_data; dims=2))
    
    return kappas, gap_vs_kappa
end

"""
    extract_index_vs_kappa(result::NamedTuple)

Extract spectral localiser index vs kappa for the given result.
Assumes single gamma, single W (ng=1, nW=1).

Returns (kappas::Vector, index_vs_kappa::Vector)
"""
function extract_index_vs_kappa(result::NamedTuple)
    kappas = result.axes["kappas"]
    index = result.index
    
    # Expect shape (1, 1, nk, nE)
    @assert size(index, 1) == 1 "Expected single gamma dimension, got $(size(index, 1))"
    @assert size(index, 2) == 1 "Expected single W dimension, got $(size(index, 2))"
    
    # Extract and squeeze to shape (nk, nE)
    index_data = dropdims(index; dims=(1, 2))
    
    # Average over energy
    index_vs_kappa = vec(mean(index_data; dims=2))
    
    return kappas, index_vs_kappa
end

"""
    extract_signature_vs_kappa(result::NamedTuple)

Extract spectral localiser signature vs kappa for the given result.
Assumes single gamma, single W (ng=1, nW=1).

Returns (kappas::Vector, signature_vs_kappa::Vector)
"""
function extract_signature_vs_kappa(result::NamedTuple)
    kappas = result.axes["kappas"]
    signature = result.signature
    
    # Expect shape (1, 1, nk, nE)
    @assert size(signature, 1) == 1 "Expected single gamma dimension, got $(size(signature, 1))"
    @assert size(signature, 2) == 1 "Expected single W dimension, got $(size(signature, 2))"
    
    # Extract and squeeze to shape (nk, nE)
    sig_data = dropdims(signature; dims=(1, 2))
    
    # Average over energy
    signature_vs_kappa = vec(mean(sig_data; dims=2))
    
    return kappas, signature_vs_kappa
end
