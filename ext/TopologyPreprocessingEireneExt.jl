module TopologyPreprocessingEireneExt

using TopologyPreprocessing
using Eirene

"""
    get_bettis(results_eirene::Dict, max_dim::Integer; min_dim::Int=1)

Calls Eirene.betticurve for each dimension from `min_dim` up to `max_dim` and
returns a Vector of Arrays{Float64,2}. Each array is (n,2): column 1 is the
filtration step, column 2 is the Betti number.
"""
function TopologyPreprocessing.get_bettis(results_eirene::Dict, max_dim::Integer; min_dim::Int = 1)
    bettis = Matrix{Float64}[]
    for d = min_dim:max_dim
        result = betticurve(results_eirene, dim = d)
        if isempty(result) && d > 1
            result = zeros(size(bettis[d-1]))
        end
        push!(bettis, result)
    end
    return bettis
end

"""
    get_barcodes(results_eirene::Dict, max_dim::Integer; min_dim::Int=1, sorted::Bool=false)

Calls Eirene.barcode for each dimension from `min_dim` up to `max_dim` and
returns a Vector of Arrays{Float64,2}. Each array is (n,2): column 1 is birth,
column 2 is death.
"""
function TopologyPreprocessing.get_barcodes(results_eirene::Dict, max_dim::Integer; min_dim::Int=1, sorted::Bool=false)
    barcodes = Matrix{Float64}[]
    for d = min_dim:max_dim
        result = barcode(results_eirene, dim=d)
        if isempty(result) && d > 1
            result = zeros(size(barcodes[d-1]))
        end
        if sorted
            result = result[sortperm(result[:, 1]), :]
        end
        push!(barcodes, result)
    end
    return barcodes
end

"""
    get_vectorized_bettis(results_eirene::Dict, max_dim::Integer; min_dim::Int=1)

Computes Betti curves for dimensions `min_dim:max_dim` and returns them as a
matrix where each column is one Betti curve (values only, no step index).
"""
function TopologyPreprocessing.get_vectorized_bettis(results_eirene::Dict, max_dim::Integer; min_dim::Int = 1)
    all_bettis = TopologyPreprocessing.get_bettis(results_eirene, max_dim, min_dim = min_dim)
    return TopologyPreprocessing.vectorize_bettis(all_bettis)
end

"""
    get_dataset_bettis_areas(dataset; min_dim=1, max_dim=3, return_matrix=true)

Runs Eirene on every matrix in `dataset`, computes Betti curves, and returns
the area under each Betti curve. Requires Eirene to be loaded.
"""
function get_dataset_bettis_areas(dataset; min_dim::Integer=1, max_dim::Integer=3, return_matrix::Bool=true)
    areas_vector = Array[]
    for data = dataset
        @info "Computing topology."
        C = eirene(data, maxdim=max_dim)
        matrix_bettis = TopologyPreprocessing.get_bettis(C, max_dim, min_dim=min_dim)
        push!(areas_vector, TopologyPreprocessing.get_area_under_betti_curve(matrix_bettis))
    end
    if return_matrix
        return vcat([areas_vector[k] for k=1:length(areas_vector)]...)
    else
        return areas_vector
    end
end

end # module TopologyPreprocessingEireneExt
