using Eirene

#%%
# ==============================
#  ======== Tested code ========


#%%
# ================================
#  ======== Untested code ========

"""
    get_barcodes(results_eirene::Dict, max_dim::Integer; min_dim::Int=1)

Calls Eirene.barcode for 'dim' in range from `min_dim` up to 'max_dim' and
stack the resulting Arrays into a vector.

The returned value is a Vector of Arrays{Float64,2}. Each array is of
different size, because of different number of detected cycles. First column
of each array contains birth step, second column contains death step.

Arrays in returned vector correspond to Betti curve dimensions form range
`min_dim` up to 'max_dim'.

"""
function get_barcodes(results_eirene::Dict, max_dim::Integer; min_dim::Int=1, sorted::Bool=false)
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


function get_birth_death_ratio(barcodes; max_dim::Integer=3)
    birth_death_ratio_π = [barcodes[k][:, 2] ./ barcodes[k][:, 1] for k in 1:max_dim]
    return birth_death_ratio_π
end

function get_barcode_lifetime(barcodes; max_dim::Integer=3)
    lifetime = [barcodes[k][:, 2] .- barcodes[k][:, 1] for k in 1:max_dim]
    return lifetime
end

#%%
"""
    get_barcode_max_lifetime(lifetimes, min_dim, max_dim)

Returns the maximal life times of barcode for all dimensions.
"""
function get_barcode_max_lifetime(lifetimes)
    total_lifetimes = length(lifetimes)
    all_max_lifetimes = zeros(total_lifetimes, 1)

    for k in 1:total_lifetimes
        all_max_lifetimes[k] = findmax(lifetimes[k])[1]
    end

    return all_max_lifetimes
end


#%%
"""
    get_barcode_max_db_ratios(lifetimes, min_dim, max_dim)

Returns the maximal life times of barcode for all dimensions.
"""
function get_barcode_max_db_ratios(db_ratos)
    total_db = length(db_ratos)
    all_max_db = zeros(total_db, 1)

    for k in 1:total_db
        all_max_db[k] = findmax(db_ratos[k])[1]
    end

    return all_max_db
end


"""
    get_normalised_barcodes(barcodes::Vector, betti_numbers::Array)

Returns the barcodes which values are within [0,1] range.
1. get corresponding bettis
2. get max number of steps
3. divide all values by total number of steps.

"""
function get_normalised_barcodes(barcodes, betti_numbers::Array; to_total_steps=true)
    if typeof(betti_numbers) == Vector
        total_steps = size(betti_numbers[1], 1)
    else
        total_steps = size(betti_numbers, 1)
    end

    if to_total_steps
        final_total_steps = max(total_steps, betti_numbers[end, 1])
    else
        final_total_steps = max(total_steps, betti_numbers[end-1, 1])
    end

    return barcodes ./ final_total_steps
end


"""
    get_normalised_barcodes_collection(barcodes_collection, bettis_collection)

Applies get_normalised_barcodes to the collection of barcodes and corresponding
betti curves.
"""
function get_normalised_barcodes_collection(barcodes_collection, bettis_collection)
    if size(barcodes_collection, 1) != size(bettis_collection, 1)
        throw(BoundsError(barcodes_collection, bettis_collection,
            "Both collections must have same number of elements",
        ))
    else
        total_collections = size(barcodes_collection, 1)
    end
    normed_collection = deepcopy(barcodes_collection)

    for k = 1:total_collections
        normed_collection[k] = get_normalised_barcodes(barcodes_collection[k], bettis_collection[k])
    end
    return normed_collection
end
