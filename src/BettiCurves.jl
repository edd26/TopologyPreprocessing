# ==============================
#  ======== Tested code ========
using Eirene

#%%
"""
    get_bettis(results_eirene::Dict, max_dim::Integer; min_dim::Int=1)

Calls Eirene.betticurve for 'dim' in range from `min_dim` up to 'max_dim' and
stack the resulting Arrays into a vector.

The returned value is a Vector of Arrays{Float64,2}. Each array is of size
(n,2), where n is the maximal number of steps taken to compute Betti curve of dimensions
ranging form `min_dim` to `max_dim`. First column of each array contains numbered steps.
Second column are the Betti curve values for corresponding step.

Arrays in returned vector correspond to Betti curve dimensions form range
`min_dim` up to 'max_dim'.

"""
function get_bettis(results_eirene::Dict, max_dim::Integer; min_dim::Int = 1)
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
# TODO add get_bettis_from_matrix, to wrap C= eirene...; get bettis

#%%

"""
  normalise_bettis(bettis::Vector)
    normalise_bettis(bettis::Array)

Normalise the number of steps for Betti curves. 'bettis' can be either vector of
arrays (each array contain Betti curve of different dimension) or an array
containing Betti curve of a single dimension.
"""
function normalise_bettis(bettis::Vector)

    @debug "Vector version"
    norm_bettis = copy(bettis)
    @debug "norm_bettis size :" size(norm_bettis)[1][1]

    max_dim = size(norm_bettis)[1]
    @debug "typeof(max_dim) :" typeof(max_dim[1])

    for d = 1:(max_dim)
        if !isempty(norm_bettis[d])
            norm_bettis[d][:, 1] /= findmax(norm_bettis[d][:, 1])[1]
        end
    end
    return norm_bettis

end

#%%
function normalise_bettis(bettis::Array)
    @debug "Array version"
    norm_bettis = copy(bettis)
    @debug "norm_bettis size :" size(norm_bettis)

    if !isempty(norm_bettis)
        norm_bettis[:, 1] /= findmax(norm_bettis[:, 1])[1]
    end
    return norm_bettis
end

#%%
# function vectorize_bettis(betti_curves::Array{Matrix{Float64,2}})
"""
    vectorize_bettis(betti_curves::Matrix{Float64})

Reshapes the 'betti_curves' from type Array{Matrices{Float64,2}} into
Matrix{Float64}.

The resulting matrix size is (n, k), where 'n' is equal to the number of
rows in each matrix, 'k' is equal to the number of matrices.

TODO: Change the name- it takse vector and returns a matrix.
TODO: get bettis could have an arguent betti_type which would determine resulting type
"""
function vectorize_bettis(betti_curves::Vector{Array{Float64,2}})

    first_betti = 1
    last_betti = size(betti_curves,1)
    return hcat([betti_curves[k][:, 2] for k = first_betti:last_betti]...)
end

#%%
@deprecate vectorize_bettis(eirene_results::Dict, maxdim::Integer, mindim::Integer) vectorize_bettis(betti_curves)

# ===
#%%
"""
  get_vectorized_bettis(results_eirene::Dict, max_dim::Integer; min_dim::Int=1)

Takes the eirene result and computes Betti curves for dimensions in range
'mindim:maxdim'. Every Betti curve is stored in successive column of the
resulting array.
TODO: this should return a matrix, where first col are indices and rest are B values (1st col is missing now)
"""
function get_vectorized_bettis(results_eirene::Dict, max_dim::Integer; min_dim::Int = 1)

    all_bettis = get_bettis(results_eirene, max_dim, min_dim = min_dim)
    bettis_vector = vectorize_bettis(all_bettis)

    return bettis_vector
end


#%%
"""
    find_max_betti(bettis_collection::Array)

Returns the highest Betti curve value from all dimensions.
"""
function find_max_betti(bettis_collection::Array)
    if typeof(bettis_collection) == Vector
        bettis_collection = vectorize_bettis(bettis_collection)
    end

    max_y_val = 0
    for betti_set in bettis_collection
        local_max = findmax(betti_set)[1]
        if local_max > max_y_val
            max_y_val = local_max
        end
    end
    return max_y_val
end

#%%

"""
    get_max_bettis(bettis)

Returns the maximal bettis of Betti curves for all dimensions.
"""
function get_max_bettis(bettis)

    all_max_bettis = findmax(bettis, dims=1)[1]

    return all_max_bettis
end

# TODO change name
# TODO check what for dim is used, change to min dim
function get_max_betti_from_collection(bettis_collection; dim = 1)
    max_betti = 0
    for betti in bettis_collection
        local_max = findmax(betti)[1]
        if (local_max > max_betti)
            max_betti = local_max
        end
    end
    return max_betti
end


#%%
@deprecate bettis_eirene(matr, maxdim; mintime = -Inf, maxtime = Inf, numofsteps = Inf, mindim = 1) get_bettis(results_eirene, max_dim; min_dim = 1)
@deprecate get_bettis_from_image(img_name, plot_params; file_path = "", plot_heatmaps = true, save_heatmaps = false, plot_betti_figrues = true) get_bettis(results_eirene, max_dim; min_dim = 1)
@deprecate get_bettis_from_image2(img_name;file_path = "",plot_heatmaps = true, save_heatmaps = false, plot_betti_figrues = true) get_bettis_from_image(img_name, plot_params; file_path = "", plot_heatmaps = true, save_heatmaps = false, plot_betti_figrues = true)
@deprecate plot_and_save_bettis2(eirene_results, plot_title::String, results_path::String; file_name = "", extension = ".png", data_size::String = "", do_save = true, extend_title = true, do_normalise = true, min_dim = 0, max_dim = 3, legend_on = true) plot_and_save_bettis(bettis, plot_title::String, results_path::String; file_name = "", extension = ".png", do_save = true, do_normalise = true, min_dim = 0, max_dim = 3, legend_on = true, kwargs...)
@deprecate get_and_plot_bettis(eirene_results; max_dim = 3, min_dim = 1, plot_title = "", legend_on = false) get_bettis(results_eirene, max_dim; min_dim = 1)

#%%
"""
  lower_ordmat_resolution(ordered_matrix::Array, total_bins::Int)

Takes ordered matrix 'input_matrix' and reduces the resolution of values in the
matrix into 'total_bins' bins.
"""
function lower_ordmat_resolution(ordered_matrix::Array, total_bins::Int)
    new_ordered_matrix = zeros(size(ordered_matrix))
    max_val = findmax(ordered_matrix)[1]
    min_val = findmin(ordered_matrix)[1]

    bin_step = max_val ÷ total_bins
    old_bins = min_val:bin_step:max_val

    for bin = 1:total_bins
        @debug "First step threshold is $(old_bins[bin])"
        indices = findall(x -> (x >= old_bins[bin]), ordered_matrix)
        new_ordered_matrix[indices] .= bin - 1
    end

    @debug "Max_val in new matrix is " findmax(new_ordered_matrix)
    @debug "And should be " total_bins - 1

    return new_ordered_matrix
end


#%%
"""
  average_bettis(bettis_matrix; up_factor=8)

Takes the average values of betti curves stored in 'bettis_matrix'.

'bettis_matrix' consist of different simulations(first index of the matrix),
different ranks (third index of the matrix). Second index of the matrices
(saples) may vary accross many simulations and for this reason, all betti curves
are upsampled by a factor of 'upsample_factor' and then the average for every
dimension is computed.
"""
function average_bettis(bettis_matrix::Matrix; up_factor = 8)

    bettis_matrix_backup = copy(bettis_matrix)

    simulations = size(bettis_matrix, 1)
    dimensions = size(bettis_matrix[1], 1)

    max_samples = 0
    for k = 1:simulations
        current_len = length(bettis_matrix[k][1][:, 1])
        if max_samples < current_len
            max_samples = current_len
        end
    end

    bettis_size = size(bettis_matrix)


    total_upsamples = (max_samples - 1) * up_factor + 1
    x_resampled = range(0, 1, step = total_upsamples)

    avg_bettis = zeros(total_upsamples, dimensions)
    std_bettis = copy(avg_bettis)
    resampled_bettis = zeros(simulations, total_upsamples, dimensions)

    # resample betti curves
    for simulation = 1:simulations, betti = 1:dimensions
        resampled_bettis[simulation, :, betti] =
            upsample_vector2(bettis_matrix[simulation][betti][:, 2], total_upsamples)
    end

    # average and std Betti
    for dimension = 1:dimensions
        avg_bettis[:, dimension] = mean(resampled_bettis[:, :, dimension], dims = 1)
        std_bettis[:, dimension] = mean(resampled_bettis[:, :, dimension], dims = 1)
    end

    return avg_bettis, std_bettis
end

#%%
function upsample_vector2(input_vector, total_upsamples)
    total_orig_samples = size(input_vector, 1) - 1

    x_vals = range(0, 1, length = total_orig_samples + 1)
    spl = Spline1D(x_vals, input_vector)

    x_upsampled = range(0, 1, length = total_upsamples)
    y_upsampled = spl(x_upsampled)

    return y_upsampled
end


#%%
"""
  upsample_vector(input_vector; upsample_factor::Int=8)

Takes an 'input_vector' and returns a vector which has 'upsample_factor' many
times more samples. New samples are interpolated with 'spl' function from
'Dierckx' package.

"""
function upsample_vector(input_vector; upsample_factor::Int = 8)
    total_orig_samples = size(input_vector, 1) - 1
    total_samples = upsample_factor * total_orig_samples + 1

    x_vals = range(0, 1, length = total_orig_samples + 1)
    spl = Spline1D(x_vals, input_vector)

    x_upsampled = range(0, 1, length = total_samples)
    y_upsampled = spl(x_upsampled)

    return y_upsampled
end


# =========--=======-========-==========-=======-
# From bettis areas
# Area under Betti curve functions
#%%
"""
    get_area_under_betti_curve(betti_curves, min_dim, max_dim)

Computes the area under Betti curves stored in 'betti_curves', where each row is
a Betti curve and each column is a value.
"""
function get_area_under_betti_curve(betti_curves::Union{Matrix{Float64}, Array{Array{Float64,2}}};do_normalised::Bool=false)
    #TODO check this part
    if size(betti_curves,2) < 2
        bettis_vector = vectorize_bettis(betti_curves)
    else
        bettis_vector = betti_curves
    end
    bettis_area = sum(bettis_vector, dims=1)

    if do_normalised
        total_steps = size(bettis_vector,1)
        bettis_area ./= total_steps
    end
    return bettis_area
end

@deprecate get_area_under_betti_curve(C, min_dim, max_dim) get_area_under_betti_curve(betti_curves; do_normalised=false)

#%%
"""
    get_dataset_bettis_areas(dataset; min_dim::Integer=1, max_dim::Integer=3, return_matrix::Bool=true)

Computes topology of every matrix in dataset, computes Betti curves for dimensions
min_dim up to max_dim and returns vector (or matrix) of areas under Betti curves.
"""
function get_dataset_bettis_areas(dataset; min_dim::Integer=1, max_dim::Integer=3, return_matrix::Bool=true)
    areas_vector = Array[]
    for data = dataset
        @info "Computing topology."
        C = eirene(data, maxdim=max_dim,)
        matrix_bettis = get_bettis(C,max_dim, min_dim=min_dim)
        push!(areas_vector, get_area_under_betti_curve(matrix_bettis))
    end
    if return_matrix
        return vcat([areas_vector[k] for k=1:10]...)
    else
        return areas_vector
    end
end
