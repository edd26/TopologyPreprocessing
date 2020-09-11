using LinearAlgebra
using StatsBase

"""
    shift_to_non_negative(matrix::Matrix)

Returns a matrix in which values are non-negative. This is done by finding the
minimal value in the input matrix and adding its absolute value to the matrix
elements.
"""
function shift_to_non_negative(matrix::Array)

    min_val = findmin(matrix)[1]
    if min_val < 0
        return matrix .-= min_val
    else
        return matrix
    end
end


"""
    normalize_to_01(matrix, norm_factor=256)

Returns a matrix which values are in range [0, 1]. If 'use_factor' is set to
'true' then values are normalized to 'norm_factor' (by default set to 256).

If the values in the input matrix are below 0, then they are shifted so that only positive numbers are in
the matrix (the values are normalized to new maximal value or norm_factor).
"""
function normalize_to_01(matrix::Array; use_factor=false, norm_factor=256)
    normalized_matrix = copy(matrix)

    min_val = findmin(normalized_matrix)[1]
    if min_val < 0
        normalized_matrix .+= abs(min_val)
    else
        normalized_matrix .-= abs(min_val)
    end

    max_val = findmax(normalized_matrix)[1]

    if use_factor
        if max_val > norm_factor
            @warn "Maximal values exceed \'norm_factor\'."
        end
        normalized_matrix = normalized_matrix./norm_factor
    else
        normalized_matrix = normalized_matrix./max_val
    end

    return normalized_matrix
end

"""
    function diagonal_symmetrize(image::Matrix; below_over_upper=false)

Takes an 'image' in the form of a matrix and return a copy which is symmetric
with respect to diagonal- values above diagonal are copied over values below the
diagonal. This can be inverted by setting 'below_over_upper=true'.

If the input matrix is not square, then square matrix is created by taking
matrix of k times 'k' elements, 'k=min(r,c)' where 'r' is number of rows and 'c'
is number of columns.
"""
# function symmetrize_image(image)
function diagonal_symmetrize(image::Matrix; below_over_upper=false)
  w, h = size(image)
  mat_size = findmin([w,h])[1]

  img = copy(image[1:mat_size, 1:mat_size])

  # Get all cartesian indices from input matrix
  matrix_indices = CartesianIndices((1:mat_size, 1:mat_size))
  # Filter out indices below diagonal
  if below_over_upper
      matrix_indices = findall(x->x[1]>x[2], matrix_indices)
  else
        matrix_indices = findall(x->x[2]>x[1], matrix_indices)
    end


  # how many elements are above diagonal
  repetition_number = Int(ceil((mat_size * (mat_size-1))/2))

  # Substitute elements
  for k=1:repetition_number
      # next_position = matrix_indices[k]
      matrix_index = matrix_indices[k]
      # ordered_matrix[matrix_index] = k
      img[matrix_index[2], matrix_index[1]] = img[matrix_index]
  end

  try
      checksquare(img)
  catch err
      if isa(err, DimensionMismatch)
          @error "Resulting matrix is not a square matrix"
          throw(err)
      end
  end
  # issymmetric(Float64.(img))
  return img
end


# =====
# matrix ordering

"""
    get_ordered_matrix(input_matrix; assign_same_values=false, force_symmetry=false,
                            small_dist_grouping=false,
                            min_dist=eps())

Takes a @input_matrix and returns ordered form of this matrix.
The ordered form is a matrix which elements represent ordering from smallest to
highest values in @input_matrix.

If @input_matrix is symmetric, then ordering happens only with upper diagonal.
Lower diagonal is symetrically copied from values above diagonal.

By default, if there is a geoup of entriess with the same value, they all are
assigned with the same ordering number. This can be changed with
@assign_same_values parameter.

Symetry ordering can be froced with @force_symmetry parameter.

By setting 'small_dist_grouping' to true, all the values that difference is
lower than 'min_dist', will be assigned with the same order number.

# Examples
```julia-repl
julia> a = [0 11 12;
            11 0 13;
            12 13 0];
julia> get_ordered_matrix(a)
3×3 Array{Int64,2}:
 0  1  2
 1  0  3
 2  3  0
```

```julia-repl
julia> b = [38 37 36 30;
            37 34 30 32;
            36 30 31 30;
            30 32 30 29]
julia> get_ordered_matrix(b; assign_same_values=false)
4×4 Array{Int64,2}:
0  6  5  2
6  0  1  4
5  1  0  3
2  4  3  0

julia> get_ordered_matrix(b; assign_same_values=true)
4×4 Array{Int64,2}:
0  4  3  1
4  0  1  2
3  1  0  1
1  2  1  0
```
"""
function get_ordered_matrix(in_matrix;
                                assign_same_values=false,
                                force_symmetry=false,
                                small_dist_grouping=false,
                                min_dist=1e-16,
                                total_dist_groups=0)
    # TODO Symmetry must be forced for matrix in which there are NaN elements- needs
    #   to be further investigated
    # TODO not working for negative only values

    # ===
    #
    # ===

    # ==
    if length(size(in_matrix)) > 2
        ord_mat = zeros(Int, size(in_matrix))

        mat_sizes = size(in_matrix)
        min_val = findmin(mat_sizes)[1]
        max_val = findmax(mat_sizes)[1]

        total_min_vals = length(findall(x->x==min_val,mat_sizes))
        total_max_vals = length(findall(x->x==max_val,mat_sizes))
        if total_min_vals > total_max_vals && total_min_vals >= 2
            mat_size = size(in_matrix,findall(x->x==min_val,mat_sizes)[1])
            slices = size(in_matrix,findall(x->x==max_val,mat_sizes)[1])
        elseif total_min_vals < total_max_vals && total_max_vals >= 2
            mat_size = size(in_matrix,findall(x->x==max_val,mat_sizes)[1])
            slices =  size(in_matrix,findall(x->x==min_val,mat_sizes)[1])
        else
            error("Can not resolve matrix size")
        end
    else
        mat_size = size(in_matrix,1)
        ord_mat = zeros(Int, mat_size, mat_size)
        mat_size = size(in_matrix,1)
        slices = 1
    end

    for m=1:slices
        ordered_matrix = ord_mat[:,:,m]
        input_matrix = in_matrix[:,:,m]

        if issymmetric(input_matrix) || force_symmetry
            symetry_order = true
        else
            symetry_order = false
            @warn "Doing non-symetric ordering"
        end

        # distance_groups !=0 && group_distances!(input_matrix, distance_groups)

        # ====
        matrix_indices = generate_indices(mat_size, symmetry_order=symetry_order)


        # Get number of elements to be ordered
        if symetry_order
            # how many elements are above diagonal
            repetition_number = Int(ceil((mat_size * (mat_size-1))/2))
        else
            # how many elements are in whole matrix
            repetition_number = Int(ceil((size(input_matrix)[1] * size(input_matrix)[1])))
        end

        # Get all values which will be sorted
        matrix_values_for_sort = input_matrix[matrix_indices]

        # Sort indices by values (highest to lowest)
        # Create a list of indices, which corresponding valeus are ordered
        sorted_indices = sort!([1:repetition_number;],
                            by=i->(matrix_values_for_sort[i],matrix_indices[i]))

        ordering_number = 0
        for k=1:repetition_number
            # global ordering_number
            next_position = sorted_indices[k]
            matrix_index = matrix_indices[next_position]

            if assign_same_values && k!=1
                prev_index = sorted_indices[k-1]
                prev_matrix_index = matrix_indices[prev_index]

                conditioin1 = input_matrix[prev_matrix_index] == input_matrix[matrix_index]
                conditioin2 = small_dist_grouping
                conditioin3 = abs(input_matrix[prev_matrix_index] - input_matrix[matrix_index]) > min_dist

                if conditioin1 || (conditioin2 && conditioin3)
                    ordered_matrix[matrix_index] = ordering_number-1
                    ordered_matrix[matrix_index[2], matrix_index[1]] = ordering_number-1
                else
                    ordered_matrix[matrix_index] = ordering_number
                    ordered_matrix[matrix_index[2], matrix_index[1]] = ordering_number
                    ordering_number+=1
                end
            else
                ordered_matrix[matrix_index[1], matrix_index[2]] = ordering_number
                ordered_matrix[matrix_index[2], matrix_index[1]] = ordering_number
                ordering_number+=1
            end
        end

        # ====
        non_zero_input = findall(x->x!=0,input_matrix)
        if isempty(ordered_matrix)
            @warn "Ordered matrix is empty"
            min_orig = findmin(input_matrix[non_zero_input])[2]
            max_new = findall(x->x==1,ordered_matrix)[1]
            @debug "Original minimal value was at position: " non_zero_input[min_orig]
            @debug "After ordering the first index value is at position: " max_new
        elseif !isempty(non_zero_input)
        else
            @warn "All values in input matrix were zeros!"
        end
        ord_mat[:,:,m] = ordered_matrix
    end # for loop

    return ord_mat
end

function get_ordered_matrix2(input_array::Array; do_slices = true, dims = 0)
    arr_size = size(input_array)
    out_arr = zeros(Int, arr_size)

    if do_slices
        # param check
        if dims > length(arr_size)
            throw(DomainError("Given dimension is greater than total size of array."))
        elseif dims>0
            throw(DomainError("Given dimension must be positive value."))
        elseif dims<=3
            throw(DomainError("Given dimension must be lower than 3."))
        end

        for dim in 1:arr_size[dims]
            if dims == 1
                out_arr[dim,:,:] = get_ordered_matrix2(input_array[dim,:,:])
            elseif dims == 2
                out_arr[:,dim,:] = get_ordered_matrix2(input_array[:,dim,:])
            elseif dims == 3
                out_arr[:,:,dim] = get_ordered_matrix2(input_array[:,:,dim])
            end
        end
    else
        out_arr = get_ordered_matrix2(input_array)
    end
end

function get_ordered_matrix2(input_array::Array)
    out_array = copy(input_array)
    arr_size = size(input_array)
    total_elements = length(input_array)

    # Collect vector of indices
    all_ind_collected = collect(reshape(generate_indices(arr_size), (length(input_array))))

    # Sort indices vector according to inpu array
    index_sorting = sort!([1:total_elements;],
                        by=i->(input_array[all_ind_collected][i],all_ind_collected[i]))

    for k = 1:total_elements
        target = index_sorting[k]
        out_array[target] = k
    end

    return out_array
end



"""
    function group_distances!(input_matrix, total_dist_groups)

Takes a matrix and rearranges values into 'total_dist_groups' number of groups.
Every group is assigned with number value from range '<0,1>'.

"""
# Care must be taken so that values from 'input_matrix' are within distance
# groups, otherwise error is thrown.
function group_distances(input_matrix::Array, total_dist_groups::Int)
    normed_matrix = normalize_to_01(input_matrix)
    target_matrix = copy(normed_matrix)

    h,w = size(input_matrix)

    if h*w < total_dist_groups
        throw(DomainError("Total number of groups exceed total number of entries in input matrix"))
    end

    total_borders = total_dist_groups+1

    range_val = collect(range(0, 1, length=total_borders))

    for k = 2:total_borders
        indices = findall(x->x>=range_val[k-1] && x<=range_val[k], normed_matrix)
        target_matrix[indices] .= range_val[k]
    end
    unique(target_matrix)

    # Sets last range to values smaller than unity, just in case this might cause trobules
    # normed_matrix[normed_matrix .> range_val[end-1]] .= 0.99
    return target_matrix
end


# function group_distances(input_matrix::Matrix, total_dist_groups::Int)
#     new_matrix = copy(input_matrix)
#     group_distances!(new_matrix, total_dist_groups)
#     return new_matrix
# end

# function normalize_distances!(input_matrix)
#     max_val = findmax(input_matrix)[1]
#     min_val = findmin(input_matrix)[1]
#
#     if min_val < 0
#         input_matrix .+= abs(min_val)
#     else
#         input_matrix .-= abs(min_val)
#     end
#     input_matrix ./= abs(max_val)
# end

"""
    function generate_indices(matrix_size; symetry_order=false)

Return all the possible indices of the matrix of size 'matrix_size'. If
'symetry_order' is set to'true', then only indices of values below diagonal are
returned.
"""
function generate_indices(matrix_size; symmetry_order=false, include_diagonal=false)
    # Get all cartesian indices from input matrix
    matrix_indices = CartesianIndices((1:matrix_size, 1:matrix_size))
    # Filter out indices below diagonal
    if symmetry_order
        if include_diagonal
            matrix_indices = findall(x->x[1]<=x[2], matrix_indices)
        else
            matrix_indices = findall(x->x[1]<x[2], matrix_indices)
        end
    else
        matrix_indices = findall(x->true, matrix_indices)
    end
    return matrix_indices
end


"""
    function generate_indices(dims::Tuple)

Generate indices for a matrix of given dimensions.
"""
function generate_indices(dims::Tuple)
    total_dims = length(dims)
    matrix_indices = CartesianIndices(dims)
end

dims = (12,4,6)

# matrix ordering
# =====


function get_high_dim_ordered_matrix(input_matrix)
    matrix_size = size(input_matrix)
    ordered_matrix_3D = zeros(Int, matrix_size)

    for slice = 1:matrix_size[1]
        ordered_matrix_3D[slice,:,:] = get_ordered_matrix(input_matrix[slice, :, :])
    end
    return ordered_matrix_3D
end


"""
    reduce_arrs_to_min_len(arrs)

Takes vector of vectors of different length and returns array of arrays which
are of the same length. Length in the output is the shortest vector length from
the input- values above this size are discarded.
"""
function reduce_arrs_to_min_len(arrs::Array)
    new_arr = copy(arrs)

    simulation = size(new_arr,1)
    min_size = Inf
    for m=1:simulation
        @debug "Simulation number" m
        current_size = size(new_arr[m],1)
        @debug "Current size: " current_size
        if convert(Float64,current_size) < min_size
            min_size = current_size
            @debug "min size changed to: " min_size
        end
    end
    # min_size = Int.(min_size)
    @debug "Concatenating"
    for m=1:simulation
        new_arr[m] = new_arr[m][1:min_size,:]
    end
    min_size = Inf
    return new_arr
end


"""
    increase_arrs_to_max_len(arrs)

Takes vector of vectors of different length and returns array of arrays which
are of the same length. Length in the output is the longest vector length from
the input- values above this size are discarded.
"""
function increase_arrs_to_max_len(arrs)
    new_arr = copy(arrs)

    simulation = size(new_arr,1)
    max_size = 0
    for m=1:simulation
        @debug "Simulation number" m
        current_size = size(new_arr[m],1)
        @debug "Current size: " current_size
        if convert(Float64,current_size) > max_size
            max_size = current_size
            @debug "min size changed to: " max_size
        end
    end
    # max_size = Int.(max_size)
    @debug "Concatenating"
    for m=1:simulation
        correct_len_arr = zeros(Int, max_size, 3)
        correct_len_arr[1:size(arrs[m],1),:] = new_arr[m][:,:]
        new_arr[m] = correct_len_arr
    end
    # min_size = Inf
    return new_arr
end
