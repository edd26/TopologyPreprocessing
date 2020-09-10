using LinearAlgebra
using StatsBase
"""
    shift_to_non_negative(matrix)

Returns a matrix in which values are non-negative. This is done by finding the
minimal value in the input matrix and adding its absolute value to the matix
elements.
"""
function shift_to_non_negative(matrix)

    min_val = findmin(matrix)[1]
    if min_val < 0
        return matrix .-= min_val
    else
        return matrix
    end
end


"""
    normalize_to_01(matrix, norm_factor=256)

Returns a matrix which values are in range [0, 1]. If the values in the input
matrix are below 0, then they are shifted so that only positive numbers are in
the matrix. If the values in the matrix of shifted matrix excced value of the
@norm_factor parameter, then the matrix is normalized to the maximal value from
the matrix.
"""
function normalize_to_01(matrix; norm_factor=256)
    normalized_matrix = copy(matrix)
    min_val = findmin(normalized_matrix)[1]
    max_val = findmax(normalized_matrix)[1]

    if min_val < 0
        normalized_matrix .+= abs(min_val)
    end

    if max_val > norm_factor
        @warn "Values normalized to maximal value, not notmalization factor."
        normalized_matrix = normalized_matrix./max_val
    else
        normalized_matrix = normalized_matrix./norm_factor
    end
    return normalized_matrix
end

"""
    function symmetrize(image)

Takes an @image and return a copy which is symmetric.
"""
function symmetrize_image(image)
  mat_size = size(image,1)

  img= copy(image)
  # Get all cartesian indices from input matrix
  matrix_indices = CartesianIndices((1:mat_size, 1:mat_size))
  # Filter out indices below diagonal
  matrix_indices = findall(x->x[1]>x[2], matrix_indices)

  # Put evrything together
  # how many elements are above diagonal
  repetition_number = Int(ceil((mat_size * (mat_size-1))/2))

  for k=1:repetition_number
      # next_position = matrix_indices[k]
      matrix_index = matrix_indices[k]
      # ordered_matrix[matrix_index] = k
      img[matrix_index[2], matrix_index[1]] = img[matrix_index]
  end
  issymmetric(Float64.(img))
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
                                distance_groups=0)
    # TODO Symmetry must be forced for matrix in which there are NaN elements- needs
    #   to be further investigated
    # TODO not working for negative only values

    # ==
    if length(size(in_matrix)) > 2
        ord_mat = zeros(Int, size(in_matrix))

        mat_sizes = size(in_matrix)
        min_val = findmin(mat_sizes)[1]
        max_val = findmax(mat_sizes)[1]

        min_vals_occurence = length(findall(x->x==min_val,mat_sizes))
        max_vals_occurence = length(findall(x->x==max_val,mat_sizes))
        if min_vals_occurence > max_vals_occurence && min_vals_occurence >= 2
            mat_size = size(in_matrix,findall(x->x==min_val,mat_sizes)[1])
            slices = size(in_matrix,findall(x->x==max_val,mat_sizes)[1])
        elseif min_vals_occurence < max_vals_occurence && max_vals_occurence >= 2
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
    # ==

        if issymmetric(input_matrix) || force_symmetry
            symetry_order = true
        else
            symetry_order = false
            @warn "Doing non-symetric ordering"
        end

        # distance_groups !=0 && group_distances!(input_matrix, distance_groups)

        # ====
        matrix_indices = generate_indices(mat_size, symetry_order)


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

function group_distances!(input_matrix, distance_groups)
    normalize_distances!(input_matrix) # Probably normalized_matrix would suffice
    distance_bins = distance_groups+1

    range_val = range(0, 1, length=distance_bins)

    for k = 2:distance_bins
        indices = findall(x->x>range_val[k-1] && x<range_val[k], input_matrix)
        input_matrix[indices] .= range_val[k]
    end
    unique(input_matrix)

    # Sets last range to values smaller than unity, just in case this might cause trobules
    input_matrix[input_matrix .> range_val[end-1]] .= 0.99
end

function normalize_distances!(input_matrix)
    max_val = findmax(input_matrix)[1]
    min_val = findmin(input_matrix)[1]

    if min_val < 0
        input_matrix .+= abs(min_val)
    else
        input_matrix .-= abs(min_val)
    end
    input_matrix ./= abs(max_val)
end

function generate_indices(matrix_size, symetry_order)
    # Get all cartesian indices from input matrix
    matrix_indices = CartesianIndices((1:matrix_size, 1:matrix_size))
    # Filter out indices below diagonal
    if symetry_order
        matrix_indices = findall(x->x[1]<x[2], matrix_indices)
    else
        matrix_indices = findall(x->true, matrix_indices)
    end
    return matrix_indices
end

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
    get_pairwise_correlation_matrix(vectorized_video, tau_max=25)

Computes pairwise correlation of the input signals accordingly to the formula
presented in paper "Clique topology reveals intrinsic geometric structure
in neural correlations" by Chad Giusti et al.

The Computations are done only for upper half of the matrix, the lower half is
a copy of upper half. Computation-wise the difference is at level of 1e-16, but
this causes that inverse is not the same as non-inverse matrix.

"""
function get_pairwise_correlation_matrix(vectorized_video, tau_max=25)
    number_of_signals = size(vectorized_video,1)
    T = size(vectorized_video,2)

    C_ij = zeros(number_of_signals,number_of_signals);
    # log_C_ij = zeros(number_of_signals,number_of_signals);

     # this is given in frames
    lags = -tau_max:1:tau_max


    for row=1:number_of_signals
        for column=row:number_of_signals
            signal_ij = vectorized_video[row,:];
            signal_ji = vectorized_video[column,:];

            # cross_corelation
            ccg_ij = crosscov(signal_ij, signal_ji, lags);
            ccg_ij = ccg_ij ./ T;

            A = sum(ccg_ij[tau_max+1:end]);
            B = sum(ccg_ij[1:tau_max+1]);
            r_i_r_j = 1;
            C_ij[row, column] = max(A, B)/(tau_max*r_i_r_j);
            C_ij[column, row] = C_ij[row, column]
            # log_C_i_j[row, column] = log10(abs(C_ij[row, column]));
        end
    end

    return C_ij
end



"""
    get_subimg_correlations(video_array, centers, sub_img_size, shift)

Computes the correlation between the subimages and subimages shifted by values
from range -@shift:@shift and returns array with frames of size
length(@centers) x length(@centers) with the number of frames equal to the
number of rames in @video_array.

Each of the subimage is center around values stored in  @centers
"""
function get_subimg_correlations(video_array, centers, sub_img_size, shift)
    half_size = ceil(Int,(sub_img_size-1)/2)
    half_range = half_size + shift
    h, w, len = get_video_dimension(video_array)
    extracted_pixels = zeros(sub_img_size, sub_img_size, len)

    for frame = 1:len
        img = video_array[frame]
        for index_x = 1:size(centers,2)
            c_x = centers[2, index_x]
            for index_y = 1:size(centers,2)
                c_y = centers[1, index_y]
                subimage = img[(c_x-half_range):(c_x+half_range),
                                (c_y-half_range):(c_y+half_range)]
                center = img[(c_x-half_size):(c_x+half_size), (c_y-half_size):(c_y+half_size)]

                for left_boundary = 1:(2*shift+1)
                    for lower_boundary = 1:(2*shift+1)
                        corelation = center .* subimage[left_boundary:left_boundary+sub_img_size-1, lower_boundary:lower_boundary+sub_img_size-1]
                        corelation = sum(corelation)
                        extracted_pixels[index_x, index_y, frame] += corelation
                    end
                end
                extracted_pixels[index_x, index_y, frame] /= 256*(sub_img_size^2)*(shift*2)^2
            end
        end
    end
    return extracted_pixels
end



"""
    reduce_arrs_to_min_len(arrs)

Takes vector of vectors of different length and returns array of arrays which
are of the same length. Length in the output is the shortest vector length from
the input- values above this size are discarded.
"""
function reduce_arrs_to_min_len(arrs)
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
are of the same length. Length in the output is the shortest vector length from
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
