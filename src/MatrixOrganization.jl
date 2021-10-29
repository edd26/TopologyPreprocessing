using LinearAlgebra,
using Plots
using Random

include("PlottingWrappers.jl")
include("PointsSubstitution.jl")

"""
	function expand_matrix(input_matrix, expansion_size, last_components;do_plot=false)

Takes 'input_matrix' (an ordering matrix used for creating cliques) and and adds
2×'expansion_size' number of rows. 'last_components' are the values in original
matrix that are added last to the clique.

Results may be be plotted by setting 'do_plot=true'.
"""
function expand_matrix(input_matrix, expansion_size, last_components;do_plot=false)
	new_comp = last_components
	matrix_size = size(input_matrix,1)
	for mat_sizes = matrix_size:2:(matrix_size+2expansion_size)
		input_matrix, new_comp = add_step_to_matrix(input_matrix, new_comp)
	end
	if do_plot
# TODO separate plotting from processing
		expand_plt_ref = plot_square_heatmap(input_matrix, 1,size(input_matrix,1);
											plt_title = "Original, size:$(matrix_size)",
											color_palete=:lightrainbow)
		display(expand_plt_ref)
	end
	return input_matrix, expand_plt_ref
end

# Shuffle matrix entries
"""
	function shuffle_matrix(input_matrix, shuffles; do_plot=false)

Takes symmetric 'input_matrix' and randomly swaps rows 'shuffles' many times.

Results may be plotted by setting 'do_plot=true'.
"""

function shuffle_matrix(input_matrix, shuffles; do_plot=false)
	matrix_size = size(input_matrix,1)
	rows = randcycle(matrix_size)
	shuffled_ord_mat = copy(input_matrix)

	for k = 1:shuffles
		# global shuffled_ord_mat, rows
		srcs, trgts = rand(rows,2)
		swap_rows!(shuffled_ord_mat, srcs, trgts)
	end

	if do_plot
# TODO separate plotting from processing
		shuff_plt_ref = plot_square_heatmap(shuffled_ord_mat, 1,size(shuffled_ord_mat,1);
											plt_title = "Shuffled, size:$(matrix_size)",
											color_palete=:lightrainbow)
		display(shuff_plt_ref)
	end
	return input_matrix, shuff_plt_ref
end



"""
	function organize_shuff_matrix(input_matrix; do_plots=false)

Reorganizes 'input_matrix' so that values highest values in a row are positioned
next to the diagonal.

Results may be plotted by setting 'do_plot=true'.
"""
function organize_shuff_matrix(input_matrix; do_plots=false)
	unscrambled_matrix = copy(input_matrix)
	matrix_size = size(input_matrix,1)
	for k = matrix_size:-2:2
		max_row_val = findmax(unscrambled_matrix[k,:])[2]
		# put to the previous last position
		swap_rows!(unscrambled_matrix, max_row_val, k-1)
		# skip 1 row and work on next one
	end
	if do_plots
# TODO separate plotting from processing
		reorganized_plt_ref = plot_square_heatmap(unscrambled_matrix, 1,size(unscrambled_matrix,1);
									plt_title = "unscrambled_matrix, size:$(matrix_size)",
									color_palete=:lightrainbow)
		display(reorganized_plt_ref)
		return unscrambled_matrix, reorganized_plt_ref
	else
		return unscrambled_matrix
	end
end


"""
		function order_max_vals_near_diagonal(input_matrix; do_plots=false, direction=:descending)

Orders values in 'input_matrix' so that values next to diagonal are descending
(by default).

TODO- not working-  Optionally, ascending order can be used by setting 'direction' to
':ascending'.

Results may be plotted by setting 'do_plot=true'.
"""
function order_max_vals_near_diagonal(input_matrix; do_plots=false, direction=:descending)
	# Find max values next to the diagonal
	matrix_size = size(input_matrix,1)

	if direction == :descending
		# ordering_function = findmax
		new_ord_value = -1
		iteration_values = matrix_size:-2:2
		# iteration_values = 2:2:matrix_size

	elseif direction == :ascending
		# ordering_function = findmin
		new_ord_value = findmax(input_matrix)[1]*2
		iteration_values = 2:2:matrix_size
	else
		# @error "Unknow ordering was given"
		throw("Unknow ordering was given")
	end

	reordered_matrix = copy(input_matrix)
	row_indices = 1:2:matrix_size
	col_indices = 2:2:matrix_size
	coord_set = [CartesianIndex(row_indices[k], col_indices[k]) for k=1:matrix_size÷2]
	diag_max_values = reordered_matrix[coord_set]

	for k = iteration_values
		max_val, max_ind = findmax(diag_max_values)
		(direction == :descending) ? (position = floor(k÷2)) : (position = floor(k÷2))
		diag_max_values[max_ind] = diag_max_values[position]
		diag_max_values[position] = new_ord_value
		max_ind *= 2

		swap_rows!(reordered_matrix, k, max_ind)
		swap_rows!(reordered_matrix, k-1, max_ind-1)
	end


	if do_plots
# TODO separate plotting from processing
		reorganized_plt_ref = plot_square_heatmap(reordered_matrix, 1,size(reordered_matrix,1);
									plt_title = "reordered_matrix, size:$(matrix_size)",
									color_palete=:lightrainbow)
		display(reorganized_plt_ref)
		return reordered_matrix, reorganized_plt_ref
	else
		return reordered_matrix
	end
end


"""
	function fine_tune_matrix(input_matrix; do_plots=false)

Check if velues next to the maximal values are organized in descending order.
"""
function fine_tune_matrix(input_matrix; do_plots=false)#, direction=:descending)
	# Find max values next to the diagonal
	matrix_size = size(input_matrix,1)
	fine_tune_matrix = copy(input_matrix)

	# if direction == :descending
	# 	# ordering_function = findmax
	# 	new_ord_value = -1
	# 	iteration_values = matrix_size:-2:2
	# 	# iteration_values = 2:2:matrix_size
	#
	# elseif direction == :ascending
	# 	# ordering_function = findmin
	# 	new_ord_value = findmax(input_matrix)[1]*2
	# 	iteration_values = 2:2:matrix_size
	# else
	# 	# @error "Unknow ordering was given"
	# 	throw("Unknow ordering was given")
	# end


	for k = 2:2:matrix_size-1
		if fine_tune_matrix[k-1,k+1] > fine_tune_matrix[k,k+1]
			swap_rows!(fine_tune_matrix, k, k-1)
		end
	end


	if do_plots
# TODO separate plotting from processing
		fine_tuned_plt_ref = plot_square_heatmap(fine_tune_matrix, 1,size(reordered_matrix,1);
										plt_title = "fine_tuned, size:$(matrix_size)",
										color_palete=:lightrainbow)
		display(fine_tuned_plt_ref)
	end
	return fine_tune_matrix, fine_tuned_plt_ref
end

# TODO separate plotting from processing
function order_max_vals_by_row_avg(input_matrix; do_plots=false)
	# Find max values next to the diagonal
	matrix_size = size(input_matrix,1)

	# Row average
	row_avg = reshape(mean(input_matrix, dims=1),(matrix_size, 1))

	# Using funciton below, because sortperm is not working on Array{Float64,2}
	sorted_rows_indexes = [findall(x->x==sort(row_avg, dims=1)[k], row_avg)[1][1] for k=1:matrix_size]
	matrix_indices = collect(range(1,matrix_size))
	# Sort indices by values (highest to lowest)
	# Create a list of indices, which corresponding valeus are ordered
	sorted_indices = sort!([1:matrix_size;],
						by=i->(sorted_rows_indexes[i],matrix_indices[i]))

	sorted_matrix = copy(input_matrix)
	for k = 1:matrix_size÷2 #iteration_values
		max_ind = sorted_indices[k]

		sorted_indices[k] = k
		sorted_indices[max_ind] = max_ind

		swap_rows!(sorted_matrix, k, max_ind)
		# swap_rows!(sorted_matrix, k-1, max_ind-1)
	end
# TODO separate plotting from processing
	reorganized_plt_ref = plot_square_heatmap(sorted_matrix, 1,size(reordered_matrix,1);
							plt_title = "reordered_matrix, size:$(matrix_size)",
							color_palete=:lightrainbow)

# TODO separate plotting from processing
		input_mat_plt_ref = plot_square_heatmap(input_matrix, 1,size(reordered_matrix,1);
									plt_title = "input_matrix, size:$(matrix_size)",
									color_palete=:lightrainbow)

		common_plot1 = plot(input_mat_plt_ref, reorganized_plt_ref, layout=(1,2),
																size=(800,400))


	# reordered_matrix = copy(input_matrix)
	# row_indices = 1:2:matrix_size
	# col_indices = 2:2:matrix_size
	# coord_set = [CartesianIndex(row_indices[k], col_indices[k]) for k=1:matrix_size÷2]
	#
	#
	# for k = iteration_values
	# 	max_val, max_ind = findmax(diag_max_values)
	# 	position = floor(k÷2)
	# 	diag_max_values[max_ind] = diag_max_values[position]
	# 	diag_max_values[position] = new_ord_value
	# 	max_ind *= 2
	#
	# 	swap_rows!(reordered_matrix, k, max_ind)
	# 	swap_rows!(reordered_matrix, k-1, max_ind-1)
	# end


	if do_plots
# TODO separate plotting from processing
		reorganized_plt_ref = plot_square_heatmap(sorted_matrix, 1,size(reordered_matrix,1);
									plt_title = "reordered_matrix, size:$(matrix_size)",
									color_palete=:lightrainbow)
		display(reorganized_plt_ref)
	end
	return reordered_matrix, reorganized_plt_ref
end


function order_max_vals_near_diagonal2(input_matrix; do_final_plot=false, do_all_plots = false, direction=:descending)
	# Find max values next to the diagonal
	matrix_size = size(input_matrix,1)
	reordered_matrix = copy(input_matrix)
	reorganized_plt_ref = []

	# for every row in matrix
	for k = 1:2:matrix_size-1
		# global reordered_matrix
# TODO separate plotting from processing
		reorganized_plt_ref_pt0 = plot_square_heatmap(reordered_matrix, 1,size(reordered_matrix,1);
								plt_title = "reordered_matrix, size:$(matrix_size)",
								color_palete=:lightrainbow)
		max_val, max_ind = findmax(reordered_matrix[k:end, k:end])
		# Take the smaller coordinate
		# (max_ind[1] < max_ind[2]) ? (target_row = max_ind[1]) : (target_row = max_ind[2])
		target_row = max_ind[1]+k-1

		reordered_matrix = swap_rows(reordered_matrix, k, target_row)
# TODO separate plotting from processing
		reorganized_plt_ref_pt1 = plot_square_heatmap(reordered_matrix, 1,size(reordered_matrix,1);
									plt_title = "reordered_matrix, size:$(matrix_size)",
									color_palete=:lightrainbow)
	    val, second_target = findmax(reordered_matrix[k,k:end])
		second_target = second_target+k-1
		reordered_matrix = swap_rows(reordered_matrix, k+1, second_target)
		# end
		#
		#
		if do_all_plots
# TODO separate plotting from processing
			reorganized_plt_ref_pt2 = plot_square_heatmap(reordered_matrix, 1,size(reordered_matrix,1);
										plt_title = "reordered_matrix, size:$(matrix_size)",
										color_palete=:lightrainbow)
			reorganized_plt_ref = plot(reorganized_plt_ref_pt0, reorganized_plt_ref_pt1, reorganized_plt_ref_pt2, layout=(1,3), size=(1400,400))
			display(reorganized_plt_ref)
		end
	end
	if do_final_plot
# TODO separate plotting from processing
		reorganized_plt_ref = plot_square_heatmap(reordered_matrix, 1,size(reordered_matrix,1);
									plt_title = "reordered_matrix, size:$(matrix_size)",
									color_palete=:lightrainbow)
		# display(reorganized_plt_ref)
	else
		reorganized_plt_ref=[]
	end
	return reordered_matrix, reorganized_plt_ref
end


"""
    function get_key_for_value(d::Dict, target_value)

Returns key of the dictionary which corresponds to the given target value.
"""
function get_key_for_value(d::Dict, target_value)

    for (key, value) in d
        if value == target_value
            return key
        end
    end
end

function order_max_vals_near_diagonal3(input_matrix, ordering; direction=:descending)
    # Find max values next to the diagonal
    matrix_size = size(input_matrix,1)
    reordered_matrix = deepcopy(input_matrix)
    new_ordering = Dict()

    # for every row in matrix
    for m = 1:2:matrix_size-1
        # global reordered_matrix

        max_val, max_ind = findmax(reordered_matrix[m:end, m:end])
        # Take the smaller coordinate
        first_target = max_ind[1]+m-1
        reordered_matrix = swap_rows(reordered_matrix, m, first_target)

        # check for duplicates
        val, second_target = findmax(reordered_matrix[m,m:end])
        second_target = second_target+m-1
        if first_target == second_target
            @debug "are same"
            second_target -= 1
        end

        reordered_matrix = swap_rows(reordered_matrix, m+1, second_target)

        # find key which initially had region = first_target
        region1 = get_key_for_value(ordering, first_target)
        region2 = get_key_for_value(ordering, second_target)

        if region1 in keys(new_ordering)
            @warn "repeated"
        end
        if region2 in keys(new_ordering)
            @warn "repeated2"
        end

        new_ordering[region1] = m
        new_ordering[region2] = m +1
        println("Replaced $(region1) fom $(ordering[region1]) to $(m)")
        println("Replaced $(region2) fom $(ordering[region2]) to $(m+1)")
    end

    return reordered_matrix, new_ordering
end

##
# """
#    matrix_poling!(input_matrix; method = "avg_pooling")
#
# Takes a matrix and changes it's values to the same value, according to 'method'.
# Possible methods are:
# - 'max_pooling'- finds maximal value and replaces all values with the maximal
# 	value.
# - 'avg_pooling'- changes values to the average value
# - 'gauss_pooling'- uses gausian kernel as weights to the values in the matrix
# """
# function matrix_poling!(input_matrix::Array; method::String = "max_pooling")
# 	if method == "max_pooling"
# 		max_val = findmax(input_matrix)[1]
# 		input_matrix .= max_val
# 	end
# 	return input_matrix
# end

# matrix_poling(mat[1:3,1:3]; method = "gauss_pooling")
function matrix_poling(input_matrix::Array; method = "avg_pooling", kernel_size=3,gauss_sigma=1)
	out_matrix = copy(input_matrix)
	if method == "max_pooling"
		max_val = findmax(out_matrix)[1]
		out_matrix .= max_val
	elseif method == "avg_pooling"
		avg_val = mean(out_matrix)
		out_matrix .= floor(Int,avg_val)
	elseif method == "gauss_pooling"
		@debug "Gauss pooling"

		# gauss_kernel = ones(kernel_size,kernel_size)
		# gauss_kernel[kernel_size÷2+1,kernel_size÷2+1] = 2
		filtering_kernel = Kernel.gaussian(gauss_sigma)
		# gauss_kernel2 = imfilter(gauss_kernel, filtering_kernel)
		# gauss_kernel3 = gauss_kernel2.-findmin(gauss_kernel2)[1]/findmax(gauss_kernel2)[1]
		# gauss_kernel3 = gauss_kernel3./sum(gauss_kernel3)
		out_matrix = imfilter(out_matrix, filtering_kernel)
		# out_matrix += out_matrix.*gauss_kernel3
		out_matrix .= floor.(Int,out_matrix)
	end
	return out_matrix
end


function subsample_matrix(square_matrix::Array; subsamp_size::Int=2, method="max_pooling")
	if !issymmetric(square_matrix)
		error("Input matrix is not square")
		return
	end

	if subsamp_size == 2
		reorganize_matrix(square_matrix; subsamp_size=subsamp_size, method=method)
		return reorganize_matrix(square_matrix; subsamp_size=subsamp_size, method=method)
	elseif subsamp_size%2 == 0
		new_matrix = reorganize_matrix(square_matrix; subsamp_size=subsamp_size, method=method)
		return reorganize_matrix(new_matrix; subsamp_size=2, method=method)
	else
		new_matrix = reorganize_matrix(square_matrix; subsamp_size=subsamp_size, method=method)
		return new_matrix
	end
end

#= Should the result be overlapping or not? Options:
	- filter it as a whole image, return diagonal
	- filter subimages- the problem is that htere will be edge effect at every border of cells
	- filter subimages and assign midde value to whole patch
	- filter whole upper diagonal matrix

Plot gaussian kernel
Should we care about
=#
function reorganize_matrix(square_matrix::Array; subsamp_size::Int=2, method="max_pooling", overlap::Int=0,gauss_sigma=1)
	if method == "gauss_pooling"
		(subsamp_size >= 3) || error("Can not do gaussian pooling for area smaller than 3x3")
	end
	@debug method
	# Subsample upper half
	square_matrix2 = Float64.(copy(square_matrix))
	total_rows, total_cols = size(square_matrix)
	size_mismatch_flag = false

	# if total_rows%2 != 0
	# 	total_rows -= 1
	# end
	# if total_cols%2 != 0
	# 	total_cols -= 1
	# end

	if method == "gauss_pooling"
		square_matrix2 = zeros(Int,size(square_matrix))
		square_matrix2[1:end-1,2:end] = UpperTriangular(square_matrix[1:end-1,2:end])

		# flip matrix
		do_matrix_flip = true
		if do_matrix_flip
			square_matrix3 = zeros(Float64,size(square_matrix))
			for row in 0:total_rows-1
				for col in 0:total_cols-1
					square_matrix3[row+1,col+1] = square_matrix[end-row,end-col]
				end
			end
			square_matrix3[1:end-1,:] = square_matrix3[2:end,:]
			square_matrix3[:,2:end] = square_matrix3[:,1:end-1]
		else
			square_matrix3 = copy(square_matrix2)
			square_matrix3[1:end-1,:] = square_matrix2[2:end,:]
			square_matrix3[:,1:end-1] = square_matrix2[:,2:end]
		end

		for row in 1:total_rows
			for col in 1:row
				square_matrix2[row,col] = square_matrix3[row,col]
			end
		end

		filtering_kernel = Kernel.gaussian(gauss_sigma)
		square_matrix2 = imfilter(square_matrix2, filtering_kernel)
		square_matrix2 .= Int.(floor.(Int,square_matrix2))
	elseif method == "row_pooling"
		function gauss_func(σ,len;μ=0)
			maxv = len÷2
			 minv= -len÷2
			if len%2 == 0
				maxv-=1
			end
			x = collect(minv:1:maxv)
			return exp.(-(((x.-μ)./σ)./2).^2)./(σ*sqrt(2π))
		end
		# Take 'subsamp_size'in total in horizontal and in vertical line from
		# current matrix element
		# subsamp_size = 5
		val_range = subsamp_size÷2
		r = (subsamp_size÷2)*2
		# total_rows = 266
		# total_cols = 266
		# row = 3

		for row = 1:1:(total_rows-1)
			for col = (row+1):1:total_cols
				if row < r  &&  col <= r
					row_range = row - 1
					col_range = val_range + (val_range-row_range÷2)
				else
					row_range = val_range
				end
				if row > total_rows-r  &&  col >= total_cols-r
					col_range = total_cols - row -1
					row_range = val_range + (val_range-col_range)
				else
					col_range = val_range
				end

				r_beg = row - row_range
				r_end = row + row_range
				c_beg = col - col_range
				c_end = col + col_range

				# if r_beg < 1 && r_end > total_rows
				if r_beg < 1
					r_end += abs(r_beg)+1
					r_beg = 1
				end
				if r_end > col
					r_beg -= abs(r_end-col)
					if r_beg <1
						r_beg=1
					end

					r_end = col-1
				end
				# end # if both

				# if c_beg < row+1 && c_end > total_cols
				if c_beg < row+1
					c_end += abs(c_beg-(row+1))
					c_beg = row+1
				end
				if c_end > total_cols
					c_beg -= abs(total_rows-c_end)
					c_end = total_cols
				end
				vrange = r_beg:r_end
				try
					square_matrix2[row,col] += sum(
												square_matrix[vrange,col]
												.* gauss_func(gauss_sigma,length(vrange))
												)
												vrange = c_beg:c_end
					square_matrix2[row,col] += sum(
													square_matrix[row,c_beg:c_end] .*
													 gauss_func(gauss_sigma,length(vrange))
													)
				catch e
					@error "Failed to compute row pooling"
					@error "row" row
					@error "col" col
					square_matrix2[row,col] = 0
					break
					# error(e)
				end

			end # for col
		end # for rows

	else
		step = subsamp_size-overlap
		for row = 1:step:(total_rows-2)
			for col = (row+subsamp_size):step:total_cols
				r_beg = row
				r_end = row+subsamp_size-1
				c_beg = col
				c_end = col+subsamp_size-1

				if r_end > total_rows || c_end > total_cols
					size_mismatch_flag = true
					continue
				end
				square_matrix2[r_beg:r_end,c_beg:c_end] =
							matrix_poling(square_matrix[r_beg:r_end,c_beg:c_end]; method=method,kernel_size=subsamp_size, gauss_sigma=gauss_sigma)
			end # for col
			size_mismatch_flag && continue
		end # for rows
	end # if method

	# Copy over lower half
	for row in 2:total_rows
		for col in 1:row-1
			square_matrix2[row,col] = square_matrix2[col,row]
		end
	end

	# keep same values on diagonal
	for row in 1:total_rows
		square_matrix2[row,row] = square_matrix[row,row]
	end

	return square_matrix2
end

function pool_matrix(square_matrix::Array; method="max_pooling")
	out_matrix = copy(square_matrix)
	pool_matrix!(square_matrix; method=method)
	return out_matrix
end

"""
	add_random_patch(input_matrix; patch_size=1, total_patches=1, locations)

Takes a matrix and replaces some values with random values. Returns a new matrix
with replaced values and indicies where replacement took place.

Values can be
replaced by setting 'patch_size' to values bigger than 1. If the input matrix
is symmetric, then output matrix will be symmetric as well (values from above
diagnoal will be copied over values from below diagonal).
"""
function add_random_patch(input_matrix::Matrix; patch_size=1, total_patches=1, locations=CartesianIndex(0))
	total_rows, total_cols = size(input_matrix)
	max_row = total_rows-patch_size+1
	max_col = total_cols-patch_size+1

	output_matrix = copy(input_matrix)
	max_val = findmax(output_matrix)[1]
	min_val = findmin(output_matrix)[1]
	matrix_type = typeof(output_matrix[1])

	if patch_size>total_rows || patch_size>total_cols
		error(DimensionMismatch,": Patch size is bigger than the matrix!")
	end

	# ===
	issymmetric(input_matrix) ? (symmetrize_matrix = true) : (symmetrize_matrix = false)

	if locations == CartesianIndex(0)
		@debug "Locations were not specified- random locations will be used"
		if symmetrize_matrix
			possible_indices = findall(x->true,UpperTriangular(output_matrix))
			possible_indices = possible_indices[findall(x->x[1]<=x[2], possible_indices)]
			possible_indices = possible_indices[findall(x->x[1]<=max_row, possible_indices)]
			possible_indices = possible_indices[findall(x->x[2]<=max_col, possible_indices)]
		else
			possible_indices = possible_indices = findall(x->true,output_matrix)
		end

		tartget_indices = possible_indices[randcycle(length(possible_indices))]

	else
		wrong_indices = findall(x->x[1]>max_row || x[2]>max_col, locations)
		if isempty(wrong_indices)
			tartget_indices = locations
			total_patches = size(locations)[1]
		else
			error(DimensionMismatch,": Given indices are bigger than the matrix dimensions!")
		end
	end

	changed_indices = CartesianIndex[]
	for replacement=1:total_patches
		row = tartget_indices[replacement][1]
		col = tartget_indices[replacement][2]
		r_range = row:row+patch_size-1
		c_range = col:col+patch_size-1

		for ind in CartesianIndices((r_range,c_range))
			push!(changed_indices,ind)
		end

		new_rand_matrix = floor.(matrix_type, rand(patch_size,patch_size) .* (max_val-min_val+1) .+ min_val)

		output_matrix[r_range,c_range] .= new_rand_matrix
	end

	if symmetrize_matrix
		# Inverse second column
		changed_indices2 = [changed_indices changed_indices]
		for ind = 1:size(changed_indices)[1]
			c_ind = changed_indices2[ind,2]
			changed_indices2[ind,2] = CartesianIndex(c_ind[2],c_ind[1])
		end

		# Copy over lower half
		for row in 2:total_rows
			for col in 1:row-1
				output_matrix[row,col] = output_matrix[col,row]
			end
		end

		@debug "Returned symmetric matrix" output_matrix
		return output_matrix, changed_indices2
	else
		return output_matrix, changed_indices
	end
end

function scramble_matrix(in_matrix::Array; k::Int=2, max_iterations=-1)
	out_matrix = copy(in_matrix)
	total_rows, total_cols = size(in_matrix)
	counter = 0
	if max_iterations < 1
		max_iterations = (total_cols*(total_cols-1))/2
	end

	for row = 1:k:total_rows-k
		# @info "row:" row
		for col=total_cols:-k:row+1
			# @info "col:" col
			if row == col-k+1
				# @info "shoulbreak"
				continue
			end
			indices = collect(CartesianIndices((row:row+k-1,col-k+1:col)))
			permut_indices = shuffle(indices)
			out_matrix[indices] .= in_matrix[permut_indices]
			counter +=1
			if counter >= max_iterations
				break
			end
		end
		if counter >= max_iterations
			break
		end
	end
	return out_matrix
end

