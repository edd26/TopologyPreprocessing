# new_component
# using Eirene
using DelimitedFiles
 	using Plots
    using LinearAlgebra
    using Images
    using Distances
    using Images
    using JLD

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
		reorganized_plt_ref = plot_square_heatmap(unscrambled_matrix, 1,size(unscrambled_matrix,1);
									plt_title = "unscrambled_matrix, size:$(matrix_size)",
									color_palete=:lightrainbow)
		display(reorganized_plt_ref)
	end
	return unscrambled_matrix, reorganized_plt_ref
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
		reorganized_plt_ref = plot_square_heatmap(reordered_matrix, 1,size(reordered_matrix,1);
									plt_title = "reordered_matrix, size:$(matrix_size)",
									color_palete=:lightrainbow)
		display(reorganized_plt_ref)
	end
	return reordered_matrix, reorganized_plt_ref
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
		fine_tuned_plt_ref = plot_square_heatmap(fine_tune_matrix, 1,size(reordered_matrix,1);
										plt_title = "fine_tuned, size:$(matrix_size)",
										color_palete=:lightrainbow)
		display(fine_tuned_plt_ref)
	end
	return fine_tune_matrix, fine_tuned_plt_ref
end

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
	reorganized_plt_ref = plot_square_heatmap(sorted_matrix, 1,size(reordered_matrix,1);
							plt_title = "reordered_matrix, size:$(matrix_size)",
							color_palete=:lightrainbow)

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

	# for every row in matrix
	for k = 1:2:matrix_size-1
		# global reordered_matrix
		reorganized_plt_ref_pt0 = plot_square_heatmap(reordered_matrix, 1,size(reordered_matrix,1);
								plt_title = "reordered_matrix, size:$(matrix_size)",
								color_palete=:lightrainbow)
		max_val, max_ind = findmax(reordered_matrix[k:end, k:end])
		# Take the smaller coordinate
		# (max_ind[1] < max_ind[2]) ? (target_row = max_ind[1]) : (target_row = max_ind[2])
		target_row = max_ind[1]+k-1

		reordered_matrix = swap_rows(reordered_matrix, k, target_row)
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
			reorganized_plt_ref_pt2 = plot_square_heatmap(reordered_matrix, 1,size(reordered_matrix,1);
										plt_title = "reordered_matrix, size:$(matrix_size)",
										color_palete=:lightrainbow)
			reorganized_plt_ref = plot(reorganized_plt_ref_pt0, reorganized_plt_ref_pt1, reorganized_plt_ref_pt2, layout=(1,3), size=(1400,400))
			display(reorganized_plt_ref)
		end
	end
	if do_final_plot
		reorganized_plt_ref = plot_square_heatmap(reordered_matrix, 1,size(reordered_matrix,1);
									plt_title = "reordered_matrix, size:$(matrix_size)",
									color_palete=:lightrainbow)
		# display(reorganized_plt_ref)
	else
		reorganized_plt_ref=[]
	end
	return reordered_matrix, reorganized_plt_ref
end


##
"""
   matrix_poling!(input_matrix; method = "max_pooling")

Takes a matrix and changes it's values to the same value, according to 'method'.
Possible methods are:
- 'max_pooling'- finds maximal value and replaces all values with the maximal
	value.
"""
# function matrix_poling!(input_matrix::Array; method::String = "max_pooling")
# 	if method == "max_pooling"
# 		max_val = findmax(input_matrix)[1]
# 		input_matrix .= max_val
# 	end
# 	return input_matrix
# end

function matrix_poling(input_matrix::Array; method = "max_pooling")
	out_matrix = copy(input_matrix)
	if method == "max_pooling"
		max_val = findmax(out_matrix)[1]
		out_matrix .= max_val
	end
	return out_matrix
end


subsamp_size=2
square_matrix = copy(sqr_matrix0)
[r_beg:r_end,c_beg:c_end])

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
		return
	end
end

function reorganize_matrix(square_matrix::Array; subsamp_size::Int=2, method="max_pooling")
	# Subsample upper half
	square_matrix2 = copy(square_matrix)
	total_rows, total_cols = size(square_matrix)

	if total_rows%2 != 0
		total_rows -= 1
	end
	if total_cols%2 != 0
		total_cols -= 1
	end


	for row = 1:subsamp_size:(total_rows-2)
		for col = (row+subsamp_size):subsamp_size:total_cols
			r_beg = row
			r_end = row+subsamp_size-1
			c_beg = col
			c_end = col+subsamp_size-1

			square_matrix2[r_beg:r_end,c_beg:c_end] = matrix_poling(square_matrix2[r_beg:r_end,c_beg:c_end])
		end
	end

	# Copy over lower half
	for row in 2:total_rows
		for col in 1:row-1
			square_matrix2[row,col] = square_matrix2[col,row]
		end
	end
	return square_matrix2
end

function pool_matrix(square_matrix::Array; method="max_pooling")
	out_matrix = copy(square_matrix)
	pool_matrix!(square_matrix; method=method)
	return out_matrix
end
