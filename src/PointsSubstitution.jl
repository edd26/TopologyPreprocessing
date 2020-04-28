# Set of functions
#
# After matrices generation, Betti curves will be generated
# For better optimization, operations on indices should be done first


# TODO Check if the threshold of values is applied here and if it has some
# 		consequence on results
using Eirene
	using Random


	include("MatrixProcessing.jl")
	include("MatrixToolbox.jl")
	include("BettiCurves.jl")
	include("PlottingWrappers.jl")


struct PlottingData
	mat_size::Int64
	dim::Int
	src_pts_number::Int
	trgt_pts_number::Int
	src_points::Vector{Int}
	trgt_points::Array{Int}
	targets::Vector{Int}

	# Constructor for input data
	function PlottingData(mat_size::Int, dim::Int, src_pts_number::Int,
							trgt_pts_number::Int)
		steps_set = [0]
		src_points = [0]
		trgt_points = [0]

		new(mat_size::Int, dim::Int, src_pts_number::Int,
								trgt_pts_number::Int, src_points,
								trgt_points, steps_set)
	end

	function PlottingData(mat_size::Int, dim::Int, src_pts_number::Int,
							trgt_pts_number::Int, src_points::Vector{Int},
							trgt_points::Array{Int}, trgt_sptep::Int)
		if trgt_sptep == 0
			steps_set = [trgt_pts_number]
		else
	        steps_set = collect(1:trgt_sptep:trgt_pts_number)
	        # Make sure that all points are used
	        isempty(findall(x->x==trgt_pts_number, steps_set)) && push!(steps_set, trgt_pts_number)
		end

		new(mat_size::Int, dim::Int, src_pts_number::Int,
								trgt_pts_number::Int, src_points::Vector{Int},
								trgt_points::Array{Int}, steps_set)
	end
end


function get_replacing_points(mat_size, src_pts_number, trgt_pts_number)
	total_pts_number = src_pts_number + src_pts_number*trgt_pts_number
	total_pts_number > mat_size && error("Too many points to substitute!")

	elements_collection = randcycle(mat_size)
    src_points = elements_collection[1:src_pts_number]
    strat_val = src_pts_number+1
    stop_val = total_pts_number
    trgt_points = elements_collection[strat_val:stop_val]
    trgt_points = reshape(trgt_points, trgt_pts_number, src_pts_number)

    return src_points, trgt_points
end

# Compute series of betti curves
function get_bettis_collection(ordered_matrices_collection; max_B_dim=3)
    bettis_collection = Array[]

    for matrix = ordered_matrices_collection
		@debug "Computing Bettis..."
		eirene_geom = eirene(matrix,maxdim=max_B_dim,model="vr")

		bettis = reshape_bettis(get_bettis(eirene_geom, max_B_dim))
		push!(bettis_collection, bettis)
    end

    return bettis_collection
end

# Plot series of betti curves with their heatmaps
function reshape_bettis(bettis)
	bettis_count = size(bettis,1)
	output_betti = zeros(size(bettis[1],1), bettis_count)

	for betti = 1:bettis_count
		output_betti[:,betti] = bettis[betti][:,2]
	end
	return output_betti
end

function get_ord_mat_collection(matrix_collection)
	mat_size = size(matrix_collection[1],1)
	ordered_mat_coll = [zeros(Int, mat_size,mat_size) for k=1:length(matrix_collection)]

	size(matrix_collection)
	for matrix = 1:length(matrix_collection)
		ordered_mat_coll[matrix] = Int.(get_ordered_matrix(matrix_collection[matrix]))
	end
	return ordered_mat_coll
end





function print_hmap_with_bettis(ordered_matrices_collection, bettis_collection,
														plot_data::PlottingData)
	num_plots = size(ordered_matrices_collection,1)
	sources = 1:(plot_data.src_pts_number)
	targets = 1:(plot_data.trgt_pts_number)
	plot_set = Any[]

    max_betti = get_max_betti_from_collection(bettis_collection;dim=1)

	index = 1
	for src = 1:size(sources,1), trgt = 1:size(targets,1)
        # index = src * trgt
        ordered_geom_gr = ordered_matrices_collection[index]
        bettis = bettis_collection[index]
        title_hmap = "trgt:$(targets[trgt])_src:$(sources[src])_r:$(rank(ordered_geom_gr))"
        title_bettis = "gr_trg=$(targets[trgt])_src=$(sources[src])_steps=$(size(bettis,1))"
        push!(plot_set, make_hm_and_betti_plot(ordered_geom_gr, bettis, title_hmap, title_bettis, max_betti))
		index +=1
	end

	return plot_set
end

function make_hm_and_betti_plot(ordered_geom_gr, bettis, title_hmap, title_bettis, max_betti)
    # @debug "src" src
    # @debug "trgt" trgt
    hmap_plot = plot_square_heatmap(ordered_geom_gr, 10,size(ordered_geom_gr,1);plt_title = title_hmap)
    plot!(yflip = true,)

    bettis_plot_ref = plot(title=title_bettis);
    max_dim = size(bettis,2)
    for p = 1:max_dim
        x_vals = collect(1:size(bettis[:,1],1))./size(bettis[:,1])

        plot!(x_vals, bettis[:,p], label="\\beta_"*string(p));
        plot!(legend=true, )
    end

    plot!(ylim=(0,max_betti))
	plot!(xlim=(0,1))
    ylabel!("Number of cycles")
    xlabel!("Steps")

    final_plot = plot(hmap_plot, bettis_plot_ref, layout = 2,
						top_margin=2mm,
						left_margin=0mm,
						bottom_margin=2mm,
						size=(600,300))
    display(final_plot)
    return final_plot
end

# TODO BUG: substitution does not work- all the plots are the same
function main_generation1()
    mat_size = 6
    dim = 80
    src_pts_number = 1
    trgt_pts_number = 2
    trgt_steps = 0

    src_points, trgt_points =
    	get_replacing_points(mat_size, src_pts_number, trgt_pts_number)

    matrix_collection =
    	get_matrix_collection(mat_size, dim, src_points, trgt_points; trgt_step=trgt_steps)

    ordered_matrices_collection = get_ord_mat_collection(matrix_collection)

    bettis_collection = get_bettis_collection(ordered_matrices_collection)


    plot_data = PlottingData(mat_size, dim, src_pts_number, trgt_pts_number, src_points, trgt_points, trgt_steps)
    # plot_data = PlottingData2(mat_size , dim, )

    plotting_data = print_hmap_with_bettis(ordered_matrices_collection,
													bettis_collection, plot_data)
end


function replace_matrix_rows(matrix, srcs, trgts)
    replacement_row = get_row(matrix, pt_src)

    new_matrix = set_row(matrix, trgt_points, replacement_row)
end

function get_row(matrix, pt_src)
    return  matrix[pt_src,:]
end

function set_row!(matrix::Array, pt_trgt::Int, replacement_row::Array)
    @debug "set_row! func"
    replacement_row[pt_trgt] = 0

    matrix[pt_trgt, :] .= replacement_row
    matrix[:, pt_trgt] .= replacement_row

    return matrix
end

function set_row(matrix::Array, pt_trgt::Int, replacement_row::Array)
    @debug "set_row func"
    new_matrix = copy(matrix)
    return set_row!(new_matrix, pt_trgt, replacement_row)
end

function set_row(matrix::Array, pt_trgt::Array, replacement_row::Array)
    @debug "set_row func"
    new_matrix = copy(matrix)
    for point in pt_trgt
        new_matrix = set_row!(new_matrix, point, replacement_row)
    end
    return new_matrix
end

function get_geom_matrix(mat_size, dim)
	# TODO change the matrix collection shape to be a matrix, not a vector
    point_cloud = generate_random_point_cloud(mat_size, dim)
    matrix_collection = generate_geometric_matrix(point_cloud)
    # matrix_collection = get_ordered_matrix(matrix_collection; assing_same_values=true)

    return matrix_collection
end

function get_rand_matrix(mat_size, dim)
    matrix_collection = generate_random_matrix(mat_size)
    matrix_collection = get_ordered_matrix(matrix_collection; assing_same_values=true)

    return matrix_collection
end

# TODO Analyse zero point behaviour
function get_dist_mat_collection(dist_matrix, src_points, trgt_points, trgt_steps; do_ordering=false)
    dist_matrix_backup = copy(dist_matrix)
    mat_size = size(dist_matrix,1)
    src_points_num = size(src_points,1)
    trgt_points_num = size(trgt_points,1)
    # ordered_mat_coll = [zeros(Int, mat_size,mat_size) for k=1:(src_points_num*trgt_points_num)]
    ordered_mat_coll = Array[]

	swapping_iterator = 0

    for srcs = 1:src_points_num
        # replacement_row = get_row(dist_matrix, src_points[srcs])

        for target = 1:trgt_points_num
            @debug "src:" src_points[srcs]
            @debug "trgt:" trgt_points[target, srcs]
            replacement_row = get_row(dist_matrix_backup, src_points[srcs])
            # dist_matrix_backup .=
			set_row!(dist_matrix_backup, trgt_points[target, srcs], replacement_row)
            # ordered_mat_coll[srcs * target] = copy(dist_matrix_backup)
			if do_ordering
				swap_rows!(dist_matrix_backup, trgt_points[target, srcs], mat_size-swapping_iterator)
				swapping_iterator +=1
			end
            push!(ordered_mat_coll, copy(dist_matrix_backup))
        end
    end

    return ordered_mat_coll
end

function get_ordered_set(distance_matrices_collection)
	result = copy(distance_matrices_collection)

	for matrix = 1:size(distance_matrices_collection,1)
		result[matrix] = get_ordered_matrix(distance_matrices_collection[matrix];assing_same_values=true )
	end
	return result
end

function matrix_analysis(test_data::PlottingData;generation_function=get_geom_matrix)
	mat_size = test_data.mat_size
	dim = test_data.dim
	src_pts_number = test_data.src_pts_number
	trgt_pts_number = test_data.trgt_pts_number
	trgt_steps = 0

	src_points, trgt_points = get_replacing_points(mat_size, src_pts_number, trgt_pts_number)
	distance_matrix = generation_function(mat_size, dim)

	distance_matrices_collection = get_dist_mat_collection(distance_matrix, src_points, trgt_points, trgt_steps)
	ordered_matrices_collection = get_ordered_set(distance_matrices_collection)
	bettis_collection = get_bettis_collection(ordered_matrices_collection)

	plot_data = PlottingData(mat_size, dim, src_pts_number, trgt_pts_number, src_points, trgt_points, trgt_steps)

	plots_set = print_hmap_with_bettis(ordered_matrices_collection,
												bettis_collection, plot_data)


	return distance_matrices_collection, ordered_matrices_collection, bettis_collection, plot_data, plots_set
end


function matrix_organization(matr, src_points, trgt_points)
	working_matrix = copy(matr)
	mat_size = size(matr, 1)
	src_number = size(src_points,1)

	swapping_sources = Any[]
	step = 0
	matrix_indices = CartesianIndices((1:mat_size, 1:mat_size))
	matrix_indices = findall(x->x[1]<x[2], matrix_indices)
	sorted_values = matr[matrix_indices]

	ordered_indices = sort!([1:mat_size;],
                        by=i->(sorted_values[i],matrix_indices[i]))

	sort(matr[:,1])
	# matr[src_points[src_pt],1]
	for src_pt = 1:src_number
		# find all source
		target_set = findall(x-> x==matr[src_points[src_pt],1], matr[:,1])

		# swap all equivalents
		for tragt = target_set
			swap_rows!(working_matrix, tragt, mat_size-step)
			step +=1
		end
	end
	return working_matrix
end

# matrix = ordered_matrices_collection[15]
function swap_rows!(matrix, src_row_num, trgt_row_num)
	src_backup = copy(matrix[src_row_num,:])
	trgt_backup = copy(matrix[trgt_row_num,:])

	matrix[src_row_num, trgt_row_num] = matrix[trgt_row_num, trgt_row_num]
	matrix[trgt_row_num, src_row_num] = matrix[src_row_num,src_row_num]

	matrix[src_row_num,src_row_num] = trgt_backup[src_row_num]
	matrix[trgt_row_num, trgt_row_num] = src_backup[trgt_row_num]

	src_backup = copy(matrix[src_row_num,:])
	matrix[src_row_num,:] .= matrix[trgt_row_num, :]
	matrix[trgt_row_num, :] .= src_backup

	matrix[:, src_row_num] = matrix[src_row_num,:]
	matrix[:, trgt_row_num] = matrix[trgt_row_num,:]
end

function swap_rows(matrix, src_row_num, trgt_row_num)
	new_matrix = copy(matrix)
	swap_rows!(new_matrix, src_row_num, trgt_row_num)
	return new_matrix
end

function ordering_matrix_analysis(test_data::PlottingData;generation_function=get_geom_matrix)
	mat_size = test_data.mat_size
	dim = test_data.dim
	src_pts_number = test_data.src_pts_number
	trgt_pts_number = test_data.trgt_pts_number
	trgt_steps = 0

	src_points, trgt_points = get_replacing_points(mat_size, src_pts_number, trgt_pts_number)
	distance_matrix = generation_function(mat_size, dim)

	distance_matrices_collection = get_dist_mat_collection(distance_matrix, src_points, trgt_points, trgt_steps)

	ordered_matrices_collection = get_ordered_set(distance_matrices_collection)
	bettis_collection = get_bettis_collection(ordered_matrices_collection)

	plot_data = PlottingData(mat_size, dim, src_pts_number, trgt_pts_number, src_points, trgt_points, trgt_steps)

	plotting_data = print_hmap_with_bettis(ordered_matrices_collection,
												bettis_collection, plot_data)


	return distance_matrices_collection, ordered_matrices_collection, bettis_collection, plot_data
end



# =================================
# Matrix modification functions

function make_matrix_steps!(input_matrix, step_number; step_size=2 )
	# input_matrix = ord_mat
	# step_number = 13
	rows = step_number:(step_number+step_size-1)
    cols = 1:step_number-1
	min_value = findmin(input_matrix[rows,cols])[1]
	input_matrix[rows,cols] .= min_value
	input_matrix[cols,rows] .= min_value
end



"""
function add_step_to_matrix(input_matrix, last_components)

Takes a symmetric matrix 'input_matrix' and appends 2 columns and 2 rows such that
resulting geometric object structure is bigger by 1 dimension. 'last_components'
determines which etries in the matrix are used for closing high dimensional simplices.
"""
function add_step_to_matrix(input_matrix, last_components)
	matrix_size = size(input_matrix,1)
	new_matrix = zeros(Int,matrix_size +2, matrix_size+2)
	new_matrix[1:matrix_size,1:matrix_size] .= input_matrix
	min_closing_component = findmin(input_matrix[last_components])[1]
	new_row1 = range(min_closing_component, length=matrix_size)
	new_row2 = range(findmax(new_row1)[1]+1, length=matrix_size)

	last = 2
	new_matrix[matrix_size+1,1:end-last] = new_row1
	new_matrix[matrix_size+2,1:end-last] = new_row2
	new_matrix[1:end-last,matrix_size+1] = new_row1
	new_matrix[1:end-last,matrix_size+2] = new_row2

	# Adjust last components
	max_new_matrix = findmax(new_matrix)[1]
	new_matrix[last_components].=input_matrix[last_components].+(max_new_matrix-min_closing_component+1)
	new_matrix[end-1,end  ] = findmax(new_matrix)[1]+1
	new_matrix[end,  end-1] = findmax(new_matrix)[1]

	new_max_val = findmax(new_matrix)[2]
	new_component = copy(last_components)
	push!(new_component, new_max_val)
	push!(new_component, CartesianIndex(new_max_val[2], new_max_val[1]))

	return new_matrix, new_component
end
