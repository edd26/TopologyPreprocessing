#=
Creted: 2020-05-04
Author: Emil Dmitruk

Structures for storing matrices used at different stages of preprocessing for
topological analysis with Eirene library.
=#

# TODO draw a diagram of all structures
# ===
struct MethodsParams
	plot_filters::Bool
	plot_heatmaps::Bool
	plot_betti_figrues::Bool
	lower_dist_mat_resolution::Bool
	lower_ord_mat_resolution::Bool
	legend_on::Bool
	do_dsiplay::Bool
	save_gabor_params::Bool
	save_subelements::Bool
	save_figure::Bool

	function MethodsParams(;plot_filters = false,
								plot_heatmaps = true,
								plot_betti_figrues = true,
								lower_dist_mat_resolution = false,
								lower_ord_mat_resolution = false,
								legend_on = true,
								do_dsiplay = false,
								save_gabor_params = false,
								save_subelements = false,
								save_figure = false)

		new(plot_filters, plot_heatmaps, plot_betti_figrues,
			lower_dist_mat_resolution, lower_ord_mat_resolution,
			legend_on, do_dsiplay, save_gabor_params, save_subelements,
			save_figure)
	end
end

struct ImageTopologyParams
	total_bins::Int
	max_size_limiter::Int
	min_B_dim::Int
	max_B_dim::Int
	file_name::String
	img_path::String
	gruping::Bool
	sub_img_size::Int
	sub_sample_size::Int
	pooling_method::String
	gabor_set::Int
	overlap::Number
	gaussian_blurr::Number

	function ImageTopologyParams(;total_bins = 5, max_size_limiter = 200,
								min_B_dim = 1, max_B_dim = 4,
								file_name = "", img_path="img/",
	    						gruping = true, sub_img_size = 33, sub_sample_size=2,
								pooling_method = "avg_pooling", gabor_set = 4, overlap = 0.0,
								gaussian_blurr=0.25)

		new(total_bins, max_size_limiter, min_B_dim, max_B_dim,
			file_name, img_path, gruping, sub_img_size, sub_sample_size,
			pooling_method, gabor_set, overlap, gaussian_blurr)
	end
end

function get_params_description(par::ImageTopologyParams)
	if par.pooling_method == "gauss_pooling"
		met = "_$(par.pooling_method)$(ceil(Int, par.gaussian_blurr*100))"
	else
		met = "_$(par.pooling_method)"
	end

	return "file_$(split(par.file_name,'.')[1])"*
            "_subimgsize$(par.sub_img_size)"*
            "_maxB_$(par.max_B_dim)"*
			"_minB_$(par.min_B_dim)"*
			"_gaborset_$(par.gabor_set)"*
			"_poolingsize_$(par.sub_sample_size)"*
            "$(met)_"*
			"_overlap_$(Int(par.overlap*10))_"
end

function get_ord_mat_from_img(par::ImageTopologyParams, met_par::MethodsParams; get_distances=false)
	@info "Current img_size" par.sub_img_size
	@info "Using file: " par.file_name
	@debug "Used params: " par.total_bins, par.gabor_set
	size_limiter = par.max_size_limiter

	# =============================== Get image masks ======================
	masks = get_filter_set(par.gabor_set, par.sub_img_size; plot_filters = false,
										save_gabor_params = met_par.save_gabor_params)

	# =============================== Get image ================================
	file_n = split(par.file_name, ".")[1]
	loaded_img = load(par.img_path*par.file_name)
	img1_gray = Gray.(loaded_img)
	img_size = size(img1_gray)

	# ================================ Process Image =======================
	# get the correlation matrix accordingly to choosen method
	local_correlations = get_local_correlations("gabor", img1_gray,img_size,
										par.sub_img_size; masks = masks,
										overlap=par.overlap)

	# ======================== Compute pairwise correlation ================
	dist_mat = pairwise(Euclidean(), local_correlations, dims=2)

	# =============================== Ordered matrix =======================
	size_limiter = size(dist_mat,1)
	if size_limiter > par.max_size_limiter
		@warn "Restricting matrix size, because matrix is too big"
		size_limiter = par.max_size_limiter
	end

	# ===
	# Matrix gupping
	met_par.lower_dist_mat_resolution && group_distances!(dist_mat, par.total_bins)

	# ===
	# Matrix ordering
	ordered_matrix = get_ordered_matrix(dist_mat[1:size_limiter,1:size_limiter];
												assign_same_values=par.gruping)
	if met_par.lower_ord_mat_resolution
		ordered_matrix = lower_ordmat_resolution(ordered_matrix, par.total_bins)
	end
	if get_distances
		return dist_mat
	else
		return ordered_matrix
	end
end

# ===
struct TopologyMatrixSet
    file_name::String

	# 2 below are not necessary, if params are includede in this structure
	sub_sample_size::Int
	pooling_method::String

	# Matrices
	ordered_matrix::Array
	reordered_matrix
		# reordered_map_ref
	pooled_matrix
	pooled_reord_matrix
	renum_pooled_orig_matrix
	renum_pooled_reord_matrix
	reordered_renum_pooled_orig_matrix

	matrix_collection
	description_vector
	ranks_collection

	params::ImageTopologyParams

	function TopologyMatrixSet(input_matrix::Array, params::ImageTopologyParams
								; description_vector)
		# TODO add parameter which describes which methods should be used
		# @warn "Using constant in structure definition- TODO: change to variable"
	    # file_name = images_set[6]
		# sub_sample_size = 2
		# pooling_method = "avg_pooling"
		# ===
		file_name = params.file_name

		reordered_matrix, reordered_map_ref =
	 			order_max_vals_near_diagonal2(input_matrix; do_final_plot=false, do_all_plots = false);

		pooled_matrix   	= reorganize_matrix(input_matrix; subsamp_size=params.sub_sample_size, method=params.pooling_method, gauss_sigma=params.gaussian_blurr)
		pooled_reord_matrix = reorganize_matrix(reordered_matrix; subsamp_size=params.sub_sample_size, method=params.pooling_method, gauss_sigma=params.gaussian_blurr)
		# =
		# gaussian_blurr = g_blurr
		# used_kernel = Kernel.gaussian(gaussian_blurr)
		# pooled_matrix   	= ceil.(Int,imfilter(input_matrix, used_kernel))
		# pooled_reord_matrix = ceil.(Int,imfilter(reordered_matrix, used_kernel))
		# =

		renum_pooled_orig_matrix  = get_ordered_matrix(pooled_matrix; assign_same_values=true)
		renum_pooled_reord_matrix = get_ordered_matrix(pooled_reord_matrix; assign_same_values=true)

		reordered_renum_pooled_orig_matrix, reordered_renum_pooled_orig_matrix_ref =
	 			order_max_vals_near_diagonal2(renum_pooled_orig_matrix; do_final_plot=false, do_all_plots = false);
		matrix_collection = Array[]
		push!(matrix_collection, input_matrix)
		push!(matrix_collection, reordered_matrix)
		push!(matrix_collection, pooled_matrix)
		push!(matrix_collection, pooled_reord_matrix)
		push!(matrix_collection, renum_pooled_orig_matrix)
		push!(matrix_collection, renum_pooled_reord_matrix)
		push!(matrix_collection, reordered_renum_pooled_orig_matrix)
		ranks_collection = zeros(Int,size(matrix_collection)[1])
		for mat = 1: size(matrix_collection)[1]
			ranks_collection[mat] = rank(matrix_collection[mat])
		end

		new(file_name, params.sub_sample_size,  params.pooling_method, input_matrix,
			reordered_matrix,
			# reordered_map_ref,
			pooled_matrix, pooled_reord_matrix,
			renum_pooled_orig_matrix, renum_pooled_reord_matrix,
			reordered_renum_pooled_orig_matrix,
			matrix_collection, description_vector, ranks_collection, params)
	end
end

# ===
struct TopologyMatrixBettisSet
	min_B_dim::Int
	max_B_dim::Int
	bettis_collection

	function TopologyMatrixBettisSet(top_mat_set::TopologyMatrixSet;min_B_dim=1, max_B_dim=3)
		bettis_collection = Any[]
		for matrix = top_mat_set.matrix_collection
 		   # ===
 		   # Persistent homology
 		   eirene_geom = eirene(matrix,maxdim=top_mat_set.params.max_B_dim,model="vr")
 		   bett_geom = get_bettis(eirene_geom, top_mat_set.params.max_B_dim, min_dim = top_mat_set.params.min_B_dim)
		   push!(bettis_collection,bett_geom)
	   end
	   new(min_B_dim, max_B_dim, bettis_collection)
   end
end

# ===
struct MatrixHeatmap
	heat_map
	matrix_property::String

	function MatrixHeatmap(in_array, description)
		hmap_len = size(in_array)[1]
		ordered_map1 = plot_square_heatmap(in_array, 5, hmap_len; plt_title=description,)
		new(ordered_map1, description)
	end
end

# ===
struct TopologyMatrixHeatmapsSet
	heatmaps::Array{MatrixHeatmap}

	heatmap_plots_set

	function TopologyMatrixHeatmapsSet(topology_matrix::TopologyMatrixSet)
		heatmaps = [MatrixHeatmap(topology_matrix.ordered_matrix,"original"),
					MatrixHeatmap(topology_matrix.reordered_matrix,"reordered"),
					MatrixHeatmap(topology_matrix.pooled_matrix,"pooled_origi"),
					MatrixHeatmap(topology_matrix.pooled_reord_matrix,"pooled_reordered"),
					MatrixHeatmap(topology_matrix.renum_pooled_orig_matrix,"renum_pooled_orig"),
					MatrixHeatmap(topology_matrix.renum_pooled_reord_matrix,"renum_pooled_reord"),
					MatrixHeatmap(topology_matrix.reordered_renum_pooled_orig_matrix,"reordered_renum_pooled_original"),
					]
		heatmap_plots_set = Any[]
		for hmap in heatmaps
			push!(heatmap_plots_set,hmap.heat_map)
		end
		new(heatmaps,heatmap_plots_set)
	end
end

# ===
struct BettiPlot
	betti_plot

	function BettiPlot(in_array; min_B_dim=1)
		betti_plot = plot_bettis2(in_array, "", legend_on=false, min_dim=min_B_dim);
		xlabel!("Steps");

		new(betti_plot)
	end
end

# ===
struct TopologyMatrixBettisPlots
	betti_plots_set
	function TopologyMatrixBettisPlots(bettis_collection::TopologyMatrixBettisSet)
		total_bettis = size(bettis_collection.bettis_collection)[1]
		betti_plots_set = Any[]
		for bett = 1:total_bettis
			betti_plot = BettiPlot(bettis_collection.bettis_collection[bett])
			push!(betti_plots_set, betti_plot.betti_plot)
		end

		new(betti_plots_set)
	end
end
