# ==============================
#  ======== Tested code ========
using Eirene
using Plots
using PlotThemes
using PlotUtils
using StatsPlots
# using Dierckx

#%%
function get_bettis(results_eirene::Dict, max_dim::Integer; min_dim::Int = 1)
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
function normalise_bettis(bettis::Vector)
    """
    	normalise_bettis(bettis::Vector)
        normalise_bettis(bettis::Array)

    Normalise the number of steps for Betti curves. 'bettis' can be either vector of
    arrays (each array contain Betti curve of different dimension) or an array
    containing Betti curve of a single dimension.
    """

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
function vectorize_bettis(betti_curves::Vector{Array{Float64,2}})
    """
        vectorize_bettis(betti_curves::Matrix{Float64})

    Reshapes the 'betti_curves' from type Array{Matrices{Float64,2}} into
    Matrix{Float64}.

    The resulting matrix size is (n, k), where 'n' is equal to the number of
    rows in each matrix, 'k' is equal to the number of matrices.

    TODO: Change the name- it takse vector and returns a matrix.
    TODO: get bettis could have an arguent betti_type which would determine resulting type
    """

    first_betti = 1
    last_betti = size(betti_curves,1)
    return hcat([betti_curves[k][:, 2] for k = first_betti:last_betti]...)
end

#%%
@deprecate vectorize_bettis(eirene_results::Dict, maxdim::Integer, mindim::Integer) vectorize_bettis(betti_curves)

# ===
#%%
function get_vectorized_bettis(results_eirene::Dict, max_dim::Integer; min_dim::Int = 1)
    """
    	get_vectorized_bettis(results_eirene::Dict, max_dim::Integer; min_dim::Int=1)

    Takes the eirene result and computes Betti curves for dimensions in range
    'mindim:maxdim'. Every Betti curve is stored in successive column of the
    resulting array.
    TODO: this should return a matrix, where first col are indices and rest are B values (1st col is missing now)
    """

    all_bettis = get_bettis(results_eirene, max_dim, min_dim = min_dim)
    bettis_vector = vectorize_bettis(all_bettis)

    return bettis_vector
end

# ==
#%%
function plot_bettis(bettis::Vector;
                        min_dim::Integer = 1,
                        use_edge_density::Bool=true,
                        betti_labels::Bool = true,
                        default_labels::Bool = true,
                        kwargs...)#; plot_size = (width=1200, height=800),
    """
    	plot_bettis(bettis; min_dim::Integer=1, betti_labels::Bool=true, default_labels::Bool=true kwargs...)

    Creates a plot for set of betti numbers stored in `bettis` and return the
    handler to the plot.

    'kwargs' are plot parameters

    Some of the possible 'kwargs' are:
    	- title::String
    	- legend:Bool
    	- size::Tuple{T, T} where {T::Number}
    	- lw::Integer or linewidth:Integer
    (for more, see plots documentation):
    TODO: min_dim is not included in all_dims variable
    TODO: add change of x label based on x values- so it is either edge density for 0:1 range values or Filtration step otherwise
    """
    max_dim = size(bettis, 1)
    all_dims = 1:max_dim

    if min_dim > max_dim
        throw(DomainError(
            min_dim,
            "\'min_dim\' must be greater that maximal dimension in \'bettis\'",
        ))
    end

    lw_pos = findfirst(x -> x == :lw || x == :linewidth, keys(kwargs))
    if !isnothing(lw_pos)
        lw = kwargs[lw_pos]
    else
        lw = 2
    end

    # Create iterator for all loops
    all_iterations = 1:(max_dim) #TODO ths can not be starting from min_dim, because it may be 0

    if use_edge_density
        for p = all_iterations
            max_step = findmax(bettis[p][:, 1])[1]
            bettis[p][:, 1] ./= max_step
        end
    end

    colors_set = get_bettis_color_palete(min_dim=min_dim)
    plot_ref = plot(; kwargs...)
    # for p = min_dim:(max_dim) #TODO ths can not be starting from min_dim, because it may be 0
    # for p = all_iterations
    for (index, p) in enumerate(min_dim:max_dim)
        args = (lc = colors_set[index], linewidth = lw)
        if betti_labels
            args = (args..., label = "β$(p)")
        end
        plot!(bettis[index][:, 1], bettis[index][:, 2]; args...)
    end

    legend_pos = findfirst(x -> x == :legend, keys(kwargs))
    if !isnothing(legend_pos)
        plot!(legend = kwargs[legend_pos])
    else
        plot!(legend = betti_labels)
    end

    x_pos = findfirst(x -> x == :xlabel, keys(kwargs))
    y_pos = findfirst(x -> x == :ylabel, keys(kwargs))
    if !isnothing(x_pos)
        xlabel!(kwargs[x_pos])
    elseif default_labels
        xlabel!("Edge density")
    end
    if !isnothing(y_pos)
        ylabel!(kwargs[y_pos])
    elseif default_labels
        ylabel!("Number of cycles")
    end

    # set tlims to integer values
    max_ylim = findmax(ceil.(Int, ylims(plot_ref)))[1]
    if max_ylim <=3
        ylims!((0, 3))
    end

    if use_edge_density
        xlims!((0, 1))
    end

    return plot_ref
end

function plot_bettis(bettis::Array;
                        min_dim::Integer = 1,
                        use_edge_density::Bool=true,
                        betti_labels::Bool = true,
                        default_labels::Bool = true,
                        normalised=true,
                        kwargs...)#; plot_size = (width=1200, height=800),
    """
    	plot_bettis(bettis::Array; min_dim::Integer=1, betti_labels::Bool=true, default_labels::Bool=true kwargs...)

    Creates a plot for set of betti numbers stored in `bettis` and return the
    handler to the plot.

    'kwargs' are plot parameters

    Some of the possible 'kwargs' are:
    	- title::String
    	- legend:Bool
    	- size::Tuple{T, T} where {T::Number}
    	- lw::Integer or linewidth:Integer
    (for more, see plots documentation):
    TODO: min_dim is not included in all_dims variable
    TODO: add change of x label based on x values- so it is either edge density for 0:1 range values or Filtration step otherwise
    """
    max_dim = size(bettis, 2)-1-min_dim
    all_dims = 1:max_dim

    if min_dim > max_dim
        throw(DomainError(
            min_dim,
            "\'min_dim\' must be greater that maximal dimension in \'bettis\'",
        ))
    end

    total_steps = size(bettis, 1)
    if normalised
        x_vals = range(0, stop=1, length=total_steps)
    else
        x_vals = range(0, stop=total_steps)
    end

    lw_pos = findfirst(x -> x == :lw || x == :linewidth, keys(kwargs))
    if !isnothing(lw_pos)
        lw = kwargs[lw_pos]
    else
        lw = 2
    end

    # if use_edge_density
    #     # for p = 1:(max_dim) #TODO ths can not be starting from min_dim, because it may be 0
    #     for (index, p) in enumerate(min_dim:max_dim)
    #         max_step = findmax(bettis[:, 1])[1]
    #         bettis[p][:, 1] ./=max_step
    #     end
    # end

    colors_set = get_bettis_color_palete(min_dim=min_dim)
    plot_ref = plot(; kwargs...)
    # for p = min_dim:(max_dim) #TODO ths can not be starting from min_dim, because it may be 0
    # for p = 1:(max_dim) #TODO ths can not be starting from min_dim, because it may be 0
    for (index, p) in enumerate(min_dim:max_dim)
        args = (lc = colors_set[index], linewidth = lw)
        if betti_labels
            args = (args..., label = "β$(p)")
        end
        plot!(x_vals, bettis[:, index]; args...)
    end

    legend_pos = findfirst(x -> x == :legend, keys(kwargs))
    if !isnothing(legend_pos)
        plot!(legend = kwargs[legend_pos])
    else
        plot!(legend = betti_labels)
    end

    x_pos = findfirst(x -> x == :xlabel, keys(kwargs))
    y_pos = findfirst(x -> x == :ylabel, keys(kwargs))
    if !isnothing(x_pos)
        xlabel!(kwargs[x_pos])
    elseif default_labels
        xlabel!("Edge density")
    end
    if !isnothing(y_pos)
        ylabel!(kwargs[y_pos])
    elseif default_labels
        ylabel!("Number of cycles")
    end

    # set tlims to integer values
    max_ylim = findmax(ceil.(Int, ylims(plot_ref)))[1]
    if max_ylim <=3
        ylims!((0, 3))
    end

    if use_edge_density
        xlims!((0, 1))
    end

    return plot_ref
end
#  ======= Untested code
# TODO add default kwargs paring function -> parse_kwargs()

function plot_all_bettis(bettis_collection;
                        min_dim::Integer = 1,
                        betti_labels::Bool = true,
                        default_labels::Bool = true,
                        normalised=true,
                        kwargs...)#; plot_size = (width=1200, height=800),
    """
    	plot_all_bettis ...
    """
    total_dims = size(bettis_collection[1],2)

    lw_pos = findfirst(x -> x == :lw || x == :linewidth, keys(kwargs))
    if !isnothing(lw_pos)
        lw = kwargs[lw_pos]
    else
        lw = 2
    end

    colors_set = get_bettis_color_palete(min_dim=min_dim)
    max_y_val = find_max_betti(bettis_collection)

    plot_ref = plot(; kwargs...)
    for b = 1:total_dims
        args = (lc = colors_set[b], linewidth = lw, alpha=0.12,label=false, ylims=(0,max_y_val))
        for bettis = bettis_collection
            betti_vals = bettis[:,b]

            total_steps = size(bettis, 1)
            x_vals = range(0, stop=1, length=total_steps)

            plot!(x_vals, betti_vals; args...)
        end
        # my_label = "β$(b)"
        # betti_vals = results_d["bettis_collection"][:hc][end]
        # x_vals = range(0, stop=1, length=size(betti_vals, 1))
        # plot!(x_vals, betti_vals; lc = colors_set[b], linewidth = 1, alpha=0.1,label=my_label, ylims=(0,max_y_val))
    end
    plot!(legend=true)

    legend_pos = findfirst(x -> x == :legend, keys(kwargs))
    if !isnothing(legend_pos)
        plot!(legend = kwargs[legend_pos])
    else
        plot!(legend = betti_labels)
    end

    x_pos = findfirst(x -> x == :xlabel, keys(kwargs))
    y_pos = findfirst(x -> x == :ylabel, keys(kwargs))
    if !isnothing(x_pos)
        xlabel!(kwargs[x_pos])
    elseif default_labels
        xlabel!("Edge density")
    end
    if !isnothing(y_pos)
        ylabel!(kwargs[y_pos])
    elseif default_labels
        ylabel!("Number of cycles")
    end

    return plot_ref
end


function find_max_betti(bettis_collection::Array)
    """
        find_max_betti(bettis_collection::Array)

    Returns the highest Betti curve value from all dimensions.
    """
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

#  ======= Untested code == end
#%%
function printready_plot_bettis(kwargs)
    """
    	printready_plot_bettis(kwargs)

    Creates a plot using 'plot_bettis' with arguments which were tested to be very
    good for using them in prints. Used arguments are:

    """
    return nothing
end


#%%
function get_bettis_color_palete(; min_dim = 1, use_set::Integer = 1)
    """
    	function get_bettis_color_palete()

    Generates vector with colours used for Betti plots. Designed for Betti plots consistency.
    """
    # TODO what does the number in the function below is used for?

    if use_set == 1
        cur_colors = [Gray(bw) for bw = 0.0:0.025:0.5]
        if min_dim == 0
            colors_set = [RGB(87 / 256, 158 / 256, 0 / 256)]
        else
            colors_set = []
        end
        max_RGB = 256
        colors_set = vcat(
            colors_set,
            [
                RGB(255 / max_RGB, 206 / max_RGB, 0 / max_RGB),
                RGB(248 / max_RGB, 23 / max_RGB, 0 / max_RGB),
                RGB(97 / max_RGB, 169 / max_RGB, 255 / max_RGB),
                RGB(163 / max_RGB, 0 / max_RGB, 185 / max_RGB),
                RGB(33 / max_RGB, 96 / max_RGB, 45 / max_RGB),
                RGB(4 / max_RGB, 0 / max_RGB, 199 / max_RGB),
                RGB(135 / max_RGB, 88 / max_RGB, 0 / max_RGB),
            ],
            cur_colors,
        )
    else
        use_set == 2
        cur_colors = get_color_palette(:auto, 1)
        cur_colors3 = get_color_palette(:lightrainbow, 1)
        cur_colors2 = get_color_palette(:cyclic1, 1)
        if min_dim == 0
            # colors_set =  [cur_colors[3], cur_colors[5], [:red], cur_colors[1]] #cur_colors[7],
            colors_set = [cur_colors3[3], cur_colors[5], cur_colors3[end], cur_colors[1]] #cur_colors[7],
        else
            colors_set = [cur_colors[5], cur_colors3[end], cur_colors[1]] #cur_colors[7],
            # colors_set =  [cur_colors[5], [:red], cur_colors[1], cur_colors[14]]
        end
        # for c =  [collect(11:25);]
        #     push!(colors_set, cur_colors2[c])
        # end
        colors_set = vcat(colors_set, [cur_colors2[c] for c in [collect(11:25);]])
    end

    return colors_set
end

# ==============================
#  ======= Untested code =======

# using Measures
# using Plots.PlotMeasures

#
# # Source: https://github.com/JuliaPlots/Plots.jl/issues/897
# function setdefaultplottingparams(;upscale=2)
#     #8x upscaling in resolution
#     fntsm = Plots.font("sans-serif", pointsize=round(12.0*upscale))
#     fntlg = Plots.font("sans-serif", pointsize=round(18.0*upscale))
#     default(titlefont=fntlg, guidefont=fntlg, tickfont=fntsm, legendfont=fntsm)
#     default(size=(800*upscale,600*upscale)) #Plot canvas size
#     default(dpi=500) #Only for PyPlot - presently broken
# end





#%%
function plot_bettis_collection(bettis_collection,
                                bett_num,
                                max_rank;
                                step = 1,
                                show_plt = true,
                                R = 0.0,
                                G = 0.4,
                                B = 1.0)
    """
    	plot_bettis_collection(bettis_collection, bett_num; step=1, show_plt=true, R=0., G=0.4, B=1.0)

    PLots collection of Betti curves of rank 'bett-num'. Every successive plot has
    lower opacity than predecessor.  'step' defines step between collection elements
    that are ploted. By default, plot is displayed after carteation. This can be
    disabled by setting 'show_plt' to false.

    Color of the plot can be set with 'R', 'G', 'B' parameters.
    """
    step > 0 || error("Steps should be natural number!")
    bettis_total = size(bettis_collection, 1)
    colors_set = zeros(Float64, bettis_total, 4)
    colors_set[:, 1] .= R
    colors_set[:, 2] .= G
    colors_set[:, 3] .= B
    max_betti = get_max_betti_from_collection(bettis_collection)
    @info "max_betti" max_betti

    x = 0
    y = bettis_total * 0.1
    va_range = collect(range(bettis_total + x, y, length = bettis_total))

    colors_set[:, 4] .= va_range / findmax(va_range)[1]
    rgba_set = RGBA[]
    for k = 1:size(colors_set, 1)
        push!(
            rgba_set,
            RGBA(colors_set[k, 1], colors_set[k, 2], colors_set[k, 3], colors_set[k, 4]),
        )
    end

    plt_reference = plot(1, title = "Betti curves collection, rank $(bett_num)", label = "")
    for b = 1:step:bettis_total
        betti = bettis_collection[b]
        x_vals_1 = (1:size(betti[:, bett_num], 1)) / size(betti[:, bett_num], 1)
        plot!(x_vals_1, betti[:, bett_num], lc = rgba_set[b], label = "rank=$(max_rank-b)")
        plot!(ylim = (0, max_betti))
    end
    xlabel!("Normalised steps")
    ylabel!("Number of cycles")
    plot!(legend = true)

    show_plt && display(plt_reference)
    return plt_reference
end

#%%
function get_max_bettis(bettis)
    """
        get_max_bettis(bettis)

    Returns the maximal bettis of Betti curves for all dimensions.
    """

    all_max_bettis = findmax(bettis, dims=1)[1]

    return all_max_bettis
end

# TODO change name
# TODO check what for dim is used, change to min dim
function get_max_betti_from_collection(bettis_collection; dim = 1)
    max_betti = 0
    for betti in bettis_collection
        # global max_betti
        local_max = findmax(betti)[1]
        if (local_max > max_betti)
            max_betti = local_max
        end
    end
    return max_betti
end


#%%
function plot_and_save_bettis(bettis,
                                plot_title::String,
                                results_path::String;
                                file_name = "",
                                extension = ".png",
                                do_save = true,
                                do_normalise = true,
                                min_dim = 0,
                                max_dim = 3,
                                legend_on = true,
                                kwargs...)
    """
    	plot_and_save_bettis(eirene_results, plot_title::String,
    								results_path::String; extension = ".png",
    								data_size::String="", do_save=true,
    								extend_title=true, do_normalise=true, max_dim=3,
    								legend_on=true)

    Plot Betti curves from 0 up to `max_dim` using `eirene_results` from Eirene library and
    returns handler for figure. Optionally, if `do_save` is set, saves the figure
    or if `do_normalise` is set, sets the steps range to be normalised to the
    horizontal axis maximal value.
    """
    bettis = get_bettis(eirene_results, max_dim)
    if do_normalise
        bettis = normalise_bettis(bettis)
    end
    plot_ref =
        plot_bettis(bettis, plot_title, legend_on = legend_on, min_dim = min_dim, kwargs...)


    if do_save
        if isempty(file_name)
            file_name = plot_title * extension
        elseif isempty(findall(x -> x == extension[2:end], split(file_name, ".")))
            #check for the extension in file name
            file_name *= extension
        end

        save_figure_with_params(
            plot_ref,
            results_path;
            extension = extension,
            prefix = split(file_name, ".")[1],
        )

    end
    return plot_ref
end


# TODO merge functions for getting betti curves
# Original function returns 2 different types of betti curves. If no default
# value parameters is given, it returns vector of matrices. If num of steps is
# given, then it return matrix maxdim x numsteps.
# """
# bettis_eirene(matr, maxdim; mintime=-Inf, maxtime=Inf, numofsteps=Inf, mindim=1)
#
# Takes the `matr` and computes Betti curves up to `maxdim`. Return matrix only
# with betti curve values
#
#
# Function taken from: https://github.com/alexyarosh/hyperbolic
# """
#%%

@deprecate bettis_eirene(matr, maxdim; mintime = -Inf, maxtime = Inf, numofsteps = Inf, mindim = 1) get_bettis(results_eirene, max_dim; min_dim = 1)

#%%
function get_bettis_from_image(img_name,
                                plot_params;
                                file_path = "",
                                plot_heatmaps = true,
                                save_heatmaps = false,
                                plot_betti_figrues = true)
    """
    	function get_bettis_from_image(img_name)

    Computes Betti curves for the image file indicated by @img_name. If the image is
    	not symmetric, then it is the elements below diagonal are copied over the
    	elmenents above the diagonal.
    """
    file_n = split(img_name, ".")[1]
    img1_gray = Gray.(load(file_path * img_name))
    img_size = size(img1_gray)

    C_ij = Float64.(img1_gray)

    if !issymmetric(C_ij)
        img1_gray = symmetrize_image(img1_gray)
        C_ij = Float64.(img1_gray)
    end
    img_size = size(C_ij, 1)
    # C_ij =-C_ij
    # C_ij .+= 1


    # ==============================================================================
    # =============================== Ordered matrix ===============================
    if size(C_ij, 1) > 80
        @warn "Running Eirene for big matrix: " img_size
        @warn "Eirene may have trobules with big matrices/images."
    end

    ordered_matrix = get_ordered_matrix(C_ij; assing_same_values = false)


    # ==============================================================================
    # ============================ Persistance homology ============================
    C = eirene(ordered_matrix, maxdim = 3, model = "vr")


    # ==============================================================================
    # ================================ Plot results ================================

# TODO separate plotting from processing
    if plot_heatmaps

        full_ordered_matrix = get_ordered_matrix(C_ij; assing_same_values = false)
        heat_map2 = plot_square_heatmap(
            full_ordered_matrix,
            10,
            img_size;
            plt_title = "Order matrix of $(file_n)",
            plot_params = plot_params,
        )

        if save_heatmaps
            heatm_details = "_heatmap_$(file_n)"
            savefig(heat_map2, heatmaps_path * "ordering" * heatm_details)
        end
    end

    if plot_betti_figrues
        plot_title = "Betti curves of $(file_n), size=$(img_size) "
        figure_name = "betti_$(file_n)_n$(img_size)"
        ref = plot_and_save_bettis(C,
                                    plot_title,
                                    figure_path,;
                                    file_name = figure_name,
                                    plot_params = plot_params,
                                    do_save = false,
                                    extend_title = false,
                                    do_normalise = false,
                                    max_dim = 3,
                                    legend_on = true,
                                    min_dim = 1)
    end
    display(img1_gray)
    display(heat_map2)
    display(ref)
end

# ===============================================
@deprecate get_bettis_from_image2(img_name;file_path = "",plot_heatmaps = true, save_heatmaps = false, plot_betti_figrues = true) get_bettis_from_image(img_name, plot_params; file_path = "", plot_heatmaps = true, save_heatmaps = false, plot_betti_figrues = true)

@deprecate plot_and_save_bettis2(eirene_results, plot_title::String, results_path::String; file_name = "", extension = ".png", data_size::String = "", do_save = true, extend_title = true, do_normalise = true, min_dim = 0, max_dim = 3, legend_on = true) plot_and_save_bettis(bettis, plot_title::String, results_path::String; file_name = "", extension = ".png", do_save = true, do_normalise = true, min_dim = 0, max_dim = 3, legend_on = true, kwargs...)

#%%
function get_and_plot_bettis(eirene_results;
                                max_dim = 3,
                                min_dim = 1,
                                plot_title = "",
                                legend_on = false)
    bettis = get_bettis(eirene_results, max_dim)
    norm_bettis = normalise_bettis(bettis)
    plot_ref =
        plot_bettis2(norm_bettis, plot_title, legend_on = legend_on, min_dim = min_dim)
    # display(plot_ref)

    return plot_ref
end


#%%
function lower_ordmat_resolution(ordered_matrix::Array, total_bins::Int)
    """
    	lower_ordmat_resolution(ordered_matrix::Array, total_bins::Int)

    Takes ordered matrix 'input_matrix' and reduces the resolution of values in the
    matrix into 'total_bins' bins.
    """
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
function average_bettis(bettis_matrix::Matrix; up_factor = 8)
    """
    	average_bettis(bettis_matrix; up_factor=8)

    Takes the average values of betti curves stored in 'bettis_matrix'.

    'bettis_matrix' consist of different simulations(first index of the matrix),
    different ranks (third index of the matrix). Second index of the matrices
    (saples) may vary accross many simulations and for this reason, all betti curves
    are upsampled by a factor of 'upsample_factor' and then the average for every
    dimension is computed.
    """

    bettis_matrix_backup = copy(bettis_matrix)

    simulations = size(bettis_matrix, 1)
    dimensions = size(bettis_matrix[1], 1)

    max_samples = 0
    for k = 1:simulations
        # global max_samples
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

    # ref = plot(range(0, 1, length=total_orig_samples), input_vector);
    # plot!(x_vals, y_upsampled);
    # display(ref)

    return y_upsampled
end




#%%
function upsample_vector(input_vector; upsample_factor::Int = 8)
    """
    	upsample_vector(input_vector; upsample_factor::Int=8)

    Takes an 'input_vector' and returns a vector which has 'upsample_factor' many
    times more samples. New samples are interpolated with 'spl' function from
    'Dierckx' package.

    """
    total_orig_samples = size(input_vector, 1) - 1
    total_samples = upsample_factor * total_orig_samples + 1

    x_vals = range(0, 1, length = total_orig_samples + 1)
    spl = Spline1D(x_vals, input_vector)

    x_upsampled = range(0, 1, length = total_samples)
    y_upsampled = spl(x_upsampled)

    # ref = plot(range(0, 1, length=total_orig_samples), input_vector);
    # plot!(x_vals, y_upsampled);
    # display(ref)

    return y_upsampled
end


# =========--=======-========-==========-=======-
# From bettis areas
# Area under Betti curve functions
#%%
function get_area_under_betti_curve(betti_curves::Union{Matrix{Float64}, Array{Array{Float64,2}}};do_normalised::Bool=false)
    """
        get_area_under_betti_curve(betti_curves, min_dim, max_dim)

    Computes the area under Betti curves stored in 'betti_curves', where each row is
    a Betti curve and each column is a value.
    """
    #TODO check this part
    if size(betti_curves,2) < 2
        bettis_vector = vectorize_bettis(betti_curves)
    else
        bettis_vector = betti_curves
    end
    # @info sum(bettis_vector, dims=1)
    bettis_area = sum(bettis_vector, dims=1)

    if do_normalised
        total_steps = size(bettis_vector,1)
        bettis_area ./= total_steps
    end
    # @info bettis_area
    return bettis_area
end

# function get_area_under_betti_curve(C, min_dim, max_dim)
#     """
#         get_area_under_betti_curve(C, min_dim, max_dim)
#
#     Computes the Betti curves and returns their area under curve.
#     """
#     all_bettis = get_bettis(C,max_dim, min_dim=min_dim)
#     bettis_vector = hcat([all_bettis[k][:,2] for k=min_dim:max_dim]...)
#     # @info sum(bettis_vector, dims=1)
#
#
#     total_steps = size(bettis_vector,1)
#
#     bettis_area = sum(bettis_vector, dims=1) ./ total_steps
#     # @info bettis_area
#     return bettis_area
# end




#%%
function get_dataset_bettis_areas(dataset; min_dim::Integer=1, max_dim::Integer=3, return_matrix::Bool=true)
    """
        get_dataset_bettis_areas(dataset; min_dim::Integer=1, max_dim::Integer=3, return_matrix::Bool=true)

    Computes topology of every matrix in dataset, computes Betti curves for dimensions
    min_dim up to max_dim and returns vector (or matrix) of areas under Betti curves.
    """
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

# struct TopologyData
#     min_dim::Integer
#     max_dim::Integer
#
#     do_normalise::Bool=true
#
#     betti_curves
#     normed_bettis
#     betti_areas::Matrix{Int}
#
#     # Constructor for input data
#     function TopologyData(my_matrix::Matrix, max_dim::Int; min_dim::Int, do_normalise::Bool=true)
#         min_dim = min_dim
#         max_dim = max_dim
#
#         @info "Computing topology for maxdim =" max_dim
#         C = eirene(my_matrix, maxdim=max_dim)
#         betti_curves = get_bettis(C, max_dim, min_dim=min_dim)
#         normed_bettis = normalise_bettis(betti_curves)
#         betti_areas = get_area_under_betti_curve(betti_curves; do_normalised=do_normalise)
#     end
# end

#%%
function get_dataset_topology(dataset;
                              min_dim::Integer=1,
                              max_dim::Integer=3,
                              get_curves::Bool=true,
                              get_areas::Bool=true,
                              get_persistence_diagrams::Bool=true,
                              do_normalise::Bool=true)
    topology_set = TopologyData[]
    for some_matrix in dataset
        resulting_topology = TopologyData(some_matrix, max_dim, min_dim=min_dim, do_normalise=do_normalise)
        push!(topology_set, resulting_topology)
    end
    return topology_set
end


#%%
function get_area_boxes(areas_matrix, min_dim::Integer, max_dim::Integer)
    """
        get_area_boxes(areas_matrix, min_dim::Integer, max_dim::Integer)

    Plots the boxplot of area under betti curves.
    """
    bplot = StatsPlots.boxplot()

    data_colors = get_bettis_color_palete()

    for (index, value) in enumerate(min_dim:max_dim)
        StatsPlots.boxplot!(bplot, areas_matrix[:,index], labels="β$(value)", color=data_colors[value])
    end

    return bplot
end

function get_bettis_collection_from_matrices(ordered_matrices_collection; max_dim::Int=3, min_dim::Int=1)
    bettis_collection = Array[]

    for matrix = ordered_matrices_collection
		@debug "Computing Bettis..."
		eirene_geom = eirene(matrix,maxdim=max_B_dim,model="vr")

		bettis = reshape_bettis(get_bettis(eirene_geom, max_B_dim))
		push!(bettis_collection, bettis)
    end

    return bettis_collection
end

# =========--=======-========-==========-=======-
# Code from Points substitution:

# Compute series of betti curves
# function get_bettis_collection(ordered_matrices_collection; max_B_dim=3)
#     bettis_collection = Array[]
#
#     for matrix = ordered_matrices_collection
# 		@debug "Computing Bettis..."
# 		eirene_geom = eirene(matrix,maxdim=max_B_dim,model="vr")
#
# 		bettis = reshape_bettis(get_bettis(eirene_geom, max_B_dim))
# 		push!(bettis_collection, bettis)
#     end
#
#     return bettis_collection
# end
#
# # Plot series of betti curves with their heatmaps
# function reshape_bettis(bettis)
# 	bettis_count = size(bettis,1)
# 	output_betti = zeros(size(bettis[1],1), bettis_count)
#
# 	for betti = 1:bettis_count
# 		output_betti[:,betti] = bettis[betti][:,2]
# 	end
# 	return output_betti
# end
#
# function get_ord_mat_collection(matrix_collection)
# 	mat_size = size(matrix_collection[1],1)
# 	ordered_mat_coll = [zeros(Int, mat_size,mat_size) for k=1:length(matrix_collection)]
#
# 	size(matrix_collection)
# 	for matrix = 1:length(matrix_collection)
# 		ordered_mat_coll[matrix] = Int.(get_ordered_matrix(matrix_collection[matrix]))
# 	end
# 	return ordered_mat_coll
# end
#
#
#
#
#
# function print_hmap_with_bettis(ordered_matrices_collection, bettis_collection,
# 														plot_data::PlottingData)
# 	num_plots = size(ordered_matrices_collection,1)
# 	sources = 1:(plot_data.src_pts_number)
# 	targets = 1:(plot_data.trgt_pts_number)
# 	plot_set = Any[]
#
#     max_betti = get_max_betti_from_collection(bettis_collection;dim=1)
#
# 	index = 1
# 	for src = 1:size(sources,1), trgt = 1:size(targets,1)
#         # index = src * trgt
#         ordered_geom_gr = ordered_matrices_collection[index]
#         bettis = bettis_collection[index]
#         title_hmap = "trgt:$(targets[trgt])_src:$(sources[src])_r:$(rank(ordered_geom_gr))"
#         title_bettis = "gr_trg=$(targets[trgt])_src=$(sources[src])_steps=$(size(bettis,1))"
#         push!(plot_set, make_hm_and_betti_plot(ordered_geom_gr, bettis, title_hmap, title_bettis, max_betti))
# 		index +=1
# 	end
#
# 	return plot_set
# end
#
# function make_hm_and_betti_plot(ordered_geom_gr, bettis, title_hmap, title_bettis, max_betti)
#     # @debug "src" src
#     # @debug "trgt" trgt
#     hmap_plot = plot_square_heatmap(ordered_geom_gr, 10,size(ordered_geom_gr,1);plt_title = title_hmap)
#     plot!(yflip = true,)
#
#     bettis_plot_ref = plot(title=title_bettis);
#     max_dim = size(bettis,2)
#     for p = 1:max_dim
#         x_vals = collect(1:size(bettis[:,1],1))./size(bettis[:,1])
#
#         plot!(x_vals, bettis[:,p], label="\\beta_"*string(p));
#         plot!(legend=true, )
#     end
#
#     plot!(ylim=(0,max_betti))
# 	plot!(xlim=(0,1))
#     ylabel!("Number of cycles")
#     xlabel!("Steps")
#
#     final_plot = plot(hmap_plot, bettis_plot_ref, layout = 2,
# 						top_margin=2mm,
# 						left_margin=0mm,
# 						bottom_margin=2mm,
# 						size=(600,300))
#     display(final_plot)
#     return final_plot
# end
#
# # TODO BUG: substitution does not work- all the plots are the same
# function main_generation1()
#     mat_size = 6
#     dim = 80
#     src_pts_number = 1
#     trgt_pts_number = 2
#     trgt_steps = 0
#
#     src_points, trgt_points =
#     	get_replacing_points(mat_size, src_pts_number, trgt_pts_number)
#
#     matrix_collection =
#     	get_matrix_collection(mat_size, dim, src_points, trgt_points; trgt_step=trgt_steps)
#
#     ordered_matrices_collection = get_ord_mat_collection(matrix_collection)
#
#     bettis_collection = get_bettis_collection(ordered_matrices_collection)
#
#
#     plot_data = PlottingData(mat_size, dim, src_pts_number, trgt_pts_number, src_points, trgt_points, trgt_steps)
#     # plot_data = PlottingData2(mat_size , dim, )
#
#     plotting_data = print_hmap_with_bettis(ordered_matrices_collection,
# 													bettis_collection, plot_data)
# end
#
#
# function get_geom_matrix(mat_size, dim)
# 	# TODO change the matrix collection shape to be a matrix, not a vector
#     point_cloud = generate_random_point_cloud(mat_size, dim)
#     matrix_collection = generate_geometric_matrix(point_cloud)
#     # matrix_collection = get_ordered_matrix(matrix_collection; assing_same_values=true)
#
#     return matrix_collection
# end
#
# function get_rand_matrix(mat_size, dim)
#     matrix_collection = generate_random_matrix(mat_size)
#     matrix_collection = get_ordered_matrix(matrix_collection; assing_same_values=true)
#
#     return matrix_collection
# end
#
# # TODO Analyse zero point behaviour
# function get_dist_mat_collection(dist_matrix, src_points, trgt_points, trgt_steps; do_ordering=false)
#     dist_matrix_backup = copy(dist_matrix)
#     mat_size = size(dist_matrix,1)
#     src_points_num = size(src_points,1)
#     trgt_points_num = size(trgt_points,1)
#     # ordered_mat_coll = [zeros(Int, mat_size,mat_size) for k=1:(src_points_num*trgt_points_num)]
#     ordered_mat_coll = Array[]
#
# 	swapping_iterator = 0
#
#     for srcs = 1:src_points_num
#         # replacement_row = get_row(dist_matrix, src_points[srcs])
#
#         for target = 1:trgt_points_num
#             @debug "src:" src_points[srcs]
#             @debug "trgt:" trgt_points[target, srcs]
#             replacement_row = get_row(dist_matrix_backup, src_points[srcs])
#             # dist_matrix_backup .=
# 			set_row!(dist_matrix_backup, trgt_points[target, srcs], replacement_row)
#             # ordered_mat_coll[srcs * target] = copy(dist_matrix_backup)
# 			if do_ordering
# 				swap_rows!(dist_matrix_backup, trgt_points[target, srcs], mat_size-swapping_iterator)
# 				swapping_iterator +=1
# 			end
#             push!(ordered_mat_coll, copy(dist_matrix_backup))
#         end
#     end
#
#     return ordered_mat_coll
# end
#
# function get_ordered_set(distance_matrices_collection)
# 	result = copy(distance_matrices_collection)
#
# 	for matrix = 1:size(distance_matrices_collection,1)
# 		result[matrix] = get_ordered_matrix(distance_matrices_collection[matrix];assing_same_values=true )
# 	end
# 	return result
# end
#
# function matrix_analysis(test_data::PlottingData;generation_function=get_geom_matrix)
# 	mat_size = test_data.mat_size
# 	dim = test_data.dim
# 	src_pts_number = test_data.src_pts_number
# 	trgt_pts_number = test_data.trgt_pts_number
# 	trgt_steps = 0
#
# 	src_points, trgt_points = get_replacing_points(mat_size, src_pts_number, trgt_pts_number)
# 	distance_matrix = generation_function(mat_size, dim)
#
# 	distance_matrices_collection = get_dist_mat_collection(distance_matrix, src_points, trgt_points, trgt_steps)
# 	ordered_matrices_collection = get_ordered_set(distance_matrices_collection)
# 	bettis_collection = get_bettis_collection(ordered_matrices_collection)
#
# 	plot_data = PlottingData(mat_size, dim, src_pts_number, trgt_pts_number, src_points, trgt_points, trgt_steps)
#
# 	plots_set = print_hmap_with_bettis(ordered_matrices_collection,
# 												bettis_collection, plot_data)
#
#
# 	return distance_matrices_collection, ordered_matrices_collection, bettis_collection, plot_data, plots_set
# end

#%%
# This does not belong here
function multiscale_matrix_testing(sample_space_dims = 3,
                                    maxsim = 5,
                                    min_B_dim = 1,
                                    max_B_dim = 3,
                                    size_start = 10,
                                    size_step = 5,
                                    size_stop = 50;
                                    do_random = true,
                                    control_saving = false,
                                    perform_eavl = false)
    """
    multiscale_matrix_testing(sample_space_dims = 3,
    									maxsim=5,
    									min_B_dim = 1,
    									max_B_dim = 3,
    									size_start = 10,
    									size_step = 5,
    									size_stop = 80; do_random=true)

    Function for testing the average number of cycles from geometric and random
    	matrices.

    It is possible to save intermidiate results- for that, @control_saving must be
    set true.

    Performance of computation of Betti curves can be monitored, if the
    @perform_eavl is set too true. Bydefault, it is set to false.
    """
    num_of_bettis = length(collect(min_B_dim:max_B_dim))

    if length(sample_space_dims) > 1
        @warn "Can not do random processing for multiple dimensions"
        do_random = false
    end

    geom_mat_results = Any[]
    if do_random
        rand_mat_results = Any[]
        result_list = [geom_mat_results, rand_mat_results]
    else
        result_list = [geom_mat_results]
    end

    for sample_space_dim in sample_space_dims
        if !do_random
            @info "Sampling space size: " sample_space_dim
        end

        repetitions = size_start:size_step:size_stop
        for space_samples in repetitions
            @info "Generating data for: " space_samples
            # ==========================================
            # ============= Generate data ==============
            # ===
            # Generate random matrix
            if do_random
                symm_mat_rand = [generate_random_matrix(space_samples) for i = 1:maxsim]
                ordered_mat_rand = [
                    get_ordered_matrix(symm_mat_rand[i]; assing_same_values = false)
                    for i = 1:maxsim
                ]
            end

            # ===
            # Generate geometric matrix
            pts_rand = [
                generate_random_point_cloud(sample_space_dim, space_samples)
                for i = 1:maxsim
            ]
            symm_mat_geom = [generate_geometric_matrix(pts_rand[i]') for i = 1:maxsim]
            ordered_mat_geom = [
                get_ordered_matrix(symm_mat_geom[i]; assign_same_values = false)
                for i = 1:maxsim
            ]

            # ======================================================================
            # ========================= Do the Betti analysis ======================
            if do_random
                set = [ordered_mat_geom, ordered_mat_rand]
            else
                set = [ordered_mat_geom]
            end
            for matrix_set in set
                @debug("Betti analysis!")
                # ===
                # Generate bettis
                many_bettis = Array[]
                if perform_eavl
                    many_timings = Float64[]
                    many_bytes = Float64[]
                    many_gctime = Float64[]
                    many_memallocs = Base.GC_Diff[]
                end

                for i = 1:maxsim
                    if i % 10 == 0
                        @info "Computing Bettis for: " i
                    end

                    if perform_eavl
                        results, timing, bytes, gctime, memallocs = @timed bettis_eirene(
                            matrix_set[i],
                            max_B_dim,
                            mindim = min_B_dim,
                        )
                        push!(many_bettis, results)
                        push!(many_timings, timing)
                        push!(many_bytes, bytes)
                        push!(many_gctime, gctime)
                        push!(many_memallocs, memallocs)
                    else
                        push!(
                            many_bettis,
                            bettis_eirene(matrix_set[i], max_B_dim, mindim = min_B_dim),
                        )
                    end
                end

                # ===
                # Get maximal number of cycles from each Betti from simulations
                max_cycles = zeros(maxsim, max_B_dim)
                for i = 1:maxsim, betti_dim = 1:max_B_dim
                    @debug("\tFindmax in bettis")
                    max_cycles[i, betti_dim] = findmax(many_bettis[i][:, betti_dim])[1]
                end

                # ===
                # Get the statistics
                avg_cycles = zeros(1, length(min_B_dim:max_B_dim))
                std_cycles = zeros(1, length(min_B_dim:max_B_dim))
                k = 1
                for betti_dim = min_B_dim:max_B_dim
                    avg_cycles[k] = mean(max_cycles[:, betti_dim])
                    std_cycles[k] = std(max_cycles[:, betti_dim])
                    k += 1
                end

                # ===
                # Put results into dictionary
                betti_statistics = Dict()
                if matrix_set == ordered_mat_geom
                    @debug("Saving ordered")
                    betti_statistics["matrix_type"] = "ordered"
                    betti_statistics["space_dim"] = sample_space_dim
                    result_list = geom_mat_results
                else
                    @debug("Saving radom")
                    betti_statistics["matrix_type"] = "random"
                    result_list = rand_mat_results
                end
                betti_statistics["space_samples"] = space_samples
                betti_statistics["simualtions"] = maxsim
                betti_statistics["min_betti_dim"] = min_B_dim
                betti_statistics["max_betti_dim"] = max_B_dim
                betti_statistics["avg_cycles"] = avg_cycles
                betti_statistics["std_cycles"] = std_cycles

                if perform_eavl
                    betti_statistics["many_timings"] = many_timings
                    betti_statistics["many_bytes"] = many_bytes
                    betti_statistics["many_gctime"] = many_gctime
                    betti_statistics["many_memallocs"] = many_memallocs
                end
                push!(result_list, betti_statistics)
            end # matrix type loop
            @debug("===============")
            if control_saving
                if do_random
                    save(
                        "multiscale_matrix_testing_$(space_samples)_$(sample_space_dim).jld",
                        "rand_mat_results",
                        rand_mat_results,
                        "geom_mat_results",
                        geom_mat_results,
                    )
                else
                    save(
                        "multiscale_matrix_testing_dimension_$(space_samples)_$(sample_space_dim).jld",
                        "geom_mat_results",
                        geom_mat_results,
                    )
                end
            end
        end # matrix_size_loop
    end # sampled space dimension

    if do_random
        return geom_mat_results, rand_mat_results
    else
        return geom_mat_results
    end
end

# function plot_betti_numbers(betti_numbers, edge_density, title="Geometric  matrix"; stop=0.6)
#     """
#     Plots Betti curves. The betti numbers should be obtained with the clique-top
#     library.
#     """
#     p1 = plot(edge_density, betti_numbers[:,1], label="beta_0", title=title, legend=:topleft) #, ylims = (0,maxy)
#     plot!(edge_density, betti_numbers[:,2], label="beta_1")
#     if size(betti_numbers,2)>2
#         plot!(edge_density, betti_numbers[:,3], label="beta_2")
#     end
#
#     return p1
# end
