using Eirene
# using Plots
# using PlotThemes
# using PlotUtils
# using Dierckx
# using StatsPlots
#%%
# ==============================
#  ======== Tested code ========


#%%
# ================================
#  ======== Untested code ========

function get_barcodes(results_eirene::Dict, max_dim::Integer; min_dim::Int = 1, sorted::Bool=false)
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
    barcodes = Matrix{Float64}[]
    for d = min_dim:max_dim
        result = barcode(results_eirene, dim = d)
        if isempty(result) && d > 1
            result = zeros(size(barcodes[d-1]))
        end
        if sorted
            result = result[sortperm(result[:,1]),:]
        end
        push!(barcodes, result)
    end
    return barcodes
end


function plot_barcodes(barcodes::Vector;
                        min_dim::Integer = 1,
                        barcodes_labels::Bool = true,
                        default_labels::Bool = true,
                        kwargs...)#; plot_size = (width=1200, height=800),
    """
    	plot_barcodes(barcodes; min_dim::Integer=1, betti_labels::Bool=true, default_labels::Bool=true kwargs...)

    Creates a plot for set of barcodesstored in `barcodes` and return the
    handler to the plot.

    'kwargs' are plot parameters

    Some of the possible 'kwargs' are:
    	- title::String
    	- legend:Bool
    	- size::Tuple{T, T} where {T::Number}
    	- lw::Integer or linewidth:Integer
    (for more, see plots documentation):
    """
    # TODO min_dim is not included in all_dims variable
    # TODO add change of x label based on x values- so it is either edge density for 0:1 range values or Filtration step otherwise
    # TODO add ordering of bars to firstly birth time, then death time

    # barcodes = all_barcodes_geom
    max_dim = size(barcodes, 1) - (1-min_dim) # TODO not sure if this is correct solution
    dims_indices = 1:length(min_dim:max_dim)
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
        lw = 8
    end

    colors_set = get_bettis_color_palete(min_dim=min_dim)
    plot_ref = plot(; kwargs...)

    all_sizes = [size(barcodes[k],1) for k = dims_indices]
    ranges_sums = vcat(0,[sum(all_sizes[1:k]) for k = dims_indices])
    y_val_ranges = [ranges_sums[k]+1:ranges_sums[k+1] for k = dims_indices]

    # for p = min_dim:(max_dim) #TODO ths can not be starting from min_dim, because it may be 0
    for (p, dim) = enumerate(min_dim:max_dim)# 1:(max_dim) #TODO ths can not be starting from min_dim, because it may be 0
        # @info p, dim
        args = (lc = colors_set[p], linewidth = lw)

        b = barcodes[p][1:1:(end-1),:]

        if dim==0
            b = sort(b, dims=1)
        end



        total_bars = size(b,1)
        y_vals = [[k, k] for k in y_val_ranges[p]]
        lc = colors_set[p]
        for k = 1:total_bars
            # TODO change label to empty one
            plot!(b[k,:], y_vals[k]; label="", lc=lc)#; args...)
        end
        if false && betti_labels
            label = "β$(dim)"
        end
        # plot!(label=label)
    end

    # display(plot_ref)

    legend_pos = findfirst(x -> x == :legend, keys(kwargs))
    if !isnothing(legend_pos)
        plot!(legend = kwargs[legend_pos])
    else
        all_labels = reshape(["β$(k)" for k in 1:max_dim], (1,max_dim))
        plot!(label = all_labels)
        plot!(legend = true)
    end

    x_pos = findfirst(x -> x == :xlabel, keys(kwargs))
    y_pos = findfirst(x -> x == :ylabel, keys(kwargs))
    if !isnothing(x_pos)
        xlabel!(kwargs[x_pos])
    elseif default_labels
        xlabel!("Birth/Death")
    end
    if !isnothing(y_pos)
        ylabel!(kwargs[y_pos])
    elseif default_labels
        ylabel!("Cycle")
    end

    return plot_ref
end


# TODO This has to be imported from other file
# function get_bettis_color_palete(; min_dim = 1, use_set::Integer = 1)
#     """
#     	function get_bettis_color_palete()
#
#     Generates vector with colours used for Betti plots. Designed for Betti plots consistency.
#     """
#     # TODO what does the number in the function below is used for?
#
#     if use_set == 1
#         cur_colors = [Gray(bw) for bw = 0.0:0.025:0.5]
#         if min_dim == 0
#             colors_set = [RGB(87 / 256, 158 / 256, 0 / 256)]
#         else
#             colors_set = RGB[]
#         end
#         colors_set = vcat(
#             colors_set,
#             [
#                 RGB(255 / 256, 206 / 256, 0 / 256),
#                 RGB(248 / 256, 23 / 256, 0 / 256),
#                 RGB(97 / 256, 169 / 256, 255 / 256),
#                 RGB(163 / 256, 0 / 256, 185 / 256),
#                 RGB(33 / 256, 96 / 256, 45 / 256),
#                 RGB(4 / 256, 0 / 256, 199 / 256),
#                 RGB(135 / 256, 88 / 256, 0 / 256),
#             ],
#             cur_colors,
#         )
#     else
#         use_set == 2
#         cur_colors = get_color_palette(:auto, 1)
#         cur_colors3 = get_color_palette(:lightrainbow, 1)
#         cur_colors2 = get_color_palette(:cyclic1, 1)
#         if min_dim == 0
#             # colors_set =  [cur_colors[3], cur_colors[5], [:red], cur_colors[1]] #cur_colors[7],
#             colors_set = [cur_colors3[3], cur_colors[5], cur_colors3[end], cur_colors[1]] #cur_colors[7],
#         else
#             colors_set = [cur_colors[5], cur_colors3[end], cur_colors[1]] #cur_colors[7],
#             # colors_set =  [cur_colors[5], [:red], cur_colors[1], cur_colors[14]]
#         end
#         # for c =  [collect(11:25);]
#         #     push!(colors_set, cur_colors2[c])
#         # end
#         colors_set = vcat(colors_set, [cur_colors2[c] for c in [collect(11:25);]])
#     end
#
#     return colors_set
# end

function get_birth_death_ratio(barcodes; max_dim::Integer=3)
    # birth_death_raio_π = [[all_barcodes_geom[k][m,2]/all_barcodes_geom[k][m,1] for m= 1:size(all_barcodes_geom[k],1)] for k in 1:max_dim]
    birth_death_ratio_π = [barcodes[k][:,2]./barcodes[k][:,1] for k in 1:max_dim]
    return birth_death_ratio_π
end

function get_barcode_lifetime(barcodes; max_dim::Integer=3)
    lifetime = [barcodes[k][:,2].-barcodes[k][:,1] for k in 1:max_dim]
    return lifetime
end

#%%
function get_barcode_max_lifetime(lifetimes)
    """
        get_barcode_max_lifetime(lifetimes, min_dim, max_dim)

    Returns the maximal life times of barcode for all dimensions.
    """
    total_lifetimes = length(lifetimes)
    all_max_lifetimes = zeros(total_lifetimes, 1)

    for k in 1:total_lifetimes
        all_max_lifetimes[k] = findmax(lifetimes[k])[1]
    end

    return all_max_lifetimes
end


function boxplot_birth_death(birth_death_ratio_π, min_dim::Integer, max_dim::Integer)
    """
        boxplot_birth_death(areas_matrix, min_dim::Integer, max_dim::Integer)

    Plots the boxplot of area under betti curves.
    """
    bplot = StatsPlots.boxplot()

    data_colors = get_bettis_color_palete()
    total_plots = size(birth_death_ratio_π,1)

    for k in 1:total_plots
        StatsPlots.boxplot!(bplot, birth_death_ratio_π[k], labels="β$(k)", color=data_colors[k])
        StatsPlots.dotplot!(bplot, birth_death_ratio_π[k],color=data_colors[k])
    end

    return bplot
end

function boxplot_lifetime(barcode_lifetime, min_dim::Integer, max_dim::Integer)
    """
        boxplot_birth_death(areas_matrix, min_dim::Integer, max_dim::Integer)

    Plots the boxplot of area under betti curves.
    """
    bplot = StatsPlots.boxplot()

    data_colors = get_bettis_color_palete()
    total_plots = size(barcode_lifetime,1)

    for k in 1:total_plots
        StatsPlots.boxplot!(bplot, barcode_lifetime[k], labels="β$(k)", color=data_colors[k])
        StatsPlots.dotplot!(bplot, barcode_lifetime[k],color=data_colors[k])
    end

    return bplot
end

#%%
function get_barcode_max_db_ratios(db_ratos)
    """
        get_barcode_max_db_ratios(lifetimes, min_dim, max_dim)

    Returns the maximal life times of barcode for all dimensions.
    """
    total_db = length(db_ratos)
    all_max_db = zeros(total_db, 1)

    for k in 1:total_db
        all_max_db[k] = findmax(db_ratos[k])[1]
    end

    return all_max_db
end


function get_normalised_barcodes(barcodes, betti_numbers::Array)
    """
        get_normalised_barcodes(barcodes::Vector, betti_numbers::Array)

    Returns the barcodes which values are within [0,1] range.
    1. get corresponding bettis
    2. get max number of steps
    3. divide all values by total number of steps.

   """
    if typeof(betti_numbers) == Vector
        total_steps = size(betti_numbers[1],1)
    else
        total_steps = size(betti_numbers,1)
    end

    return barcodes ./ total_steps
end


function get_normalised_barcodes_collection(barcodes_collection, bettis_collection)
    """
        get_normalised_barcodes_collection(barcodes_collection, bettis_collection)

    Applies get_normalised_barcodes to the collection of barcodes and corresponding
    betti curves.
    """
    if size(barcodes_collection,1) != size(bettis_collection,1)
        throw(BoundsError(barcodes_collection, bettis_collection,
            "Both collections must have same number of elements",
        ))
    else
        total_collections  = size(barcodes_collection,1)
    end
    normed_collection = deepcopy(barcodes_collection)

    for k = 1:total_collections
        normed_collection[k]=get_normalised_barcodes(barcodes_collection[k], bettis_collection[k])
    end
    return normed_collection
end


function plot_bd_diagram(barcodes; dims=1:length(barcodes), use_js::Bool=false,
                                                class_sizes=[],
                                                class_labels=[],
                                                kwargs...)
    """
    	plot_bd_diagram(barcodes;
                            dims::Range,
                            use_js::Bool=false,
                            kwargs...)

    Creates a birth/death diagram from `barcodes` and returns the handlers to
    the plots.

    By default, dims is set to range '1:length(barcodes)', which plots all of
    the diagrams. If set to an integer, plots only 1 dimension.

    If 'use_js' is set to true, plotly backend is used for plotting.

    'kwargs' are plot parameters

    Some of the possible 'kwargs' are:
    	- title::String
    	- legend:Bool
    	- size::Tuple{T, T} where {T::Number}
    	- lw::Integer or linewidth:Integer
    (for more, see plots documentation):
    TODO dims are defined as if the dimensions always starts at 1- this has to be changed
    """
    # TODO max min should be ready to use from input data- might be better to have better structures as an inupt
    max_dim = size(barcodes, 1)
    min_dim = findmin(dims)[1]
    all_dims = min_dim:max_dim

    if findmax(dims)[1] > max_dim
        throw(DomainError(
            min_dim,
            "\'dims\' must be less than maximal dimension in \'bettis\'",
        ))
    end

    lw_pos = findfirst(x -> x == :lw || x == :linewidth, keys(kwargs))
    if !isnothing(lw_pos)
        lw = kwargs[lw_pos]
    else
        lw = 2
    end

    colors_set = TopologyPreprocessing.get_bettis_color_palete(min_dim=1)

    if use_js
        plotly()
    else
        gr()
    end

    plot_ref = plot(;xlims=(0,1), ylims=(0,1), kwargs...)

    for p in dims
        # colors_set[p]
        my_vec = barcodes[p]

        # TODO class size is not a default and basic bd diagram property- should be factored out to other function
        if class_labels != [] && class_sizes != []
            labels = ["class/size $(class_labels[k])/$(class_sizes[class_labels[k]])" for k in 1:size(class_labels,1)]
        elseif  class_sizes == []
            labels = ["class: $(k)" for k in 1:size(my_vec,1)]
        else
            labels = ["class/size $(k)/$(class_sizes[k])" for k in 1:size(my_vec,1)]
        end

        args = (color = colors_set[p],
                linewidth = lw,
                label="β$(p)",
                aspect_ratio=1,
                size = (600,600),
                legend = :bottomright,
                hover = labels,
                kwargs...)
        if class_labels != [] class_sizes != []
            for x = class_labels
                plot!(my_vec[Int(x), 1], my_vec[Int(x), 2], seriestype = :scatter; args...)
            end
        end

        plot!(my_vec[:, 1], my_vec[:, 2], seriestype = :scatter; args...)
    end

    # Add diagonal
    all_births = vcat([barcodes[d][:,1] for d in dims]...)
    all_deaths = vcat([barcodes[d][:,2] for d in dims]...)
    max_x = findmax(all_births)[1]
    max_y = findmax(all_deaths)[1]
    plot!([0, max_y], [0, max_y], label = "")

    return plot_ref
end

#%%
function plot_all_bd_diagrams(barcodes_collection;
                                min_dim::Integer = 1,
                                betti_labels::Bool = true,
                                default_labels::Bool = true,
                                all_legend=false,
                                my_alpha=0.12,
                                aspect_ratio=1,
                                base_w=600,
                                base_h=600,
                                kwargs...)
    """
    	plot_all_bd_diagrams(barcodes_collection;
                            min_dim::Integer=1,
                            betti_labels::Bool=true,
                            default_labels::Bool=true,
                            all_legend=false,
                            my_alpha=0.12,
                            aspect_ratio=1,
                            kwargs...)

    Creates a set of birth/death diagrams from `barcodes_collection`
         and returns a dictionary with the handlers to the plots.

    'kwargs' are plot parameters

    Some of the possible 'kwargs' are:
    	- title::String
    	- legend:Bool
    	- size::Tuple{T, T} where {T::Number}
    	- lw::Integer or linewidth:Integer
    (for more, see plots documentation):
    """
    total_dims = size(barcodes_collection[1],1)

    lw_pos = findfirst(x -> x == :lw || x == :linewidth, keys(kwargs))
    if !isnothing(lw_pos)
        lw = kwargs[lw_pos]
    else
        lw = 2
    end

    title_pos = findfirst(x -> x == :title, keys(kwargs))
    if !isnothing(title_pos)
        my_title = kwargs[title_pos]
    else
        my_title = "Birth death diagram"
    end

    colors_set = TopologyPreprocessing.get_bettis_color_palete(min_dim=min_dim)
    plot_dict = Dict()

    for b = 1:total_dims
        args = (lc = colors_set[b],
                linewidth = lw,
                label=false,
                aspect_ratio=aspect_ratio,
                size = (base_w,base_h),
                kwargs...)
        plot_dict["β$(b)"] = scatter(;xlims=(0,1), ylims=(0,1), dpi=300, args...)
        for bars = barcodes_collection
            barcode = bars[b]

            scatter!(barcode[:, 1], barcode[:, 2],
                    markeralpha=my_alpha,
                    markercolor = colors_set[b],
                    dpi=300)
        end
        plot!(legend=all_legend)
        plot!(title=(my_title*", β$(b)"))

        # legend_pos = findfirst(x -> x == :legend, keys(kwargs))
        # if !isnothing(legend_pos)
        #     plot!(legend = kwargs[legend_pos])
        # else
        #     plot!(legend = betti_labels)
        # end

        x_pos = findfirst(x -> x == :xlabel, keys(kwargs))
        y_pos = findfirst(x -> x == :ylabel, keys(kwargs))
        if !isnothing(x_pos)
            xlabel!(kwargs[x_pos])
        elseif default_labels
            xlabel!("Birth")
        end
        if !isnothing(y_pos)
            ylabel!(kwargs[y_pos])
        elseif default_labels
            ylabel!("Death")
        end
    end

    return plot_dict
end


## ===-
# Simpler plotting

function plot_simple_bd_diagram(barcodes; dims=1:length(barcodes), max_bd=0, use_js::Bool=false,
                                                kwargs...)
    """
    	plot_bd_diagram(barcodes;
                            dims::Range,
                            use_js::Bool=false,
                            kwargs...)

    Creates a birth/death diagram from `barcodes` and returns the handlers to
    the plots.

    By default, dims is set to range '1:length(barcodes)', which plots all of
    the diagrams. If set to an integer, plots only 1 dimension.

    If 'use_js' is set to true, plotly backend is used for plotting.

    'kwargs' are plot parameters

    Some of the possible 'kwargs' are:
    	- title::String
    	- legend:Bool
    	- size::Tuple{T, T} where {T::Number}
    	- lw::Integer or linewidth:Integer
    (for more, see plots documentation):
    TODO dims are defined as if the dimensions always starts at 1- this has to be changed
    """
    # TODO max min should be ready to use from input data- might be better to have better structures as an inupt
    max_dim = size(barcodes, 1)
    min_dim = findmin(dims)[1]
    all_dims = min_dim:max_dim

    if findmax(dims)[1] > max_dim
        throw(DomainError(
            min_dim,
            "\'dims\' must be less than maximal dimension in \'bettis\'",
        ))
    end

    lw_pos = findfirst(x -> x == :lw || x == :linewidth, keys(kwargs))
    if !isnothing(lw_pos)
        lw = kwargs[lw_pos]
    else
        lw = 2
    end

    colors_set = TopologyPreprocessing.get_bettis_color_palete(min_dim=1)

    if use_js
        plotly()
    else
        gr()
    end

    plot_ref = plot(;kwargs...)

    for p in dims
        # colors_set[p]
        my_vec = barcodes[p]

        args = (color = colors_set[p],
                linewidth = lw,
                aspect_ratio=1,
                size = (600,600),
                legend = :bottomright,
                kwargs...)

        # scatter!(my_vec[:, 1], my_vec[:, 2], args...)
        plot!(my_vec[:, 1], my_vec[:, 2], seriestype = :scatter; args...)
    end

    # Add diagonal

    if max_bd>0
        max_x = max_bd
        max_y = max_bd
        plot!([0, max_y], [0, max_y], label = "")
    else
        all_births = vcat([barcodes[d][:,1] for d in dims]...)
        all_deaths = vcat([barcodes[d][:,2] for d in dims]...)
        max_x = findmax(all_births)[1]
        max_y = findmax(all_deaths)[1]
        plot!([0, max_y], [0, max_y], label = "")
    end

    return plot_ref
end
