module TopologyPreprocessingPlotsExt

using TopologyPreprocessing
using Plots
using StatsPlots

# ==============================
# Colour palette
# ==============================

"""
  function get_bettis_color_palete()

Generates vector with colours used for Betti plots. Designed for Betti plots consistency.
"""
function TopologyPreprocessing.get_bettis_color_palete(; min_dim = 1, use_set::Integer = 1)

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
            colors_set = [cur_colors3[3], cur_colors[5], cur_colors3[end], cur_colors[1]]
        else
            colors_set = [cur_colors[5], cur_colors3[end], cur_colors[1]]
        end
        colors_set = vcat(colors_set, [cur_colors2[c] for c in [collect(11:25);]])
    end

    return colors_set
end

# ==============================
# Betti curve plots
# ==============================

"""
  plot_bettis(bettis; min_dim::Integer=1, betti_labels::Bool=true, default_labels::Bool=true kwargs...)

Creates a plot for set of betti numbers stored in `bettis` and return the handler to the plot.
"""
function TopologyPreprocessing.plot_bettis(bettis::Vector;
                        min_dim::Integer = 1,
                        use_edge_density::Bool=true,
                        betti_labels::Bool = true,
                        default_labels::Bool = true,
                        kwargs...)
    max_dim = size(bettis, 1)
    all_dims = 1:max_dim

    if min_dim > max_dim
        throw(DomainError(
            min_dim,
            "\'min_dim\' must be greater that maximal dimension.",
        ))
    end

    lw_pos = findfirst(x -> x == :lw || x == :linewidth, keys(kwargs))
    if !isnothing(lw_pos)
        lw = kwargs[lw_pos]
    else
        lw = 2
    end

    all_iterations = 1:(max_dim)

    if use_edge_density
        for p = all_iterations
            max_step = findmax(bettis[p][:, 1])[1]
            bettis[p][:, 1] ./= max_step
        end
    end

    colors_set = get_bettis_color_palete(min_dim=min_dim)
    plot_ref = plot(; kwargs...)
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

    max_ylim = findmax(ceil.(Int, ylims(plot_ref)))[1]
    if max_ylim <=3
        ylims!((0, 3))
    end

    if use_edge_density
        xlims!((0, 1))
    end

    return plot_ref
end

function TopologyPreprocessing.plot_bettis(bettis::Array;
                        dims_range=1:size(bettis,2),
                        use_edge_density::Bool=true,
                        betti_labels::Bool = true,
                        default_labels::Bool = true,
                        normalised=true,
                        kwargs...)
    max_dim = dims_range[end]
    min_dim = dims_range[1]

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

    colors_set = get_bettis_color_palete(min_dim=min_dim)
    plot_ref = plot(; kwargs...)
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

    max_ylim = findmax(ceil.(Int, ylims(plot_ref)))[1]
    if max_ylim <=3
        ylims!((0, 3))
    end

    if use_edge_density
        xlims!((0, 1))
    end

    return plot_ref
end

"""
    plot_all_bettis ...
"""
function plot_all_bettis(bettis_collection;
                        min_dim::Integer = 1,
                        betti_labels::Bool = true,
                        default_labels::Bool = true,
                        normalised=true,
                        kwargs...)
    total_dims = size(bettis_collection[1],2)

    lw_pos = findfirst(x -> x == :lw || x == :linewidth, keys(kwargs))
    if !isnothing(lw_pos)
        lw = kwargs[lw_pos]
    else
        lw = 2
    end

    colors_set = get_bettis_color_palete(min_dim=min_dim)
    max_y_val = TopologyPreprocessing.find_max_betti(bettis_collection)

    plot_ref = plot(; kwargs...)
    for b = 1:total_dims
        args = (lc = colors_set[b], linewidth = lw, alpha=0.12,label=false, ylims=(0,max_y_val))
        for bettis = bettis_collection
            betti_vals = bettis[:,b]

            total_steps = size(bettis, 1)
            x_vals = range(0, stop=1, length=total_steps)

            plot!(x_vals, betti_vals; args...)
        end
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

"""
  plot_bettis_collection(bettis_collection, bett_num; step=1, show_plt=true, R=0., G=0.4, B=1.0)

PLots collection of Betti curves of rank 'bett-num'. Every successive plot has
lower opacity than predecessor.
"""
function plot_bettis_collection(bettis_collection,
                                bett_num,
                                max_rank;
                                step = 1,
                                show_plt = true,
                                R = 0.0,
                                G = 0.4,
                                B = 1.0)
    step > 0 || error("Steps should be natural number!")
    bettis_total = size(bettis_collection, 1)
    colors_set = zeros(Float64, bettis_total, 4)
    colors_set[:, 1] .= R
    colors_set[:, 2] .= G
    colors_set[:, 3] .= B
    max_betti = TopologyPreprocessing.get_max_betti_from_collection(bettis_collection)
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

"""
  printready_plot_bettis(kwargs)

Creates a plot using 'plot_bettis' with arguments which were tested to be very
good for using them in prints.
"""
function printready_plot_bettis(kwargs)
    return nothing
end

"""
  plot_and_save_bettis(bettis, plot_title, results_path; ...)

Plot Betti curves from 0 up to `max_dim` and optionally save the figure.
"""
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
    bettis = TopologyPreprocessing.get_bettis(eirene_results, max_dim)
    if do_normalise
        bettis = TopologyPreprocessing.normalise_bettis(bettis)
    end
    plot_ref =
        plot_bettis(bettis, plot_title, legend_on = legend_on, min_dim = min_dim, kwargs...)

    if do_save
        if isempty(file_name)
            file_name = plot_title * extension
        elseif isempty(findall(x -> x == extension[2:end], split(file_name, ".")))
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

"""
    get_area_boxes(areas_matrix, min_dim::Integer, max_dim::Integer)

Plots the boxplot of area under betti curves.
"""
function get_area_boxes(areas_matrix, min_dim::Integer, max_dim::Integer)
    bplot = StatsPlots.boxplot()

    data_colors = get_bettis_color_palete()

    for (index, value) in enumerate(min_dim:max_dim)
        StatsPlots.boxplot!(bplot, areas_matrix[:,index], labels="β$(value)", color=data_colors[value])
    end

    return bplot
end

# ==============================
# Barcode plots
# ==============================

function sort_barcodes!(barcodes, min_dim, max_dim)
    sorted_barcodes = copy(barcodes)
    for (dim_index, dim) = enumerate(min_dim:max_dim)
        permutation = sortperm(barcodes[dim_index][:, 1])
        sorted_barcodes[dim_index] = barcodes[dim_index][permutation, :]
    end
    barcodes = sorted_barcodes
end

"""
  plot_barcodes!(barcodes, plot_ref; min_dim::Integer=1, kwargs...)

Adds barcode lines to an existing plot.
"""
function TopologyPreprocessing.plot_barcodes!(barcodes::Vector, plot_ref;
    min_dim::Integer=1,
    barcodes_labels::Bool=true,
    default_labels::Bool=true,
    sort_by_birth::Bool=false,
    alpha::Float64=1.0,
    kwargs...)

    max_dim = size(barcodes, 1) - (1 - min_dim)

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

    dims_indices = 1:length(min_dim:max_dim)
    all_sizes = [size(barcodes[k], 1) for k = dims_indices]
    ranges_sums = vcat(0, [sum(all_sizes[1:k]) for k = dims_indices])
    y_val_ranges = [ranges_sums[k]+1:ranges_sums[k+1] for k = dims_indices]

    if sort_by_birth
        sort_barcodes!(barcodes, min_dim, max_dim)
    end

    total_cycles = sum([size(barcodes[k], 1) for (k, dim) in enumerate(min_dim:max_dim)])
    for (p, dim) = enumerate(min_dim:max_dim)
        args = (lc=colors_set[p], linewidth=lw)

        b = barcodes[p][:, :]
        all_infs = findall(x -> isinf(x), b)
        for k in all_infs
            b[k] = 1
        end

        total_bars = size(b, 1)
        y_vals = [[k, k] for k in y_val_ranges[p]]
        lc = colors_set[p]
        for k = 1:total_bars
            plot!(b[k, :], y_vals[k]; label="", lc=lc, alpha=alpha)
        end
    end

    legend_pos = findfirst(x -> x == :legend, keys(kwargs))
    if !isnothing(legend_pos)
        plot!(legend=kwargs[legend_pos])
    else
        all_labels = reshape(["β$(k)" for k in 1:max_dim], (1, max_dim))
        plot!(label=all_labels)
        plot!(legend=true)
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

    ylims!(0, total_cycles + 2)

    return plot_ref
end

function TopologyPreprocessing.plot_barcodes(barcodes::Vector;
    kwargs...)
    plot_ref = plot(; kwargs...)
    TopologyPreprocessing.plot_barcodes!(barcodes, plot_ref;
        kwargs...)
end

"""
    boxplot_birth_death(birth_death_ratio_π, min_dim::Integer, max_dim::Integer)

Plots the boxplot of birth/death ratios.
"""
function TopologyPreprocessing.boxplot_birth_death(birth_death_ratio_π, min_dim::Integer, max_dim::Integer)
    bplot = StatsPlots.boxplot()

    data_colors = get_bettis_color_palete()
    total_plots = size(birth_death_ratio_π, 1)

    for k in 1:total_plots
        StatsPlots.boxplot!(bplot, birth_death_ratio_π[k], labels="β$(k)", color=data_colors[k])
        StatsPlots.dotplot!(bplot, birth_death_ratio_π[k], color=data_colors[k])
    end

    return bplot
end

"""
    boxplot_lifetime(barcode_lifetime, min_dim::Integer, max_dim::Integer)

Plots the boxplot of barcode lifetimes.
"""
function TopologyPreprocessing.boxplot_lifetime(barcode_lifetime, min_dim::Integer, max_dim::Integer)
    bplot = StatsPlots.boxplot()

    data_colors = get_bettis_color_palete()
    total_plots = size(barcode_lifetime, 1)

    for k in 1:total_plots
        StatsPlots.boxplot!(bplot, barcode_lifetime[k], labels="β$(k)", color=data_colors[k])
        StatsPlots.dotplot!(bplot, barcode_lifetime[k], color=data_colors[k])
    end

    return bplot
end

"""
  plot_bd_diagram(barcodes; dims, use_js::Bool=false, kwargs...)

Creates a birth/death diagram from `barcodes`.
"""
function plot_bd_diagram(barcodes; dims=1:length(barcodes),
    use_js::Bool=false,
    class_sizes=[],
    class_labels=[],
    normilised_diagonal::Bool=true,
    alpha=0.4,
    kwargs...)
    min_dim = dims[1]
    max_dim = dims[end]

    if max_dim > length(barcodes)
        throw(DimensionMismatch("Can not have dims range larger than barcodes length"))
    end

    colors_set = get_bettis_color_palete(min_dim=min_dim)

    if use_js
        plotly()
    else
        gr()
    end

    plot_ref = plot(; xlims=(0, 1), ylims=(0, 1), kwargs...)
    if normilised_diagonal
        max_coord = 1
    else
        max_x = max([k for k in vcat([barcodes[d][:, 1] for (d, dim) in dims |> enumerate]...) if !isinf(k)]...)
        max_y = max([k for k in vcat([barcodes[d][:, 2] for (d, dim) in dims |> enumerate]...) if !isinf(k)]...)
        max_coord = max(max_x, max_y)
    end

    scaling_factor = 1.05
    min_val = -0.05
    plot!([0, scaling_factor * max_coord], [0, scaling_factor * max_coord], label="")
    xlims!(min_val, scaling_factor * max_coord)
    ylims!(min_val, scaling_factor * max_coord)

    for (p, dim) in enumerate(dims)
        my_vec = barcodes[p]

        if class_labels != [] && class_sizes != []
            labels = ["class/size $(class_labels[k])/$(class_sizes[class_labels[k]])" for k in 1:size(class_labels, 1)]
        elseif class_sizes == []
            labels = ["class: $(k)" for k in 1:size(my_vec, 1)]
        else
            labels = ["class/size $(k)/$(class_sizes[k])" for k in 1:size(my_vec, 1)]
        end

        args = (color=colors_set[p],
            label="β$(dim)",
            aspect_ratio=1,
            size=(600, 600),
            legend=:bottomright,
            framestyle=:origin,
            alpha=alpha,
            kwargs...)
        if class_labels != []
            class_sizes != []
            for x = class_labels
                plot!(my_vec[Int(x), 1], my_vec[Int(x), 2], seriestype=:scatter; args...)
            end
        end

        plot!(my_vec[:, 1], my_vec[:, 2], seriestype=:scatter; args...)
    end

    xlabel!("birth")
    ylabel!("death")

    return plot_ref
end

"""
  plot_all_bd_diagrams(barcodes_collection; min_dim::Integer=1, kwargs...)

Creates a set of birth/death diagrams from `barcodes_collection` and returns
a dictionary with handlers to the plots.
"""
function plot_all_bd_diagrams(barcodes_collection;
    min_dim::Integer=1,
    betti_labels::Bool=true,
    default_labels::Bool=true,
    all_legend=false,
    my_alpha=0.12,
    aspect_ratio=1,
    base_w=600,
    base_h=600,
    kwargs...)
    total_dims = size(barcodes_collection[1], 1)

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

    colors_set = get_bettis_color_palete(min_dim=min_dim)
    plot_dict = Dict()

    for b = 1:total_dims
        args = (lc=colors_set[b],
            linewidth=lw,
            label=false,
            aspect_ratio=aspect_ratio,
            size=(base_w, base_h),
            kwargs...)
        plot_dict["β$(b)"] = scatter(; xlims=(0, 1), ylims=(0, 1), dpi=300, args...)
        for bars = barcodes_collection
            barcode = bars[b]

            scatter!(barcode[:, 1], barcode[:, 2],
                markeralpha=my_alpha,
                markercolor=colors_set[b],
                dpi=300)
        end
        plot!(legend=all_legend)
        plot!(title=(my_title * ", β$(b)"))

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

"""
  plot_simple_bd_diagram(barcodes; dims, max_bd=0, use_js::Bool=false, kwargs...)

Creates a simplified birth/death diagram.
"""
function plot_simple_bd_diagram(barcodes; dims=1:length(barcodes), max_bd=0, use_js::Bool=false,
    kwargs...)
    max_dim = size(barcodes, 1)
    min_dim = findmin(dims)[1]

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

    colors_set = get_bettis_color_palete(min_dim=1)

    if use_js
        plotly()
    else
        gr()
    end

    plot_ref = plot(; kwargs...)

    for (p, dim) in dims |> enumerate
        my_vec = barcodes[p]

        args = (color=colors_set[p],
            linewidth=lw,
            aspect_ratio=1,
            size=(600, 600),
            legend=:bottomright,
            kwargs...)

        plot!(my_vec[:, 1], my_vec[:, 2], seriestype=:scatter; args...)
    end

    if max_bd > 0
        max_x = max_bd
        max_y = max_bd
        plot!([0, max_y], [0, max_y], label="")
    else
        all_births = vcat([barcodes[d][:, 1] for d in dims]...)
        all_deaths = vcat([barcodes[d][:, 2] for d in dims]...)
        max_x = findmax(all_births)[1]
        max_y = findmax(all_deaths)[1]
        plot!([0, max_y], [0, max_y], label="")
    end

    return plot_ref
end

# ==============================
# Heatmap / layout utilities (from PlottingWrappers)
# ==============================

"""
    plot_square_heatmap(matrix, tick_step, tick_end; plt_title, ...)

Takes matrix and plots it as a heatmap.
"""
function plot_square_heatmap(matrix, tick_step, tick_end;
                                plt_title="", yflip_matrix=true,
                                 plot_params= (dpi=300,
                                                size=(900,800),
                                                lw=1,
                                                thickness_scaling=1,
                                                top_margin= 0,
                                                left_margin=[0 0],
                                                bottom_margin= 0
                                                ),
                                color_palete=:lightrainbow,
                                add_labels=true)
    heat_map = heatmap(matrix, color=color_palete,
                    title=plt_title,
                    size=plot_params.size, dpi=plot_params.dpi,
                    ticks=0:tick_step:tick_end)
    yflip_matrix && plot!(yflip = true)

    if add_labels
        xlabel!("Matrix index")
        ylabel!("Matrix index")
    end

    return heat_map
end

"""
row_plot(bd_plots::Dict; base_h=600, base_w=600, kwargs...)

Plots all plots from the input dictionary in a single row layout.
"""
function row_plot(bd_plots::Dict; base_h = 800, base_w = 800,
                    top_margin= 10mm,
                    left_margin=[10mm 10mm],
                    bottom_margin= 10mm,
                    kwargs...)
    total_dims = length(bd_plots)
    all_plts = tuple()
    for k = 1:total_dims
        all_plts = (all_plts..., bd_plots["β$(k)"])
    end

    nice_plot = plot(all_plts...,
                        layout=(1,total_dims),
                        size=(total_dims*base_w,base_h),
                        thickness_scaling=2,
                        margin=2mm,
                        kwargs...)
    return nice_plot
end

end # module TopologyPreprocessingPlotsExt
