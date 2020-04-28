"""
Save figures.
"""
function save_figure(plot_ref, results_path, plot_title; extension=".png" )
    full_path = results_path*plot_title*extension
    savefig(plot_ref, full_path)
    @info "File saved under: " full_path
end


"""
Save betti curves.
"""
function save_betti(plot_ref, results_path, plot_title)
    full_title =  "betti_curves_"*plot_title;
    save_figure(plot_ref, results_path, full_title)
end

"""
Save figures with set of parameters given as 'kwargs'.
"""
function save_figure_with_params(plot_reference, results_path; extension=".png", prefix="", kwargs... )
    plot_title = ""
    kwargs_kyes = keys(kwargs)
    kwargs_vals = values(kwargs)

    total_params = size(collect(kwargs_vals),1)

    for param = 1:total_params
        plot_title *= "$(kwargs_kyes[param])_$(kwargs_vals[param])_"
    end

    full_path = results_path*prefix*plot_title*extension

    # savefig(plot_ref, full_path)
    @info "File saved as: " full_path
end
