using Eirene
using Plots
include("clique_top_Julia/CliqueTop.jl")
include("VideoProcessing.jl")
include("MatrixToolbox.jl")
include("Settings.jl")

function testing_pariwise_corr()
    do_clique_top = test_params["do_clique_top"]
    do_eirene =     test_params["do_eirene"]
    save_figures = test_params["save_figures"]
    plot_betti_figrues = test_params["plot_betti_figrues"]
    plot_vectorized_video = test_params["plot_vectorized_video"]
    size_limiter = test_params["size_limiter"]
    ind_distrib = test_params["ind_distrib"]
    videos_set = test_params["videos_set"]
    tau_max_set = test_params["tau_max_set"]
    points_per_dim_set = test_params["points_per_dim_set"]
    shifts_set = test_params["shifts_set"]
    patch_params = test_params["patch_params"]
    video_path = test_params["video_path"]
    results_path = test_params["results_path"]
    videos = test_params["videos_names"]


    do_local_corr = false
    do_local_grad = false

    if ind_distrib == "local_corr"
        shift_set = test_params["shift_set"]
        sub_img_size_set = [9]
        do_local_corr = true
        do_local_grad = false
        @debug "Doing local correlation" do_local_corr
    elseif ind_distrib == "local_grad"
        shift_set = [1]
        sub_img_size_set = test_params["sub_img_size_set"]
        do_local_corr = false
        do_local_grad = true
        @debug "Doing local gradient" do_local_grad
    else
        shift_set = [1]
        sub_img_size_set = [9]
        do_local_corr = false
        do_local_grad = false
    end
    @info "Using following distribution: " test_params["ind_distrib"]

    @debug "All videos are: " videos_names
    @debug "Video set is : " videos_set
    for video in videos_set
        choice = videos_names[video]
        @info "Selected video: " choice
        @debug "Path and choice is:" video_path*choice

        video_array = get_video_array_from_file(video_path*choice)
        @info "Array extracted."

        video_dimensions = get_video_dimension(video_array)
        for points_per_dim in points_per_dim_set
            for shift in shift_set, sub_img_size in sub_img_size_set
                if do_local_corr
                    centers = get_local_centers(points_per_dim, video_dimensions, shift, sub_img_size)

                    extracted_pixels_matrix = get_subimg_correlations(video_array, centers, sub_img_size, shift)
                elseif do_local_grad
                    centers = get_local_centers(points_per_dim, video_dimensions, shift, sub_img_size)

                    extracted_pixels_matrix = get_local_gradients(video_array, centers, sub_img_size)
                else
                    indicies_set = get_video_mask(points_per_dim, video_dimensions,  distribution=ind_distrib, patch_params)
                    extracted_pixels_matrix = extract_pixels_from_video(video_array, indicies_set, video_dimensions)
                end
                @info "Pixels extracted."

                vectorized_video = vectorize_video(extracted_pixels_matrix)
                @info "Video is vectorized, proceeding to Pairwise correlation."

                for tau in tau_max_set
                    ## Compute pairwise correlation
                    C_ij = get_pairwise_correlation_matrix(vectorized_video, tau)

                    # set the diagonal to zero
                    for diag_elem in 1:size(C_ij,1)
                        C_ij[diag_elem,diag_elem] = 0
                    end
                    @info "Pairwise correlation finished, proceeding to persistance homology."

                    # Compute persistance homology with CliqueTopJulia
                    size_limiter = test_params["size_limiter"]
                    @debug "using size limiter = " size_limiter

                    if size_limiter > size(C_ij,1)
                        @warn "Used size limiter is larger than matrix dimension: " size_limiter size(C_ij,1)
                        @warn "Using maximal size instead"
                        size_limiter = size(C_ij,1)
                    end

                    @debug "do_clique_top: " do_clique_top
                    @debug "test_params['do_clique_top']: " test_params["do_clique_top"]
                    if do_clique_top
                        @debug pwd()
                        @time c_ij_betti_num, edge_density, persistence_intervals, unbounded_intervals = compute_clique_topology(C_ij[1:size_limiter, 1:size_limiter], edgeDensity = 0.6)
                    end

                    @debug "do_eirene: " do_eirene
                    if do_eirene
                        C = eirene(C_ij[1:size_limiter, 1:size_limiter],maxdim=3,model="vr")
                    end

                    # ---------------------------------------------------------
                    # Plot results
                    @debug "Proceeding to plotting."
                    if plot_vectorized_video
                        vector_plot_ref = heatmap(vectorized_video, color=:grays)
                        if save_figures
                            name = split(choice, ".")[1]
                            name = "vec_" * name * "_sz$(size_limiter)_p$(points_per_dim)_tau$(tau).png"
                            savefig(vector_plot_ref, name)
                        end #save vec
                    end #plot vec

                    if plot_betti_figrues && do_clique_top
                        betti_numbers = c_ij_betti_num
                        title = "Betti curves for pairwise corr. matrix"
                        p1 = plot_betti_numbers(c_ij_betti_num, edge_density, title);

                        heat_map1 = heatmap(C_ij,  color=:lightrainbow, title="Pariwise Correlation matrix, number of points: $(points_per_dim)");

                        betti_plot_clq_ref = plot(p1, heat_map1, layout = (2,1))

                        if save_figures
                            saving_figures(betti_plot_clq_ref, results_cliq_path, choice, points_per_dim, tau, size_limiter)
                        end#save fig
                    end #plot cliq

                    if plot_betti_figrues && do_eirene
                        p1, heat_map1 = plot_eirene_betti_curves(C, C_ij)
                        betti_plot_ei_ref = plot(p1, heat_map1, layout = (2,1))

                        if save_figures
                            saving_figures(betti_plot_ei_ref, results_cliq_path, choice, points_per_dim, tau, size_limiter)
                        end#save fig
                    end #plot eirene
                end #for tau
            end #for shift
        end #for points_per_dim
    end #for video set
    @info "Finished testing"
end #func
