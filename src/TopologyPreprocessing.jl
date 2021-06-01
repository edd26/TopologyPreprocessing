module TopologyPreprocessing
    # include("AverageBettis.jl")

    # include("DirOperations.jl")
    # include("GeometricSampling.jl")
    # include("GifGenerator.jl")
    #
    #
    # include("MatrixToolbox.jl")
    # include("PlottingWrappers.jl")
    # include("SavingFigures.jl")
    # include("VideoProcessing.jl")

    export matrix_poling,
            subsample_matrix,
            add_random_patch

    export shift_to_non_negative,
            normalize_to_01,
            diagonal_symmetrize,
            group_distances,
            generate_indices,
            reduce_arrs_to_min_len,
            reduce_arrs_to_max_len

    export get_ordered_matrix,
            group_distances,
            generate_indices,
            arr_to_vec,
            cartesianInd_to_vec,
            sort_index_by_values,
            set_values!

    export get_bettis,
            normalise_bettis,
            get_vectorized_bettis,
            plot_bettis,
            get_bettis_color_palete

    export get_barcodes,
            plot_barcodes,
            get_birth_death_ratio,
            get_barcode_lifetime,
            get_barcode_max_lifetime,
            boxplot_birth_death,
            boxplot_lifetime,
            get_barcode_max_db_ratios

    include("MatrixProcessing.jl")
    include("MatrixOrganization.jl")
    include("BettiCurves.jl")
    include("BarCodes.jl")
    include("ImageProcessing.jl")
    include("PointsSubstitution.jl")


end # module
