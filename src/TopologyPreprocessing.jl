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

    # MatrixOrganization.jl
    export matrix_poling,
            subsample_matrix,
            add_random_patch

    # MatrixProcessing.jl
    export shift_to_non_negative,
            normalize_to_01,
            diagonal_symmetrize,
            group_distances,
            generate_indices,
            reduce_arrs_to_min_len,
            increase_arrs_to_max_len,
            get_ordered_matrix,
            group_distances,
            generate_indices,
            arr_to_vec,
            cartesianInd_to_vec,
            sort_indices_by_values,
            set_values!

    # BettiCurves.jl
    export get_bettis,
            normalise_bettis,
            get_vectorized_bettis,
            plot_bettis,
            get_bettis_color_palete

    # BarCodes.jl
    export get_barcodes,
            plot_barcodes,
            plot_barcodes!,
            get_birth_death_ratio,
            get_barcode_lifetime,
            get_barcode_max_lifetime,
            boxplot_birth_death,
            boxplot_lifetime,
            get_barcode_max_db_ratios




    include("MatrixOrganization.jl")
    include("MatrixProcessing.jl")
    include("BettiCurves.jl")
    include("Barcodes.jl")

    # include("ImageProcessing.jl")
    # include("PointsSubstitution.jl")


end # module
