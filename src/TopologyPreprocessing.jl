module TopologyPreprocessing
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
export normalise_bettis,
    vectorize_bettis,
    get_area_under_betti_curve

# BettiCurves.jl — Eirene extension (loaded when Eirene is available)
export get_bettis,
    get_vectorized_bettis,
    get_normalised_barcodes,
    plot_bettis,
    get_bettis_color_palete

# BarCodes.jl — Eirene extension (loaded when Eirene is available)
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

# Eirene-dependent stubs — implemented in ext/TopologyPreprocessingEireneExt.jl
function get_bettis end
function get_barcodes end
function get_vectorized_bettis end

# Plotting stubs — implemented in ext/TopologyPreprocessingPlotsExt.jl
function plot_bettis end
function get_bettis_color_palete end
function plot_barcodes end
function plot_barcodes! end
function boxplot_birth_death end
function boxplot_lifetime end

end # module
