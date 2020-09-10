module TopologyPreprocessing

include("AverageBettis.jl")
include("BettiCurves.jl")
include("DirOperations.jl")
include("GeometricSampling.jl")
include("GifGenerator.jl")
include("ImageProcessing.jl")
include("MatrixOrganization.jl")
include("MatrixProcessing.jl")
include("MatrixToolbox.jl")
include("PlottingWrappers.jl")
include("PointsSubstitution.jl")
include("SavingFigures.jl")
include("VideoProcessing.jl")

export matrix_poling,
        subsample_matrix,
        add_random_patch

end # module
