using SafeTestsets

@safetestset "MatrixProcessing tests" begin include("matrixProcessing_tests.jl") end

@safetestset "MatrixOrganization tests" begin include("matrixOrganisation_tests.jl") end

@safetestset "BettiCurves tests" begin include("bettiCurves_tests.jl") end
