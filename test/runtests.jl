using SafeTestsets

# @testset "TopologyPreprocessing.jl" begin
#     # Write your own tests here.
# end

@safetestset "MatrixOrganization tests" begin include("matrixOrganisation_tests.jl") end
