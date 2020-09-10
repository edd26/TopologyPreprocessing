using TopologyPreprocessing
using Test
using LinearAlgebra

# using MatrixOrganization

@testset "MatrixProcessing.jl -> " begin
    # Ints
    positive_matrix = [1 2 3;
                       4 5 6]
    negative_matrix = [1 2 -3;
                       -4 5 6]
    @test shift_to_non_negative(positive_matrix) == positive_matrix
    @test isempty(findall(x->x<0, shift_to_non_negative(negative_matrix)))

    # Floats
    positive_matrix2 = [1. 2 3; 4 5 6]
    negative_matrix2 = [1 2 -3; -4 5 6]
    @test shift_to_non_negative(positive_matrix2) == positive_matrix2
    @test isempty(findall(x->x<0, shift_to_non_negative(negative_matrix2)))

    positive_3d_matrix = rand(2,3,4)
    negative_3d_matrix = rand(2,3,4).*2 .-1
    @test shift_to_non_negative(positive_3d_matrix) == positive_3d_matrix
    @test isempty(findall(x->x<0, shift_to_non_negative(negative_3d_matrix)))

    # =====================

    @test findmin(normalize_to_01(negative_matrix))[1] == 0
    @test findmax(normalize_to_01(negative_matrix))[1] == 1

    powers_of_2 = [0 ; [2^k for k in 0:2:8]]
    powers_of_3 = [0 ; [3^k for k in 0:2:8]]
    @test findmin(normalize_to_01(powers_of_2, use_factor=true))[1] == 0
    @test findmax(normalize_to_01(powers_of_2, use_factor=true))[1] == 1

    @test findmin(normalize_to_01(powers_of_3, use_factor=true))[1] == 0
    @test findmax(normalize_to_01(powers_of_3, use_factor=true))[1] > 1

    @test findmin(normalize_to_01(powers_of_3, use_factor=true, norm_factor=3^8))[1] == 0
    @test findmax(normalize_to_01(powers_of_3, use_factor=true, norm_factor=3^8))[1] == 1

    # =====================

    square_matrix = [1 2 3;
                     4 5 6;
                     7 8 9]
    @test issymmetric(diagonal_symmetrize(square_matrix))
    @test LinearAlgebra.checksquare(diagonal_symmetrize(square_matrix)) == 3

    @test issymmetric(diagonal_symmetrize(positive_matrix))
    @test LinearAlgebra.checksquare(diagonal_symmetrize(positive_matrix)) == 2

    @test diagonal_symmetrize(square_matrix)[end,1] == square_matrix[1,end]
    @test diagonal_symmetrize(square_matrix)[2,1] == square_matrix[1,2]

    @test diagonal_symmetrize(square_matrix, below_over_upper=true)[1,end] == square_matrix[end,1]
    @test diagonal_symmetrize(square_matrix, below_over_upper=true)[1,2] == square_matrix[2,1]

    # =====================

    

end
