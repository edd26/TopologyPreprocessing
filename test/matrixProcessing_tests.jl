using TopologyPreprocessing
using Test
using LinearAlgebra

# using MatrixOrganization

@testset "MatrixProcessing.jl -> basics" begin
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
end

@testset "MatrixProcessing.jl -> distances and indices" begin
    # Ints
    positive_matrix = [1 2 3;4 5 6]
    positive_matrix2 = [1. 2 3; 4 5 6]
    positive_3d_matrix = rand(2,3,4)

    negative_matrix = [1 2 -3;-4 5 6]
    negative_matrix2 = [1 2 -3; -4 5 6]
    negative_3d_matrix = rand(2,3,4).*2 .-1

    negative_6d_matrix = rand(2,3,4,5,6,7).*2 .-1

    powers_of_2 = [0 ; [2^k for k in 0:2:8]]
    powers_of_3 = [0 ; [3^k for k in 0:2:8]]

    square_matrix = [1 2 3;
                     4 5 6;
                     7 8 9]

    # ========================================================

    @test length(unique(group_distances(square_matrix,1))) == 1
    @test length(unique(group_distances(square_matrix,2))) == 2
    @test length(unique(group_distances(square_matrix,9))) == 9
    @test_throws DomainError group_distances(square_matrix,20)

    @test length(unique(group_distances(positive_matrix,1))) == 1
    @test length(unique(group_distances(positive_matrix,2))) == 2
    @test length(unique(group_distances(positive_matrix,6))) == 6
    @test_throws DomainError group_distances(positive_matrix,20)

    @test length(unique(group_distances(positive_3d_matrix,2))) == 2
    @test length(unique(group_distances(negative_3d_matrix,2))) == 2
    @test length(unique(group_distances(negative_6d_matrix,2))) == 2

    group_distances(square_matrix,2)

    # ========================================================

    @test length(generate_indices(3))==3^2
    @test length(generate_indices(223))==223^2

    n=3
    @test length(generate_indices(n, symmetry_order=true, include_diagonal=false)) == ((n^2 - n) ÷ 2)
    n=9
    @test length(generate_indices(n, symmetry_order=true, include_diagonal=false)) == ((n^2 - n) ÷ 2)
    index_matrix1 = generate_indices(n, symmetry_order=true, include_diagonal=false)
    @test length(findall(x-> x == CartesianIndex(n,n),index_matrix1)) == 0
    @test length(findall(x-> x == CartesianIndex(n-2,n-1),index_matrix1)) == 1
    @test length(findall(x-> x == CartesianIndex(n-1,n-2),index_matrix1)) == 0

    n=223
    @test length(generate_indices(n, symmetry_order=true)) == ((n^2 - n) ÷ 2)

    n=3
    index_matrix2 = generate_indices(n, symmetry_order=true, include_diagonal=true)
    @test length(index_matrix2) == (((n-1)*n)÷2 +n)
    @test findmax(index_matrix2)[1] == CartesianIndex(n,n)
    @test length(findall(x-> x == CartesianIndex(1,n),index_matrix2)) == 1
    @test length(findall(x-> x == CartesianIndex(1,1),index_matrix2)) == 1
    @test length(findall(x-> x == CartesianIndex(n,n),index_matrix2)) == 1

    n=9
    @test length(generate_indices(n, symmetry_order=true, include_diagonal=true)) == (((n-1)*n)÷2 +n)
    n=223
    @test length(generate_indices(n, symmetry_order=true, include_diagonal=true)) == (((n-1)*n)÷2 +n)


end


@testset "MatrixProcessing.jl -> arrays of arrays" begin
    # Need to be taken from Bettis
    # arr_of_arrs

    # reduce_arrs_to_min_len(arr_of_arrs)

end

@testset "MatrixProcessing.jl -> matrix ordering" begin
    power_matrix = zeros(Int, (2,3,4))
    for k = 1:4
        power_matrix[:,:,k] = reshape([(k+1)^n for n =1:6], (2,3))
    end

    ordered_power_matrix = copy(power_matrix)
    ordered_power_matrix[:, :, 1] =[1  6  12;
                                    3  8  13]
    ordered_power_matrix[:, :, 2] =[2  11  17;
                                        7  15  20]
    ordered_power_matrix[:, :, 3] = [4  14  21;
                                        9  18  23]
    ordered_power_matrix[:, :, 4] =[5  16  22;
                                    10  19  24]

    square_matrix1 = [1 2 3;
                     4 5 6;
                     7 8 9]

     square_matrix2 = [1 4 7;
                        2 5 8;
                        3 6 9]

    @test sum(get_ordered_matrix2(square_matrix1) .== square_matrix1) == 9
    @test sum(get_ordered_matrix2(square_matrix2) .== square_matrix2) == 9

    @test get_ordered_matrix2(power_matrix) == ordered_power_matrix

    @test !isempty(get_ordered_matrix(power_matrix) == ordered_power_matrix)




end
