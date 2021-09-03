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
    @test length(generate_indices(n, symmetry_order=true, include_diagonal=false)) == ((n^2 - n) ÷ 2)
    @test length(generate_indices(n, symmetry_order=true, include_diagonal=true)) == ((n^2 - n) ÷ 2 + n)

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

    n=9
    @test length(generate_indices((n,n), symmetry_order=true, include_diagonal=true)) == (((n-1)*n)÷2 +n)
    n=223
    @test length(generate_indices((n,n), symmetry_order=true, include_diagonal=true)) == (((n-1)*n)÷2 +n)
end


# @testset "MatrixProcessing.jl -> arrays of arrays" begin
#     # Need to be taken from Bettis
#     # arr_of_arrs
#
#     # reduce_arrs_to_min_len(arr_of_arrs)
#
# end

@testset "MatrixProcessing.jl -> matrix ordering helping functions" begin
    let testing_matrix0 = Array{Float64,2}(undef, 2, 3),
        testing_matrix1 = [1 2 3; 4 5 6; 7 8 9],
        testing_matrix2 = ones((2,3,4))
        testing_matrix3 = [1 2 3; 4 5 6]
        testing_matrix4 = [1, 4, 7, 2, 5, 8, 3, 6, 9]

        @test arr_to_vec(testing_matrix0) isa Vector
        @test length(testing_matrix0) == length(arr_to_vec(testing_matrix0))

        @test arr_to_vec(testing_matrix1) isa Vector
        @test length(testing_matrix1) == length(arr_to_vec(testing_matrix1))
        @test arr_to_vec(testing_matrix1) == [1, 4, 7, 2, 5, 8, 3, 6, 9]

        @test arr_to_vec(testing_matrix2) isa Vector
        @test length(testing_matrix2) == length(arr_to_vec(testing_matrix2))

        @test arr_to_vec(testing_matrix3) isa Vector
        @test length(testing_matrix3) == length(arr_to_vec(testing_matrix3))
        @test arr_to_vec(testing_matrix3) == [1, 4, 2, 5, 3, 6]

        @test arr_to_vec(testing_matrix4) isa Vector
        @test length(testing_matrix4) == length(arr_to_vec(testing_matrix4))
        @test arr_to_vec(testing_matrix4) == [1, 4, 7, 2, 5, 8, 3, 6, 9]
    end

    let testing_matrix1 = CartesianIndices((2,2)),
        testing_matrix2 = CartesianIndices((9,1))

        @test cartesianInd_to_vec(testing_matrix1) isa Vector
        @test length(cartesianInd_to_vec(testing_matrix1)) == length(testing_matrix1)

        @test cartesianInd_to_vec(testing_matrix2) isa Vector
        @test length(cartesianInd_to_vec(testing_matrix2)) == length(testing_matrix2)
    end

    let vals_matrix1a = [1 2 3; 4 5 6],
        vals_matrix1b = [6 5 4; 3 2 1],
        vals_matrix1c = [4 5 6; 1 2 3]
        let index_matrix1a = [CartesianIndex(1,1), CartesianIndex(1,2), CartesianIndex(1,3), CartesianIndex(2,1), CartesianIndex(2,2), CartesianIndex(2,3)]
            @test length(sort_index_by_values(vals_matrix1a, index_matrix1a)) == length(vals_matrix1a)
            @test sort_index_by_values(vals_matrix1a, index_matrix1a) == collect(1:6)

            @test length(sort_index_by_values(vals_matrix1b, index_matrix1a)) == length(vals_matrix1b)
            @test sort_index_by_values(vals_matrix1b, index_matrix1a) == collect(6:-1:1)

            @test length(sort_index_by_values(vals_matrix1c, index_matrix1a)) == length(vals_matrix1c)
            @test sort_index_by_values(vals_matrix1c, index_matrix1a) == [4, 5, 6, 1, 2, 3]
        end

        let index_matrix1b = CartesianIndices((2,3))
            @test_throws TypeError sort_index_by_values(vals_matrix1a, index_matrix1b)
        end
    end

    let vals_matrix2 = [1 2 3; 4 5 6; 7 8 9],
        index_matrix2a = CartesianIndices((3,3)),
        index_matrix2b =   [CartesianIndex(1,1) CartesianIndex(1,2) CartesianIndex(1,3);
                            CartesianIndex(2,1) CartesianIndex(2,2) CartesianIndex(2,3);
                            CartesianIndex(3,1) CartesianIndex(3,2) CartesianIndex(3,3)],
        index_matrix2c =   [CartesianIndex(1,1), CartesianIndex(1,2), CartesianIndex(1,3),
                            CartesianIndex(2,1), CartesianIndex(2,2), CartesianIndex(2,3),
                            CartesianIndex(3,1), CartesianIndex(3,2), CartesianIndex(3,3)]

        @test_throws TypeError sort_index_by_values(vals_matrix2, index_matrix2a)
        @test_throws TypeError sort_index_by_values(vals_matrix2, index_matrix2b)

        @test sort_index_by_values(vals_matrix2, index_matrix2c) isa Vector
        @test sort_index_by_values(vals_matrix2, index_matrix2c) == 1:9
        @test length(sort_index_by_values(vals_matrix2, index_matrix2c)) == length(vals_matrix2)
    end

    let vals_matrix3 = [1, 4, 7, 2, 5, 8, 3, 6, 9],
        index_matrix3a = CartesianIndices((9,1)),
        index_matrix3b = CartesianIndices((9,)),
        index_matrix3c = [1, 4, 7, 2, 5, 8, 3, 6, 9]

        @test_throws TypeError sort_index_by_values(vals_matrix3, index_matrix3a)
        @test_throws TypeError sort_index_by_values(vals_matrix3, index_matrix3b)

        @test sort_index_by_values(vals_matrix3, index_matrix3c) isa Vector
        @test sort_index_by_values(vals_matrix3, index_matrix3c) == 1:9
        @test length(sort_index_by_values(vals_matrix3, index_matrix3c)) == length(vals_matrix3)
    end

    let target_coords1 = CartesianIndex(2,3),
        target_value = -20

        let some_matrix = [1 2 3; 4 5 6; 7 8 9]

            set_values!(some_matrix, target_coords1, target_value; do_symmetry=false)
            @test some_matrix[target_coords1] == target_value
        end

        let some_matrix = [1 2 3; 4 5 6; 7 8 9]

            another_matrix = set_values!(some_matrix, target_coords1, target_value; do_symmetry=false)
            @test some_matrix[target_coords1] == target_value
            @test another_matrix[target_coords1] == target_value
            @test another_matrix === some_matrix
        end

        let some_matrix = [1 2 3; 4 5 6; 7 8 9]

            another_matrix = set_values!(some_matrix, target_coords1, target_value; do_symmetry=true)
            @test some_matrix[target_coords1] == target_value
            @test some_matrix[target_coords1[1],target_coords1[2]] == target_value
            @test some_matrix[target_coords1[2],target_coords1[1]] == target_value
            @test another_matrix === some_matrix
        end

        let some_matrix = [1 2 3; 4 5 6; 7 8 9],
            some_matrix2 = [1 2 3; 4 5 6; 7 8 9],
            target_coords2 = CartesianIndex(8,9)

            @test_throws BoundsError set_values!(some_matrix, target_coords2, target_value; do_symmetry=false)
            @test some_matrix == some_matrix2
        end

        let some_matrix = [1 2 3; 4 5 6; 7 8 9],
            target_coords3 = CartesianIndices((2,2))

            @test_throws MethodError set_values!(some_matrix, target_coords3, target_value; do_symmetry=false)
        end
    end
end

@testset "MatrixProcessing.jl -> matrix ordering" begin

    """
        check_for_min_val_position(input_matrix::Matrix; force_symmetry=false, assign_same_values=false)

    Takes an 'input_matrix' and checks if its ordered form has the minimum value
    in the same position as the minimum value of 'input_matrix'.

    """
    function check_for_min_val_position(input_matrix::Matrix; force_symmetry=false, assign_same_values=false)
        # check if min val in uniquely value matrix is in the same position
        ordered_matrix = get_ordered_matrix(input_matrix, force_symmetry=force_symmetry, assign_same_values=assign_same_values)
        if issymmetric(input_matrix) || force_symmetry
            off_diag_ind = generate_indices(size(input_matrix),include_diagonal=false)
        else
            off_diag_ind = generate_indices(size(input_matrix),include_diagonal=true)
        end
        min_val = findmin(input_matrix[off_diag_ind])[1]
        # min_orig = findmin(input_matrix[off_diag_ind])[2]
        all_min_input = findall(x->x==min_val,input_matrix[off_diag_ind])
        all_min_ordered = findall(x->x==1,ordered_matrix[off_diag_ind])
        return off_diag_ind[all_min_input] == off_diag_ind[all_min_ordered]
    end

    # ==
    let power_matrix = zeros(Int, (2,3,4))
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

        @test !isempty(get_ordered_matrix(power_matrix) == ordered_power_matrix)
        @test get_ordered_matrix(power_matrix) == ordered_power_matrix

        # @test_throws MethodError get_ordered_matrix(power_matrix)
    end

    # ==
    let square_matrix1 = [1 2 3; 4 5 6; 7 8 9],
        square_matrix2 = [1 4 7; 2 5 8; 3 6 9]
        @test get_ordered_matrix(10 .*square_matrix1) == (square_matrix1 )
        @test get_ordered_matrix(10 .*square_matrix2) == (square_matrix2 )
        @test sum(get_ordered_matrix(square_matrix1) .== square_matrix1 ) == 9
        @test sum(get_ordered_matrix(square_matrix2) .== square_matrix2 ) == 9
        @test get_ordered_matrix(10 .*square_matrix1, force_symmetry=true) == get_ordered_matrix(square_matrix1, force_symmetry=true)
        @test get_ordered_matrix(10 .*square_matrix2, force_symmetry=true) == get_ordered_matrix(square_matrix2, force_symmetry=true)

        # check if min val in uniquely value matrix is in the same position
        let input_matrix = 10square_matrix1
            @test check_for_min_val_position(input_matrix)
        end
    end

    # ==
    let square_matrix_same_vals1 = [1 2 3; 3 4 5; 6 7 8],
        square_matrix_same_vals2 = [1 3 6; 2 4 7; 3 5 8]

        @test get_ordered_matrix(10 .*square_matrix_same_vals1, assign_same_values=true) == (square_matrix_same_vals1 )
        @test get_ordered_matrix(10 .*square_matrix_same_vals2, assign_same_values=true) == (square_matrix_same_vals2 )

        # forcing symmetry test
        some_ord_mat = get_ordered_matrix(10 .*square_matrix_same_vals1, force_symmetry=true, assign_same_values=false)
        # remove 1, because 0 is not off diagonal
        @test length(unique(some_ord_mat))-1 == (size(square_matrix_same_vals1,1)*(size(square_matrix_same_vals1,1)-1))/2
    end

    # ==
    let square_matrix_same_vals3 = [1 2 3; 3 4 3; 5 6 7],
        square_matrix_same_vals4 = [1 3 3; 2 4 6; 3 3 7]
        @test get_ordered_matrix(10 .*square_matrix_same_vals3, force_symmetry=true, assign_same_values=true) ==
                            get_ordered_matrix(square_matrix_same_vals3, force_symmetry=true, assign_same_values=true)
        @test get_ordered_matrix(10 .*square_matrix_same_vals4, force_symmetry=true, assign_same_values=true) ==
                            get_ordered_matrix(square_matrix_same_vals4, force_symmetry=true, assign_same_values=true)
    end

    # ==================
    # Tests on symmetric matrices
    let test_mat_b1 = [ 1  1  1  4  5  9 ;
                        1  1  2  3  6  10;
                        1  2  1  3  7  11;
                        4  3  3  1  8  12;
                        5  6  7  8  1  13;
                        9  10 11 12 13 1 ;]

        # no same values
        let test_mat_b1_ord1 = [  0   1   2   6   7  11;
                                  1   0   3   4   8  12;
                                  2   3   0   5   9  13;
                                  6   4   5   0  10  14;
                                  7   8   9  10   0  15;
                                 11  12  13  14  15   0],
            test_mat_b1_indices = generate_indices(size(test_mat_b1))

            ord_mat_b1_1 = get_ordered_matrix(10 .*test_mat_b1, force_symmetry=false, assign_same_values=false)
            @test issymmetric(ord_mat_b1_1)
            @test ord_mat_b1_1[test_mat_b1_indices] == test_mat_b1_ord1[test_mat_b1_indices]

            ord_mat_b1_2 = get_ordered_matrix(10 .*test_mat_b1, force_symmetry=true, assign_same_values=false)
            @test issymmetric(ord_mat_b1_2)
            @test ord_mat_b1_2[test_mat_b1_indices] == test_mat_b1_ord1[test_mat_b1_indices]
        end

        # assign same values
        let test_mat_b1_ord2 = [ 0  1  1  4  5  9 ;
                                 1  0  2  3  6  10;
                                 1  2  0  3  7  11;
                                 4  3  3  0  8  12;
                                 5  6  7  8  0  13;
                                 9 10  11 12 13 0 ;]

            let ord_mat_b1_3 = get_ordered_matrix(10 .*test_mat_b1, force_symmetry=false, assign_same_values=true)
                @test issymmetric(ord_mat_b1_3)
                # Removed, because ordered diagonal is all 1: @test ord_mat_b1_3 == (test_mat_b1 )
                @test ord_mat_b1_3 == test_mat_b1_ord2
            end

            let ord_mat_b1_4 = get_ordered_matrix(10 .*test_mat_b1, force_symmetry=true, assign_same_values=true)
                @test issymmetric(ord_mat_b1_4)
                # Removed, because ordered diagonal is all 1: @test ord_mat_b1_4 == (test_mat_b1 )
                @test ord_mat_b1_4 == test_mat_b1_ord2
            end
        end

        let input_matrix = 10 .*test_mat_b1
            @test check_for_min_val_position(input_matrix; force_symmetry=false, assign_same_values=true)
        end
    end

    # ==
    # Non-symmetric matrix test
    let test_mat_b2 = [ 1  1  3  4  5  9 ;
                        1  1  2  3  6  10;
                        14 2  1  3  7  11;
                        4  15 3  1  8  12;
                        5  5  7  8  1  13;
                        9  10 11 12 13 1 ;]

        let ord_mat_b2_1 = get_ordered_matrix(10 .*test_mat_b2, force_symmetry=false, assign_same_values=false)
            @test !issymmetric(ord_mat_b2_1)
            @test findall(x->x<=8,ord_mat_b2_1) == findall(x->x==1,test_mat_b2) # all values with one are first 8 values used for ordering
            @test length(unique(ord_mat_b2_1)) == length(test_mat_b2) #check if all values are unique
        end

        let test_mat_b2_ord1 = [0  1  6  4  7  11;
                                1  0  2  5  8  12;
                                6  2  0  3  9  13;
                                4  5  3  0 10  14;
                                7  8  9 10  0  15;
                                11 12 13 14 15 0 ;]

            test_mat_b2_indices = generate_indices(size(test_mat_b2))
            filter!(x->x!=CartesianIndex(1,3), test_mat_b2_indices)
            filter!(x->x!=CartesianIndex(1,4), test_mat_b2_indices)
            filter!(x->x!=CartesianIndex(2,4), test_mat_b2_indices)
            filter!(x->x!=CartesianIndex(3,4), test_mat_b2_indices)
            filter!(x->x!=CartesianIndex(3,1), test_mat_b2_indices)
            filter!(x->x!=CartesianIndex(4,1), test_mat_b2_indices)
            filter!(x->x!=CartesianIndex(4,2), test_mat_b2_indices)
            filter!(x->x!=CartesianIndex(4,3), test_mat_b2_indices)

            ord_mat_b2_2 = get_ordered_matrix(10 .*test_mat_b2, force_symmetry=true, assign_same_values=false)
            @test issymmetric(ord_mat_b2_2)
            @test ord_mat_b2_2[test_mat_b2_indices] == test_mat_b2_ord1[test_mat_b2_indices]
            # forcing symmetry test:
            @test length(unique(ord_mat_b2_2))-1 == (size(test_mat_b2,1)*(size(test_mat_b2,1)-1))/2
        end

        # TODO to fix tests, add 1 to all values off diagonal for the input matrix
        let test_mat_b2_ord2 = [0  0  2  3  4  8;
                                0  0  1  2  5  9;
                                13 1  0  2  6  10;
                                3  14 2  0  7  11;
                                4  4  6  7  0  12;
                                8  9  10 11 12 0 ;]
            ord_mat_b2_3 = get_ordered_matrix(10 .*test_mat_b2, force_symmetry=false, assign_same_values=true)
            @test_skip !issymmetric(ord_mat_b2_3)
            @test_skip  ord_mat_b2_3 == test_mat_b2_ord2
        end

        # TODO to fix tests, add 1 to all values off diagonal for the input matrix
        let test_mat_b2_ord3 = [0  0  2   3   4  8
                                0  0  1   2   5  9
                                2  1  0   2   6  10
                                3  2  2   0   7  11
                                4  5  6   7   0  12
                                8  9  10  11  12 0 ;]
            ord_mat_b2_4 = get_ordered_matrix(10 .*test_mat_b2, force_symmetry=true, assign_same_values=true)
            @test_skip issymmetric(ord_mat_b2_4)
            @test_skip ord_mat_b2_4 == test_mat_b2_ord3
        end
    end

    # ==
    let test_mat_b3 = [ 1  1  3  4  5  9 ;
                        1  1  2  3  6  10;
                        3  2  1  3  7  11;
                        4  3  3  1  8  12;
                        5  6  7  8  1  13;
                        9  10 11 12 13 1 ;]

        let test_mat_b3_ord1 = [ 0  0  5  3  6  10;
                                    0  0  1  4  7  11;
                                    5  1  0  2  8  12;
                                    3  4  2  0  9  13;
                                    6  7  8  9  0  14;
                                    10 11 12 13 14 0 ;],
          test_mat_b3_indices = generate_indices(size(test_mat_b3))
            filter!(x->x!=CartesianIndex(1,3), test_mat_b3_indices)
            filter!(x->x!=CartesianIndex(1,4), test_mat_b3_indices)
            filter!(x->x!=CartesianIndex(2,4), test_mat_b3_indices)
            filter!(x->x!=CartesianIndex(3,4), test_mat_b3_indices)
            filter!(x->x!=CartesianIndex(3,1), test_mat_b3_indices)
            filter!(x->x!=CartesianIndex(4,1), test_mat_b3_indices)
            filter!(x->x!=CartesianIndex(4,2), test_mat_b3_indices)
            filter!(x->x!=CartesianIndex(4,3), test_mat_b3_indices)

            ord_mat_b3_1 = get_ordered_matrix(10 .*test_mat_b3, force_symmetry=false, assign_same_values=false)
            @test issymmetric(ord_mat_b3_1)
            @test_skip ord_mat_b3_1[test_mat_b3_indices] == test_mat_b3_ord1[test_mat_b3_indices]

            ord_mat_b3_2 = get_ordered_matrix(10 .*test_mat_b3, force_symmetry=true, assign_same_values=false)
            @test issymmetric(ord_mat_b3_2)
            @test_skip ord_mat_b3_2[test_mat_b3_indices] == test_mat_b3_ord1[test_mat_b3_indices]
        end

        # TODO remove tests that do not add anything new for testing and are just another similar case
        let test_mat_b3_ord2 = [ 0  0  2  3  4  8 ;
                                 0  0  1  2  5  9 ;
                                 2  1  0  2  6  10;
                                 3  2  2  0  7  11;
                                 4  5  6  7  0  12;
                                 8  9  10 11 12 0 ;]

            ord_mat_b3_3 = get_ordered_matrix(10 .*test_mat_b3, force_symmetry=false, assign_same_values=true)
            @test issymmetric(ord_mat_b3_3)
            @test_skip ord_mat_b3_3 == (test_mat_b3 .-1)
            @test_skip ord_mat_b3_3 == test_mat_b3_ord2

            ord_mat_b3_4 = get_ordered_matrix(10 .*test_mat_b3, force_symmetry=true, assign_same_values=true)
            @test issymmetric(ord_mat_b3_4)
            @test_skip ord_mat_b3_4 == (test_mat_b3 .-1)
            @test_skip ord_mat_b3_4 == test_mat_b3_ord2
        end

        let input_matrix = 10 .*test_mat_b3
            @test check_for_min_val_position(input_matrix; force_symmetry=false, assign_same_values=true)
        end
    end

    # ==
    let test_mat_b4 = [  1  1  41 4  5  9  13 17 25 33;
                        1  1  2  42 6  10 14 18 26 34;
                        41 2  1  3  7  11 15 19 27 35;
                        4  42 3  1  8  12 16 20 28 36;
                        5  6  7  8  1  21 43 24 29 37;
                        9  10 11 12 21 1  22 44 30 38;
                        13 14 15 16 43 22 1  23 31 39;
                        17 18 19 20 24 44 23 1  32 40;
                        25 26 27 28 29 30 31 32 1  45;
                        33 34 35 36 37 38 39 40 45 1;]

        @test_skip get_ordered_matrix(10 .*test_mat_b4, force_symmetry=false, assign_same_values=false) == (test_mat_b4 .-1)
        @test_skip get_ordered_matrix(10 .*test_mat_b4, force_symmetry=false, assign_same_values=true) == (test_mat_b4 .-1)

        let ord_mat = get_ordered_matrix(10 .*test_mat_b4, force_symmetry=true, assign_same_values=false)
            @test issymmetric(ord_mat)
            @test_skip ord_mat == (test_mat_b4 .-1)
        end

        let ord_mat = get_ordered_matrix(10 .*test_mat_b4, force_symmetry=true, assign_same_values=true)
            @test issymmetric(ord_mat)
            @test_skip  ord_mat == (test_mat_b4 .-1)
        end

        let input_matrix = 10test_mat_b4
            @test check_for_min_val_position(input_matrix)
        end
    end

    # ==
    let test_mat_b5 = -[1  1  3  4  5  9 ;
                        1  1  2  3  6  10;
                        14 2  1  3  7  11;
                        4  15 3  1  8  12;
                        5  5  7  8  1  13;
                        9  10 11 12 13 1 ;]
        let ord_mat_b5_1 = get_ordered_matrix(10 .*test_mat_b5, force_symmetry=false, assign_same_values=false)
            @test !issymmetric(ord_mat_b5_1)
            @test_skip findall(x->x>=28,ord_mat_b5_1) == findall(x->x==-1,test_mat_b5) # all values with one are first 8 values used for ordering
            @test length(unique(ord_mat_b5_1)) == length(test_mat_b5) #check if all values are unique
        end

        let test_mat_b5_ord1 = [ 0  14  9 11  8  4
                                 14 0  13 10  7  3
                                 9  13 0  12  6  2
                                 11 10 12 0   5  1
                                 8  7  6  5   0  0
                                 4  3  2  1   0  0 ;]

            test_mat_b5_indices = generate_indices(size(test_mat_b5))
            filter!(x->x!=CartesianIndex(1,3), test_mat_b5_indices)
            filter!(x->x!=CartesianIndex(1,4), test_mat_b5_indices)
            filter!(x->x!=CartesianIndex(2,4), test_mat_b5_indices)
            filter!(x->x!=CartesianIndex(3,4), test_mat_b5_indices)
            filter!(x->x!=CartesianIndex(3,1), test_mat_b5_indices)
            filter!(x->x!=CartesianIndex(4,1), test_mat_b5_indices)
            filter!(x->x!=CartesianIndex(4,2), test_mat_b5_indices)
            filter!(x->x!=CartesianIndex(4,3), test_mat_b5_indices)

            ord_mat_b5_2 = get_ordered_matrix(10 .*test_mat_b5, force_symmetry=true, assign_same_values=false)

            @test issymmetric(ord_mat_b5_2)
            @test_skip ord_mat_b5_2[test_mat_b5_indices] == test_mat_b5_ord1[test_mat_b5_indices]
            # forcing symmetry test:
            @test length(unique(ord_mat_b5_2))-1 == (size(test_mat_b5,1)*(size(test_mat_b5,1)-1))/2
        end

        let test_mat_b5_ord2 = [14 14 12 11 10 6;
                                14 14 13 12 9  5;
                                1  13 14 12 8  4;
                                11 0  12 14 7  3;
                                10 10 8  7  14 2;
                                6  5  4  3  2  14],
            ord_mat_b5_3 = get_ordered_matrix(10 .*test_mat_b5, force_symmetry=false, assign_same_values=true)

            @test !issymmetric(ord_mat_b5_3)
            @test_skip ord_mat_b5_3 == test_mat_b5_ord2
        end

        let test_mat_b5_ord3 = [0  12 10 9  8  4;
                                12 0  11 10 7  3;
                                10 11 0  10 6  2;
                                9  10 10 0  5  1;
                                8  7  6  5  0  0;
                                4  3  2  1  0  0]
            ord_mat_b5_4 = get_ordered_matrix(10 .*test_mat_b5, force_symmetry=true, assign_same_values=true)
            @test issymmetric(ord_mat_b5_4)
            @test_skip ord_mat_b5_4 == test_mat_b5_ord3
        end
    end
    # ==================
end
