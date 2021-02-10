using TopologyPreprocessing
using Test

# using MatrixOrganization

@testset "MatrixOrganization.jl -> matrix pooling" begin
	in_vector = [1 2 3 4 3 2 1]
	in_matrix = [1 2 3; 5 6 7; 8 9 0]
	in_matrix2 = [1 2 3; 5 6 7; 8 9 0]

	sqr_matrix0 = [1 2 3; 2 6 7; 3 7 0]
	resulting_matrix0_2 = [1 2 3; 2 6 7; 3 7 0]

	sqr_matrix1 = [	0	1	13	4	5	9;
					1	0	2	14	6	10;
					13	2	0	3	7	11;
					4	14	3	0	8	12;
					5	6	7	8	0	15;
					9	10	11	12	15	0]
	resulting_matrix1_2 = [0   1  14  14  10  10;
						 1   0  14  14  10  10;
						 14  14   0   3  12  12;
						 14  14   3   0  12  12;
						 10  10  12  12   0  15;
						 10  10  12  12  15   0]

	 sqr_matrix2 = [ 0	1	113	4	5	9	13	17	81	82	83	84;
						1	0	2	114	6	10	14	18	85	86	87	88;
						113	2	0	3	7	11	15	19	89	90	91	92;
						4	114	3	0	8	12	16	20	93	94	95	96;
						5	6	7	8	0	21	115	24	29	30	31	32;
						9	10	11	12	21	0	22	116	33	34	35	36;
						13	14	15	16	115	22	0	23	37	38	39	40;
						17	18	19	20	24	116	23	0	41	42	43	44;
						81	85	89	93	29	33	37	41	0	25	117	28;
						82	86	90	94	30	34	38	42	25	0	26	118;
						83	87	91	95	31	35	39	43	117	26	0	27;
						84	88	92	96	32	36	40	44	28	118	27	0]
		result_matrix2_2 = [   0    1  114  114  10  10   18   18   86   86   88   88;
								  1    0  114  114  10  10   18   18   86   86   88   88;
								114  114    0    3  12  12   20   20   94   94   96   96;
								114  114    3    0  12  12   20   20   94   94   96   96;
								 10   10   12   12   0  21  116  116   34   34   36   36;
								 10   10   12   12  21   0  116  116   34   34   36   36;
							  	 18   18   20   20 116 116    0   23   42   42   44   44;
								 18   18   20   20 116 116   23    0   42   42   44   44;
								 86   86   94   94  34  34   42   42    0   25  118  118;
								 86   86   94   94  34  34   42   42   25    0  118  118;
								 88   88   96   96  36  36   44   44  118  118    0   27;
								 88   88   96   96  36  36   44   44  118  118   27    0]
		result_matrix2_3 = [ 0	1	113	114	114	114	89	89	89	92	92	92;
								1	0	2	114	114	114	89	89	89	92	92	92;
								113	2	0	114	114	114	89	89	89	92	92	92;
								114	114	114	0	8	12	116	116	116	96	96	96;
								114	114	114	8	0	21	116	116	116	96	96	96;
								114	114	114	12	21	0	116	116	116	96	96	96;
								89	89	89	116	116	116	0	23	37	117	117	117;
								89	89	89	116	116	116	23	0	41	117	117	117;
								89	89	89	116	116	116	37	41	0	117	117	117;
								92	92	92	96	96	96	117	117	117	0	26	118;
								92	92	92	96	96	96	117	117	117	26	0	27;
								92	92	92	96	96	96	117	117	117	118	27	0]
		result_matrix2_4 = [ 0	1	114	114	20	20	20	20	96	96	96	96;
								1	0	114	114	20	20	20	20	96	96	96	96;
								114	114	0	3	20	20	20	20	96	96	96	96;
								114	114	3	0	20	20	20	20	96	96	96	96;
								20	20	20	20	0	21	116	116	44	44	44	44;
								20	20	20	20	21	0	116	116	44	44	44	44;
								20	20	20	20	116	116	0	23	44	44	44	44;
								20	20	20	20	116	116	23	0	44	44	44	44;
								96	96	96	96	44	44	44	44	0	25	118	118;
								96	96	96	96	44	44	44	44	25	0	118	118;
								96	96	96	96	44	44	44	44	118	118	0	27;
								96	96	96	96	44	44	44	44	118	118	27	0]

	@test matrix_poling(in_vector) !== in_vector
	@test matrix_poling(in_vector, method="max_pooling") == [4 4 4 4 4 4 4]
	# @test matrix_poling!(in_vector) === in_vector
	# @test matrix_poling!(in_vector) == [4 4 4 4 4 4 4]

	@test matrix_poling(in_matrix) !== in_matrix
	@test matrix_poling(in_matrix, method="max_pooling") == 9 .* ones(size(in_matrix))
	# @test matrix_poling!(in_matrix) === in_matrix
	# @test matrix_poling!(in_matrix) == 9 .* ones(size(in_matrix))

	@test matrix_poling(in_matrix2[1:2,1:2], method="max_pooling") == 6 .* ones(size(in_matrix2[1:2,1:2]))
	@test matrix_poling(in_matrix2[1:2,1:2], method="max_pooling") != in_matrix2[1:2,1:2]

	# ====
	# Subsampling matrix

	# Function is supposed to work only on upper half, and here the upper half is too small, so there are no operations
	@test subsample_matrix(sqr_matrix0, subsamp_size=2, method="max_pooling") == resulting_matrix0_2

	@test subsample_matrix(sqr_matrix1, subsamp_size=2, method="max_pooling") == resulting_matrix1_2

	@test subsample_matrix(sqr_matrix2, subsamp_size=2, method="max_pooling") == result_matrix2_2
	@test subsample_matrix(sqr_matrix2, subsamp_size=3, method="max_pooling") == result_matrix2_3
	@test subsample_matrix(sqr_matrix2, subsamp_size=4, method="max_pooling") == result_matrix2_4
end

@testset "MatrixOrganization.jl -> add_random_patch" begin
	# TODO set seed for add_random_path
	# TODO Seed has to be set for this test
	in_vector = [1, 2, 3, 4, 3, 2, 1]

	sqr_matrix0 = [ 1 2 3;
	 				2 6 7;
					3 7 0]
	sqr_matrix1 = [1 2 3  4  5;
					 2 1 6  7  8;
					 3 6 1  9  10;
					 4 7 9  1  11;
					 5 8 10 11 1]
	sqr_matrix2 = [	0	1	13	4	5	9;
					1	0	2	14	6	10;
					13	2	0	3	7	11;
					4	14	3	0	8	12;
					5	6	7	8	0	15;
					9	10	11	12	15	0]
	 sqr_matrix3 = [ 0	1	113	4	5	9	13	17	81	82	83	84;
 					1	0	2	114	6	10	14	18	85	86	87	88;
 					113	2	0	3	7	11	15	19	89	90	91	92;
 					4	114	3	0	8	12	16	20	93	94	95	96;
 					5	6	7	8	0	21	115	24	29	30	31	32;
 					9	10	11	12	21	0	22	116	33	34	35	36;
 					13	14	15	16	115	22	0	23	37	38	39	40;
 					17	18	19	20	24	116	23	0	41	42	43	44;
 					81	85	89	93	29	33	37	41	0	25	117	28;
 					82	86	90	94	30	34	38	42	25	0	26	118;
 					83	87	91	95	31	35	39	43	117	26	0	27;
 					84	88	92	96	32	36	40	44	28	118	27	0]
	function get_unchanged_indices(input_matrix,ind)
		indices = CartesianIndices(size(input_matrix))
		indices = findall(x->x!=ind,indices)
		for i = ind
			indices = indices[findall(x->x!=i,indices)]
		end
		return indices
	end

    out_m, ind = add_random_patch(sqr_matrix0)
	indices = get_unchanged_indices(sqr_matrix0,ind)
	@test size(ind) == (1,2)
	@test sqr_matrix0[indices] == out_m[indices]

	big_sqr_matrix0 = sqr_matrix0 .*100
	out_m, ind = add_random_patch(big_sqr_matrix0, patch_size=1,total_patches=2)
	indices = get_unchanged_indices(big_sqr_matrix0,ind)
	@test size(ind) == (2,2)
	@test big_sqr_matrix0[indices] == out_m[indices]
	@test sum(big_sqr_matrix0[ind] .!= out_m[ind]) == 4
	@test sum(big_sqr_matrix0[ind] .== out_m[ind]) == 0

	out_m, ind = add_random_patch(sqr_matrix1, patch_size=1,total_patches=2)
	indices = get_unchanged_indices(sqr_matrix1,ind)
	@test size(ind) == (2,2)
	@test sqr_matrix1[indices] == out_m[indices]

	# TODO those 2 tests fails when random value is the equal to one that is replaced
	# @test sum(sqr_matrix1[ind] .!= out_m[ind]) == 4
	# @test sum(sqr_matrix1[ind] .== out_m[ind]) == 0

	# ===
	input_matrix = sqr_matrix1
	function test_adding_rand_patch(input_matrix, t_patches,p_size)
		out_m, ind = add_random_patch(input_matrix, patch_size=p_size, total_patches=t_patches)
		indices = get_unchanged_indices(input_matrix,ind)
		@test size(ind) == (t_patches*p_size^2,2)
		@test input_matrix[indices] == out_m[indices]
		# For values from range, tests below does not make sense:
		# @test sum(input_matrix[ind] .!= out_m[ind]) == length(ind)
		# @test sum(input_matrix[ind] .== out_m[ind]) == 0
	end
	t_patches = 1
	p_size = 2
		test_adding_rand_patch(sqr_matrix0, t_patches,p_size)
		test_adding_rand_patch(sqr_matrix1, t_patches,p_size)
		test_adding_rand_patch(sqr_matrix2, t_patches,p_size)
		test_adding_rand_patch(sqr_matrix3, t_patches,p_size)

	t_patches = 1
	p_size = 3
		test_adding_rand_patch(sqr_matrix0, t_patches,p_size)
		test_adding_rand_patch(sqr_matrix1, t_patches,p_size)
		test_adding_rand_patch(sqr_matrix2, t_patches,p_size)
		test_adding_rand_patch(sqr_matrix3, t_patches,p_size)

	t_patches = 1
	p_size = 4
	correct_error = 3
    # TODO change this into test_throws
	try
		add_random_patch(sqr_matrix0, patch_size=p_size, total_patches=t_patches)
	catch err
		my_error = 0
		if isa(err, DomainError)
			println("DomainError")
			my_error = 1
		elseif isa(err, DimensionMismatch)
			println("DimensionMismatch")
			my_error = 2
		else
			println("Unknow error")
			my_error = 3
		end
		@test my_error == correct_error
	end

	test_adding_rand_patch(sqr_matrix1, t_patches, p_size)
	test_adding_rand_patch(sqr_matrix2, t_patches, p_size)
	test_adding_rand_patch(sqr_matrix3, t_patches, p_size)

	t_patches = 3
	p_size = 5
	test_adding_rand_patch(sqr_matrix2, t_patches,p_size)
	test_adding_rand_patch(sqr_matrix3, t_patches,p_size)

	# ===
	# locations
	locations1 = [CartesianIndex(1,2), CartesianIndex(2,3)]
	locations2 = [CartesianIndex(1,2), CartesianIndex(99,3)]


	t_patches = 1
	p_size = 1
	out_m, ind = add_random_patch(sqr_matrix1, patch_size=p_size, total_patches=t_patches,locations=locations1)
		indices = get_unchanged_indices(sqr_matrix1,ind)
		@test ind[:,1] == locations1
		@test size(ind) == (size(locations1)[1]*p_size^2,2)
		@test sqr_matrix1[indices] == out_m[indices]
		# The number of
		@test sum(sqr_matrix1[locations1] .!= out_m[locations1]) == length(locations1)
		@test sum(sqr_matrix1[locations1] .== out_m[locations1]) == 0

		correct_error = 0
		try
			out_m, ind = add_random_patch(sqr_matrix1, patch_size=p_size, total_patches=t_patches,locations=locations2)
		catch err
			# global correct_error
			if isa(err, DomainError)
				correct_error = 1
			else
				correct_error = 2
			end
		finally
			# global correct_error
			@test correct_error == 2
		end
	# TODO test for index below diagonal
	# TODO too many indices
end
