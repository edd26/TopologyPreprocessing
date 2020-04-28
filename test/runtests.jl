using TopologyPreprocessing
using Test

@testset "TopologyPreprocessing.jl" begin
    # Write your own tests here.
end

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
	@test matrix_poling(in_vector) == [4 4 4 4 4 4 4]
	# @test matrix_poling!(in_vector) === in_vector
	# @test matrix_poling!(in_vector) == [4 4 4 4 4 4 4]

	@test matrix_poling(in_matrix) !== in_matrix
	@test matrix_poling(in_matrix) == 9 .* ones(size(in_matrix))
	# @test matrix_poling!(in_matrix) === in_matrix
	# @test matrix_poling!(in_matrix) == 9 .* ones(size(in_matrix))

	@test matrix_poling(in_matrix2[1:2,1:2]) == 6 .* ones(size(in_matrix2[1:2,1:2]))
	@test matrix_poling(in_matrix2[1:2,1:2]) != in_matrix2[1:2,1:2]

	# ====
	# Subsampling matrix

	# Function is supposed to work only on upper half, and here the upper half is too small, so there are no operations
	@test subsample_matrix(sqr_matrix0, subsamp_size=2, method="max_pooling") == resulting_matrix0_2

	@test subsample_matrix(sqr_matrix1, subsamp_size=2, method="max_pooling") == resulting_matrix1_2

	@test subsample_matrix(sqr_matrix2, subsamp_size=2, method="max_pooling") == result_matrix2_2
	@test subsample_matrix(sqr_matrix2, subsamp_size=3, method="max_pooling") == result_matrix2_3
	@test subsample_matrix(sqr_matrix2, subsamp_size=4, method="max_pooling") == result_matrix2_4
end
