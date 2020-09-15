using TopologyPreprocessing
using Test

# using MatrixOrganization

@testset "BettiCurves.jl" begin
    # get_bettis
    let sample_distance_matrix1 = [0  1  25 4  5  9  13 17;
                                   1  0  2  26 6  10 14 18;
                                   25 2  0  3  7  11 15 19;
                                   4  26 3  0  8  12 16 20;
                                   5  6  7  8  0  21 27 24;
                                   9  10 11 12 21 0  22 28;
                                   13 14 15 16 27 22 0  23;
                                   17 18 19 20 24 28 23 0 ],
        sample_distance_matrix2 = [1  1  41 4  5  9  13 17 25 33;
                                  1  1  2  42 6  10 14 18 26 34;
                                  41 2  1  3  7  11 15 19 27 35;
                                  4  42 3  1  8  12 16 20 28 36;
                                  5  6  7  8  1  21 43 24 29 37;
                                  9  10 11 12 21 1  22 44 30 38;
                                  13 14 15 16 43 22 1  23 31 39;
                                  17 18 19 20 24 44 23 1  32 40;
                                  25 26 27 28 29 30 31 32 1  45;
                                  33 34 35 36 37 38 39 40 45 1;]
        for matrix = [sample_distance_matrix1, sample_distance_matrix2]
            for max_B_dim = 1:4
                eirene_results = eirene(matrix, model="vr", maxdim = max_B_dim)
                all_bettis = get_bettis(eirene_results, max_B_dim)

                @test length(all_bettis) == max_B_dim
                @test all_bettis isa Vector{Array{Float64,2}}
            end

            for max_B_dim = 1:4, min_B_dim = 1:3
                if min_B_dim > max_B_dim
                    @debug "Continue at " min_B_dim, max_B_dim
                    continue
                end
               eirene_results = eirene(matrix, model="vr", maxdim = max_B_dim)
               all_bettis = get_bettis(eirene_results, max_B_dim, min_dim=min_B_dim)

               @test length(all_bettis) == max_B_dim - (min_B_dim-1)
               @test all_bettis isa Vector{Array{Float64,2}}
            end
        end
    end

    # normalise_bettis
    let sample_distance_matrix1 = [0  1  25 4  5  9  13 17;
                                   1  0  2  26 6  10 14 18;
                                   25 2  0  3  7  11 15 19;
                                   4  26 3  0  8  12 16 20;
                                   5  6  7  8  0  21 27 24;
                                   9  10 11 12 21 0  22 28;
                                   13 14 15 16 27 22 0  23;
                                   17 18 19 20 24 28 23 0 ],
        sample_distance_matrix2 = [1  1  41 4  5  9  13 17 25 33;
                                  1  1  2  42 6  10 14 18 26 34;
                                  41 2  1  3  7  11 15 19 27 35;
                                  4  42 3  1  8  12 16 20 28 36;
                                  5  6  7  8  1  21 43 24 29 37;
                                  9  10 11 12 21 1  22 44 30 38;
                                  13 14 15 16 43 22 1  23 31 39;
                                  17 18 19 20 24 44 23 1  32 40;
                                  25 26 27 28 29 30 31 32 1  45;
                                  33 34 35 36 37 38 39 40 45 1;]

        #  as betticurve results
        for matrix = [sample_distance_matrix1, sample_distance_matrix2]
            for max_B_dim = 1:4
                if min_B_dim > max_B_dim
                    @debug "Continue at " min_B_dim, max_B_dim
                    continue
                end
                eirene_results = eirene(matrix, model="vr", maxdim = max_B_dim)
                betti_result = betticurve(eirene_results, dim=max_B_dim)

                normed_all_bettis = normalise_bettis(betti_result)
                @test typeof(normed_all_bettis) == typeof(betti_result)
                @test length(normed_all_bettis) != max_B_dim
                @test size(normed_all_bettis) == size(betti_result)
                @test normed_all_bettis isa Array{Float64,2}

                # Betti values are unchanged:
                @test normed_all_bettis[:,2] == betti_result[:,2]
                # Max val is 1
                @test findmax(normed_all_bettis[:,1])[1] == 1.
            end
        end

        #  as get_bettis results
        for matrix = [sample_distance_matrix1, sample_distance_matrix2]
            for max_B_dim = 1:4, min_B_dim = 1:3
                if min_B_dim > max_B_dim
                    @debug "Continue at " min_B_dim, max_B_dim
                    continue
                end
                eirene_results = eirene(matrix, model="vr", maxdim = max_B_dim)
                all_bettis = get_bettis(eirene_results, max_B_dim)

                normed_all_bettis = normalise_bettis(all_bettis)
                @test typeof(normed_all_bettis) == typeof(all_bettis)
                @test length(normed_all_bettis) == max_B_dim
                @test normed_all_bettis isa Vector{Array{Float64,2}}

                # Betti values are unchanged:
                @test normed_all_bettis[max_B_dim][:,2] == all_bettis[max_B_dim][:,2]
                # Max val is 1
                @test findmax(normed_all_bettis[max_B_dim][:,1])[1] == 1.
            end
        end
    end

end
