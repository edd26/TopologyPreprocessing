using TopologyPreprocessing
using Test
# using Eirene

#%%
@testset "BettiCurves.jl" begin
    sample_distance_matrix1 = [0  1  25 4  5  9  13 17;
                                  1  0  2  26 6  10 14 18;
                                  25 2  0  3  7  11 15 19;
                                  4  26 3  0  8  12 16 20;
                                  5  6  7  8  0  21 27 24;
                                  9  10 11 12 21 0  22 28;
                                  13 14 15 16 27 22 0  23;
                                  17 18 19 20 24 28 23 0 ]
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
    # get_bettis
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

    # normalise_bettis
    #  as betticurve results
    for matrix = [sample_distance_matrix1, sample_distance_matrix2]
        for max_B_dim = 1:4, min_B_dim = 1:3
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

    # get_vectorized_bettis
    let max_B_dim = 5,
        min_B_dim = 1,
        eirene_results = eirene(sample_distance_matrix1, model="vr", maxdim = max_B_dim)

        let eirene_bettis = get_bettis(eirene_results, max_B_dim, min_dim=min_B_dim),
            vectorized_bettis = get_vectorized_bettis(eirene_results, max_B_dim, min_dim=min_B_dim)

            @test size(vectorized_bettis)[2] == max_B_dim - (min_B_dim-1)
            for d in min_B_dim:max_B_dim
                @test vectorized_bettis[:,d] == eirene_bettis[d][:,2]
            end
        end
    end
end


@testset "BettiCurves.jl -> plot bettis" begin
    # TODO remove tests which test Plots.plot function and not plot_bettis functionality
    sample_distance_matrix1 = [0  1  25 4  5  9  13 17;
                                  1  0  2  26 6  10 14 18;
                                  25 2  0  3  7  11 15 19;
                                  4  26 3  0  8  12 16 20;
                                  5  6  7  8  0  21 27 24;
                                  9  10 11 12 21 0  22 28;
                                  13 14 15 16 27 22 0  23;
                                  17 18 19 20 24 28 23 0 ]
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

    # plot_bettis tests for get_bettis:
    let max_B_dim = 5,
        min_B_dim = 1,
        eirene_results = eirene(sample_distance_matrix1, model="vr", maxdim = max_B_dim)
        all_bettis = get_bettis(eirene_results, max_B_dim)

        p = plot_bettis(all_bettis);
        @test length(p.series_list) == max_B_dim-(min_B_dim-1)
        @test p.attr[:plot_title] == ""

        @test_throws DomainError plot_bettis(all_bettis, min_dim=max_B_dim+1)

        for (dim_index, dim)= enumerate(min_B_dim:max_B_dim)
            @test p.series_list[dim_index][:label] == "β$(dim)"
            if !isnan(all_bettis[dim_index][:,1][1])
                @test p.series_list[dim_index][:x] == all_bettis[dim_index][:,1]
                @test p.series_list[dim_index][:y] == all_bettis[dim_index][:,2]
            end
        end


        # for dim = min_B_dim:max_B_dim
        #     p = plot_bettis(all_bettis, min_dim = dim);
        #     @test length(p.series_list) == max_B_dim-(dim-1)
        # end


        let p1 = plot_bettis(all_bettis, betti_labels=false)
            for (dim_index, dim)= enumerate(min_B_dim:max_B_dim)
                @test p1.series_list[dim][:label] == "y$(dim)"
                if !isnan(all_bettis[dim_index][:,1][1])
                    @test p1.series_list[dim][:x] == all_bettis[dim][:,1]
                    @test p1.series_list[dim][:y] == all_bettis[dim][:,2]
                end
            end
        end

        let lw=4,
            p1 = plot_bettis(all_bettis, betti_labels=true, lw=lw)
            for (dim_index, dim)= enumerate(min_B_dim:max_B_dim)
                @test p1.series_list[dim_index][:label] == "β$(dim)"
                @test p1.series_list[dim_index][:linewidth] == lw
            end
        end

        let plt_title = "test_title",
            p1 = plot_bettis(all_bettis, title=plt_title, lw=9, xlabel="2")
            @test_skip p1.attr[:plot_title] == plt_title # why plot-title is not returning the title?
            for dim = min_B_dim:max_B_dim
                @test p1.series_list[dim][:label] == "β$(dim)"
            end
        end

        let plt_title = "test_title",
            lw = 9,
            p1 =  plot_bettis(all_bettis, title=plt_title, lw=lw, xlabel="2", default_labels=false)
            @test_skip p1.attr[:plot_title] == plt_title # why plot-title is not returning the title?
            for dim = min_B_dim:max_B_dim
                @test p1.series_list[dim][:label] == "β$(dim)"
                @test p1.series_list[dim][:linewidth] == lw
                # @test for xlabel
                # @test for no  label
            end
        end

    end
    # plot_bettis tests for get_vectorized_bettis:
    let max_B_dim = 5,
        min_B_dim = 1,
        eirene_results = eirene(sample_distance_matrix1, model="vr", maxdim = max_B_dim)
        all_bettis = get_vectorized_bettis(eirene_results, max_B_dim)

    end
end


@testset "BettiCurves.jl -> area under betti curves" begin

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
                                     33 34 35 36 37 38 39 40 45 1;],
            max_B_dim = 5,
            min_B_dim = 1
    #==
    checks if the size is anyhow changed during proces;
    checks is the values Array{Matrix} and reshaped matrix are the same

    ==#
        for min_B_dim in [1, 2, 3, 4, 5]
            eirene_results1 =
                eirene(sample_distance_matrix1, model = "vr", maxdim = max_B_dim)
            eirene_results2 =
                eirene(sample_distance_matrix2, model = "vr", maxdim = max_B_dim)

            bettis_collection = [
                get_bettis(eirene_results1, max_B_dim),
                get_bettis(eirene_results2, max_B_dim),
            ]

            for bettis_col in bettis_collection
                total_vecs = length(bettis_col)
                vec_len, vec_width = size(bettis_col[1])

                reshaped_betti = TopologyPreprocessing.vectorize_bettis(bettis_col)
                @test vec_len .== size(reshaped_betti, 1)
                @test total_vecs .== size(reshaped_betti, 2)
                for k = 1:total_vecs
                    @test reshaped_betti[:, k] == bettis_col[k][:, 2]
                end
            end
        end
        #==
        checks if get vectorized bettis has same values as get_bettis
        ==#
        for min_B_dim in [1, 2, 3, 4, 5]
            eirene_results1 =
                eirene(sample_distance_matrix1, model = "vr", maxdim = max_B_dim)
            eirene_results2 =
                eirene(sample_distance_matrix2, model = "vr", maxdim = max_B_dim)

            bettis_collection = [
                get_bettis(eirene_results1, max_B_dim),
                get_bettis(eirene_results2, max_B_dim),
            ]
            vec_bett_collection = [
                get_vectorized_bettis(eirene_results1, max_B_dim),
                get_vectorized_bettis(eirene_results2, max_B_dim),
            ]

            for index = 1:length(bettis_collection)
                bettis_col = bettis_collection[index]
                vec_bettis_col = vec_bett_collection[index]

                total_vecs = length(bettis_col)
                for k = 1:total_vecs
                    @test vec_bettis_col[:, k] == bettis_col[k][:, 2]
                end
            end
        end
    end
end
