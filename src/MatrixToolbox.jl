using Distances
using DataFrames
using Random
using LightGraphs
using DelimitedFiles

export generate_random_point_cloud,
        generate_geometric_matrix,
        generate_shuffled_matrix,
        generate_random_matrix,
        generate_matrix_ordering,
        generate_set_of_graphs,
        plot_betti_numbers,
        save_matrix_to_file;


"""
Returns a random matrix of size @number_of_points x @dimensions in which every
column is a point and every n-th row is a position in the n-th dimension.
"""
generate_random_point_cloud(number_of_points = 12, dimensions=2) =
                                    rand(Float64, dimensions, number_of_points)


"""
Return a matrix which stores the pariwise distances between every point in the
@random_points matrix.
"""
function generate_geometric_matrix(random_points)
    geometric_matrix = pairwise(Euclidean(), random_points, dims=2)
    return geometric_matrix
end

"""
Returns a symetric matrix with randomly permuted valuse from the @input_matrix.
"""
function generate_shuffled_matrix(input_matrix)
    matrix_size = size(input_matrix,1)

    indicies_collection = findall(x->x>0, input_matrix)
    rand!(indicies_collection, indicies_collection)
    shuffeled_matrix = copy(input_matrix)

    # Swap the elements
    n=1
    for k in 1:matrix_size
        for m in k+1:matrix_size
            a = indicies_collection[n][1]
            b = indicies_collection[n][2]
            shuffeled_matrix[k,m] = input_matrix[a,b]
            shuffeled_matrix[m,k] = input_matrix[b,a]

            shuffeled_matrix[a,b] = input_matrix[k,m]
            shuffeled_matrix[b,a] = input_matrix[m,k]

            n +=1
        end
    end
    return shuffeled_matrix
end

"""
Returns matrix with random values which are symmetric accros the diagonal. The
matrix has @matrix_size rows and @matrix_size columns.
"""
function generate_random_matrix(matrix_size)
    elemnts_above_diagonal = Int((matrix_size^2-matrix_size)/2)
    random_matrix = zeros(matrix_size, matrix_size)
    set_of_random_numbers = rand(elemnts_above_diagonal)

    h = 1
    for k in 1:matrix_size
        for m in k+1:matrix_size
            random_matrix[k,m] = set_of_random_numbers[h]
            random_matrix[m,k] = set_of_random_numbers[h]
            h += 1
        end
    end

    return random_matrix
end


"""
Returns set of graphs generated from the @matrix_ordering. In every succesive
graph, single connection between points is added.

NOTE: the function does not take the coordinates of the numbered vertices.
"""
function generate_set_of_graphs(matrix_size, matrix_ordering)
    vetrices = matrix_size
    edges = matrix_ordering
    num_of_edges = size(edges)[2]

    set_of_graphs = [a=Graph(vetrices) for a=1:num_of_edges]
    edges_counter = zeros(Int, num_of_edges)
    edge_density =  zeros(num_of_edges)

    k=1
    for k in range(1,stop=num_of_edges)~
        add_edge!(set_of_graphs[k], edges[1,k], edges[2,k]);
        edges_counter[k] = ne(set_of_graphs[k])
        edge_density[k] = edges_counter[k]/binomial(matrix_size,2)

        if k<num_of_edges # if is used to eliminate copying at last iteration
            set_of_graphs[k+1] = copy(set_of_graphs[k])
        end
    end
    return set_of_graphs, edge_density
end

"""
Plots Betti curves. The betti numbers should be obtained with the clique-top
library.
"""
function plot_betti_numbers(betti_numbers, edge_density, title="Geometric  matrix"; stop=0.6)
    p1 = plot(edge_density, betti_numbers[:,1], label="beta_0", title=title, legend=:topleft) #, ylims = (0,maxy)
    plot!(edge_density, betti_numbers[:,2], label="beta_1")
    if size(betti_numbers,2)>2
        plot!(edge_density, betti_numbers[:,3], label="beta_2")
    end

    return p1
end

"""
Saves given @matrix to the csv file with the name @filename. If there is no path
added to the @filename, then file saved is in local folder.
"""
function save_matrix_to_file(matrix, filename)
    open(filename, "w") do io
        writedlm(io,  matrix, ',')
    end
end


# =====
# Copied form Julia learning repo

"""
Returns ordering of the @geometric_matrix given as an input. If value @ascending
is set to true, the values are number from the lowest value, to the highest. If
false, the values are numbered from highest to the lowest.
"""
function generate_matrix_ordering(geometric_matrix, ascending = true)
    matrix_size = size(geometric_matrix, 2)
    elemnts_above_diagonal = Int((matrix_size^2-matrix_size)/2)
    matrix_ordering = zeros(Int, 2,elemnts_above_diagonal)

    A = copy(geometric_matrix)

    (ascending) ? (method=findmax) : (method=findmin)

    for element in 1:elemnts_above_diagonal
    #     Find maximal distance
        minimal_value = method(A)
    #     Get the coordinates (only 2 dimensions, because it is distance matrix)
        matrix_ordering[1,element] = Int(minimal_value[2][1])
        matrix_ordering[2,element] = Int(minimal_value[2][2])
    #
    # #     Zero minval in A (above and below diagonal) so next minval can be found
        A[matrix_ordering[1,element], matrix_ordering[2,element]] = 0.0
        A[matrix_ordering[2,element], matrix_ordering[1,element]] = 0.0
    end

    # change from min to max order to the max to min order (? necessary ?)
    if ascending
        matrix_ordering = matrix_ordering[:,end:-1:1]
    end

    return matrix_ordering
end

"""
    function get_geometric_matrix(points, dimensions; save_as_file=false)

Created a point cloud with 'points' number of points from 'dimension'
dimensional eucidean unit cube and computes distances between the points.

Distance matrix may be saved to csv file by setting 'save_as_file' to 'true'.
"""
function get_geometric_matrix(points, dimensions; save_as_file=false)
    point_cloud = generate_random_point_cloud(points,dimensions)
    geom_mat = generate_geometric_matrix(point_cloud)

    if save_as_file
        open("geometric_matrix_points$(points)_dims$(dimensions).csv", "w") do io
           writedlm(io, geom_mat, ',')
       end
   end
   return geom_mat
end
