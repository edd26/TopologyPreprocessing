
# ============================================================================
# exported from MatrixProcessing on 10.09.2020
"""
    get_pairwise_correlation_matrix(vectorized_video, tau_max=25)

Computes pairwise correlation of the input signals accordingly to the formula
presented in paper "Clique topology reveals intrinsic geometric structure
in neural correlations" by Chad Giusti et al.

The Computations are done only for upper half of the matrix, the lower half is
a copy of upper half. Computation-wise the difference is at level of 1e-16, but
this causes that inverse is not the same as non-inverse matrix.

"""
function get_pairwise_correlation_matrix(vectorized_video, tau_max=25)
    number_of_signals = size(vectorized_video,1)
    T = size(vectorized_video,2)

    C_ij = zeros(number_of_signals,number_of_signals);
    # log_C_ij = zeros(number_of_signals,number_of_signals);

     # this is given in frames
    lags = -tau_max:1:tau_max


    for row=1:number_of_signals
        for column=row:number_of_signals
            signal_ij = vectorized_video[row,:];
            signal_ji = vectorized_video[column,:];

            # cross_corelation
            ccg_ij = crosscov(signal_ij, signal_ji, lags);
            ccg_ij = ccg_ij ./ T;

            A = sum(ccg_ij[tau_max+1:end]);
            B = sum(ccg_ij[1:tau_max+1]);
            r_i_r_j = 1;
            C_ij[row, column] = max(A, B)/(tau_max*r_i_r_j);
            C_ij[column, row] = C_ij[row, column]
            # log_C_i_j[row, column] = log10(abs(C_ij[row, column]));
        end
    end

    return C_ij
end




"""
    get_subimg_correlations(video_array, centers, sub_img_size, shift)

Computes the correlation between the subimages and subimages shifted by values
from range -@shift:@shift and returns array with frames of size
length(@centers) x length(@centers) with the number of frames equal to the
number of rames in @video_array.

Each of the subimage is center around values stored in  @centers
"""
# TODO Check if this is the same as some of the functions from the ImageProcessing
function get_subimg_correlations(video_array, centers, sub_img_size, shift)
    half_size = ceil(Int,(sub_img_size-1)/2)
    half_range = half_size + shift
    h, w, len = get_video_dimension(video_array)
    extracted_pixels = zeros(sub_img_size, sub_img_size, len)

    for frame = 1:len
        img = video_array[frame]
        for index_x = 1:size(centers,2)
            c_x = centers[2, index_x]
            for index_y = 1:size(centers,2)
                c_y = centers[1, index_y]
                subimage = img[(c_x-half_range):(c_x+half_range),
                                (c_y-half_range):(c_y+half_range)]
                center = img[(c_x-half_size):(c_x+half_size), (c_y-half_size):(c_y+half_size)]

                for left_boundary = 1:(2*shift+1)
                    for lower_boundary = 1:(2*shift+1)
                        corelation = center .* subimage[left_boundary:left_boundary+sub_img_size-1, lower_boundary:lower_boundary+sub_img_size-1]
                        corelation = sum(corelation)
                        extracted_pixels[index_x, index_y, frame] += corelation
                    end
                end
                extracted_pixels[index_x, index_y, frame] /= 256*(sub_img_size^2)*(shift*2)^2
            end
        end
    end
    return extracted_pixels
end
