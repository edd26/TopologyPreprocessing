# import Makie
import VideoIO
using StatsBase
using Images
using ImageFeatures
# using TestImages
# using ImageDraw
using CoordinateTransformations
# using Makie
# using Logging

export get_video_array_from_file,
        get_video_dimension,
        get_video_mask,
        extract_pixels_from_video,
        vectorize_video,
        get_pairwise_correlation_matrix,
        get_average_from_tiles,
        rotate_img_around_center,
        rotate_vid_around_center,
        export_images_to_vid,
        rotate_and_save_video,
        get_local_correlations,
        get_local_centers,
        get_local_gradients,
        normalize_to_01,
        shift_to_non_negative,
        plotimg;

"""
 get_video_array_from_file(video_name)

Returns array to which video frames are copied. Frames are in grayscale.

Function opens stream, then loads the file and gets all the frames from a
video.
 """
function get_video_array_from_file(video_name)
   video_streamer = VideoIO.open(video_name) # candle
   video_array = Vector{Array{UInt8}}(undef,0);
   video_file = VideoIO.openvideo(video_streamer, target_format=VideoIO.AV_PIX_FMT_GRAY8)

   while !eof(video_file)
      push!(video_array,reinterpret(UInt8, read(video_file)))
   end
   close(video_file)

   return video_array
end

"""
    get_video_dimension(video_array)

Returns the tuple which contains width, height and the number of the frames of
array in whcih video was loaded.
"""
function get_video_dimension(video_array)
   v_hei = size(video_array[1],1)
   v_wid = size(video_array[1],2)
   v_len = size(video_array,1)

   video_dim_tuple = (video_height=v_hei, video_width=v_wid, video_length=v_len)

   return video_dim_tuple
end

"""
    get_video_mask(points_per_dim, video_dimensions; distribution="uniform", sorted=true, x=1, y=1)

Returns matrix of size @points_per_dim x 2 in which indicies of video frame are
stored. The indicies are chosen based one the @distribution argument. One option
is uniform distribution, the second is random distribution.

Uniform distribution: distance between the points in given dimension is the
 even, but vertical distance may be different from horizontal distance between points. This depends on the size of a frame in a image.

Random distribution: the distance between the points is not constant, because
the points are chosen randomly in the ranges 1:horizontal size of frame,
1:vertical size of frame. The returned values may be sorted in ascending order,
if @sorted=true.
"""
function get_video_mask(points_per_dim, video_dimensions;
                            distribution="uniform", sorted=true, patch_params)
    video_height, video_width,  = video_dimensions
    y=patch_params["y"]
    spread=patch_params["spread"]

    if x == 1
        x = Int64(floor(video_width/2))
        @warn "Given x is to close to the border. Seeting the value to " x
    elseif x < Int64(ceil(points_per_dim/2))
        x = Int64(ceil(points_per_dim/2))
        @warn "Given x is to close to the border. Seeting the value to " x
    elseif x > video_width-Int64(ceil(points_per_dim/2))
        x = video_width - Int64(ceil(points_per_dim/2))
        @warn "Given x is to close to the border. Seeting the value to " x
    end

    if y == 1
        y = Int64(floor(video_height/2))
        @warn "Given y is to close to the border. Seeting the value to " y
    elseif y < Int64(ceil(points_per_dim/2))
        y = Int64(ceil(points_per_dim/2))
        @warn "Given y is to close to the border. Seeting the value to " y
    elseif y > video_height-Int64(ceil(points_per_dim/2))
        y = video_height - Int64(ceil(points_per_dim/2))
        @warn "Given y is to close to the border. Seeting the value to " y
    end

    if spread*points_per_dim+x > video_width || spread*points_per_dim+y > video_height
        @warn "Given patch parameters might result in indicies exceeding frame size."
    end

    if distribution == "uniform"
        columns = points_per_dim
        rows = points_per_dim

        # +1 is used so that the number of points returned is as requested
        row_step = Int64(floor(video_height/rows))
        column_step = Int64(floor(video_width/columns))

        (video_height/row_step != points_per_dim) ? row_step+=1 : row_step
        (video_width/column_step !=
                                points_per_dim) ? column_step+=1 : video_width

        vertical_indicies = collect(1:row_step:video_height)
        horizontal_indicies = collect(1:column_step:video_width)

        vertical_indicies = reshape(vertical_indicies, (1,points_per_dim))
        horizontal_indicies = reshape(horizontal_indicies, (1,points_per_dim))

        indicies_set = [vertical_indicies; horizontal_indicies]
    elseif distribution == "random"
        vertical_indicies = rand(1:video_height,1, points_per_dim)
        horizontal_indicies = rand(1:video_width,1, points_per_dim)

        if sorted
            vertical_indicies = sort(vertical_indicies[1,:])
            horizontal_indicies = sort(horizontal_indicies[1,:])

            vertical_indicies = reshape(vertical_indicies, (1,points_per_dim))
            horizontal_indicies =
                              reshape(horizontal_indicies, (1,points_per_dim))
        end
        indicies_set = [vertical_indicies; horizontal_indicies]
    elseif distribution == "patch"
        indicies_set = [collect(1:spread:(spread*points_per_dim)).+x collect(1:spread:(spread*points_per_dim)).+y]'
    end

   return indicies_set
end

"""
    extract_pixels_from_video(video_array, indicies_set, video_dim_tuple)

Takes every frame from @video_array and extracts pixels which indicies are in
@indicies_set, thus creating video with only chosen indicies.
"""
function extract_pixels_from_video(video_array, indicies_set, video_dim_tuple)
   rows = size(indicies_set,2)
   columns = size(indicies_set,2)
   video_length = video_dim_tuple[3]

   extracted_pixels = zeros(rows, columns, video_length)
   for frame_number in 1:video_length
      extracted_pixels[:,:,frame_number] =
                video_array[frame_number][indicies_set[1,:],indicies_set[2,:]]
   end

   return extracted_pixels
end


"""
    vectorize_video(video)

Rearrenges the video so that set of n frames (2D matrix varying in
time) the set of vectors is returned, in which each row is a pixel, and each
column is the value of the pixel in n-th frame.
"""
function vectorize_video(video)
    video_length = size(video, 3)
    rows = size(video,1)
    columns = size(video,2)

    number_of_vectors = rows*columns

    vectorized_video = zeros(number_of_vectors, video_length)

    index = 1;
    for row=1:rows
        for column=1:columns
            vectorized_video[index,:] = video[row, column,:];
            index = index+1;
        end
    end

    return vectorized_video
end

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
    get_average_from_tiles(extracted_pixels_matrix, N)

Fnction takes a 3D array in which video is stored and splits every frame into
non overlaping tiles of size NxN. If size of @extracted_pixels_matrix is not
square of N, then only N^2 x N^2 matrix will be used for averaging.
"""
function get_average_from_tiles(extracted_pixels_matrix, N)
    # N = size(extracted_pixels,1)
    num_frames = size(extracted_pixels_matrix,3)
    mask_matrix = ones(N, N)
    result_matrix = zeros(N, N, num_frames)
    col_index = 1
    row_index = 1

    for frame = 1:num_frames
        for col = 1:N:N^2
            for row = 1:N:N^2
                result_matrix[mod(col,N), mod(row,N), frame] =
                        dot(extracted_pixels_matrix[col:(col+N-1),
                            row:(row+N-1), frame], mask_matrix) ./N^2
                row_index += 1
            end
            col_index += 1
        end
    end
    return result_matrix
end


"""
    rotate_img_around_center(img, angle = 5pi/6)

Function rotates a single image (or a frame) around the center of the image by
@angle radians.
"""
function rotate_img_around_center(img, angle = 5pi/6)
  θ = angle
  rot = recenter(RotMatrix(θ), [size(img)...] .÷ 2)  # a rotation around the center
  x_translation = 0
  y_translation = 0
  tform = rot ∘ Translation(y_translation, x_translation)
  img2 = warp(img, rot, axes(img))

  return img2
end



"""
    rotate_vid_around_center(img, rotation = 5pi/6)

Function rotates a video around the center of the image by @rotation radians and
 the outout into matrix.
"""
function rotate_vid_around_center(src_vid_path,src_vid_name; rotation = 5pi/6)
    video_array = []
    video_src_strm = VideoIO.open(src_vid_path*src_vid_name)
    video_src = VideoIO.openvideo(video_src_strm,
                                        target_format=VideoIO.AV_PIX_FMT_GRAY8)

    while !eof(video_src)
      img = read(video_src)
      img = rotate_img_around_center(img, rotation)

      push!(video_array,img)
    end
    close(video_src)

  return video_array
end


"""
    export_images_to_vid(video_array, dest_file)

Exports set of images stored in @video_array to the dest_file.

"""
function export_images_to_vid(video_array, dest_file)
    @debug "Exporting set of images to file"
    fname = dest_file

    video_dimensions = get_video_dimension(video_array)
    h = video_dimensions.video_height
    w = video_dimensions.video_width
    nframes = video_dimensions.video_length
    overwrite=true
    fps=30
    options = ``
    ow = overwrite ? `-y` : `-n`

    open(`ffmpeg
            -loglevel warning
            $ow
            -f rawvideo
            -pix_fmt rgb24
            -s:v $(h)x$(w)
            -r $fps
            -i pipe:0
            $options
            -vf "transpose=0"
            -pix_fmt yuv420p
            $fname`, "w") do out
        for i = 1:nframes
            write(out, convert.(RGB{N0f8}, clamp01.(video_array[i])))
        end
    end
    @debug "Video was saved"
end


"""
    rotate_and_save_video(src_vid_path, src_vid_name, dest_vid_name;
                                                                rotation=5pi/6)

Fuction opens the @src_vid_name file, collects all the frames and then rotates
the frame aroung the center and saves new video as @dest_vid_name at
@src_vid_path.

Function was tested for following extensions;
    .mov

A solution for writing to a video file was taken from:
https://discourse.julialang.org/t/creating-a-video-from-a-stack-of-images/646/7
"""
function rotate_and_save_video(src_vid_path, src_vid_name, dest_vid_name, rotation=5pi/6)
    @debug src_vid_path src_vid_name dest_vid_name

    if !isfile(src_vid_path*src_vid_name)
        @warn "Source file at given path does not exist. Please give another name."
        return
    elseif isfile(src_vid_path*dest_vid_name)
        @warn "File with destination video name at src_video_path already exists. Please give another name."
        return
    end

    video_array = rotate_vid_around_ceter(src_vid_path, src_vid_name)
    @debug "Video was rotated"

    export_images_to_exist_vid(video_array, src_vid_path*dest_vid_name)
    @info "The file was created:\n  $fname"
end




"""
    get_local_correlations(video_array, centers, sub_img_size, shift)

Computes the correlation between the subimages and subimages shifted by values
from range -@shift:@shift and returns array with frames of size
length(@centers) x length(@centers) with the number of frames equal to the
number of rames in @video_array.

Each of the subimage is center around values stored in  @centers
"""
function get_local_correlations(video_array, centers, sub_img_size, shift)
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


function get_local_centers(points_per_dim, video_dimensions, shift=0, sub_img_size=0 )
    /# TODO Applied teproray solution here, so it works only for local gradients
    # start = 0
    # (points_per_dim>shift) ? start_ind = ceil(Int, points_per_dim/2)+ shift :
    #                         start=shift
    start_ind = ceil(Int, sub_img_size/2)
    min_va,  = findmin(video_dimensions)
    last_ind = min_va - start_ind

    set = broadcast(floor, Int, range(start_ind, stop=last_ind,  length=points_per_dim))
    centers = [set set]'
    return centers
end


"""
    get_local_gradients(video_array, centers, sub_img_size)

Computes the gradients in the subimage, takes the mean of sum of absolute
values of both hotizontal and vertical gradients as a representative of a
subimage.
"""
function get_local_gradients(video_array, centers, sub_img_size)
    @debug "Entering get_local_gradients"
    half_size = ceil(Int,(sub_img_size-1)/2)
    half_range = half_size
    h, w, len = get_video_dimension(video_array)
    extracted_pixels = zeros(sub_img_size, sub_img_size, len)

    @debug "starting frame procesing"
    for frame = 1:len
        img = video_array[frame]
        img_grad = imgradients(img, KernelFactors.ando3, "replicate")
        img_grad_abs = map(abs, img_grad[1]) + map(abs, img_grad[2])
        for index_x = 1:size(centers,2)
            c_x = centers[2, index_x]
            for index_y = 1:size(centers,2)
                c_y = centers[1, index_y]
                sub_img  = img_grad_abs[(c_x-half_size):(c_x+half_size),
                                        (c_y-half_size):(c_y+half_size)]

                extracted_pixels[index_x, index_y, frame] =mean(sub_img)
            end
        end
        # @debug "Next frame" frame
    end
    return extracted_pixels
end

"""
    normalize_to_01(matrix, norm_factor=256)

Returns a matrix which values are in range [0, 1]. If the values in the input
matrix are below 0, then they are shifted so that only positive numbers are in
the matrix. If the values in the matrix of shifted matrix excced value of the
@norm_factor parameter, then the matrix is normalized to the maximal value from
the matrix.
"""
function normalize_to_01(matrix, norm_factor=256)
    normalized_matrix = copy(matrix)
    min_val = findmin(normalized_matrix)[1]
    max_val = findmax(normalized_matrix)[1]

    if min_val < 0
        normalized_matrix .-= min_val
    end

    if max_val > norm_factor
        @warn "Values normalized to maximal value, not notmalization factor."
        normalized_matrix = normalized_matrix./max_val
    else
        normalized_matrix = normalized_matrix./norm_factor
    end
    return normalized_matrix
end

"""
    shift_to_non_negative(matrix)

Returns a matrix in which values are non-negative. This is done by finding the
minimal value in the input matrix and adding its absolute value to the matix
elements.
"""
function shift_to_non_negative(matrix)

    min_val = findmin(matrix)[1]
    if min_val < 0
        return matrix .-= min_val
    else
        return matrix
    end
end


"""
    plotimg(matrix_to_plot)

Display an image as a plot. The values from the input matrix are adjusted to the
value range of [0, 1].

If @cut_off is true then the matrix values above 256 are set to 256 and then all
values are normalized to the value 256. If @cut_off is false, then values are
normalized to maximal value.
"""
function plotimg(matrix_to_plot, cut_off=false)
    matrix_type = typeof(matrix_to_plot)
    min_val = findmin(matrix_to_plot)[1]
    int_types_arr = [Matrix{UInt8}; Matrix{UInt16}; Matrix{UInt32};
                    Matrix{UInt64}; Matrix{UInt128}; Matrix{Int8};
                    Matrix{Int16}; Matrix{Int32}; Matrix{Int64};
                    Matrix{Int128}]
    float_types_arr = [Matrix{Float16} Matrix{Float32} Matrix{Float64}]

    if min_val<0
        matrix_to_plot = shift_to_non_negative(matrix_to_plot)
    end

    max_val = findmax(matrix_to_plot)[1]

    if max_val > 256 && cut_off
        matrix_to_plot[findall(x -> x>256, matrix_to_plot)] = 256
    end

    if in(matrix_type, int_types_arr)
        matrix_to_plot = normalize_to_01(matrix_to_plot)
    elseif in(matrix_type, float_types_arr)
        matrix_to_plot = normalize_to_01(matrix_to_plot, max_val)
    end

    return colorview(Gray, matrix_to_plot)
end
