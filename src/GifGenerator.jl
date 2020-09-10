using Luxor, Colors, FileIO

# export generate_gif,
#         dotted_plane,
#         backdrop,
#         load_gif_to_array;

"""
Generate a gif at @full_path using the function @frame_function. Image size and
 length sotred in the tuple @gif_dims.
"""
function generate_gif(frame_function, full_path, gif_dims)
    w, h, len = gif_dims
    len -= 1
    demo = Movie(w,h, "test")

    animate(demo, [
        Scene(demo, backdrop, 0:len),
        Scene(demo, frame_function, 0:len, easingfunction=easeinoutcubic)
                    ], creategif=true, pathname=full_path,  usenewffmpeg=false)
    @info "Generation finished."
end


"""
Generate a frame filled whith dots moving toward or out form the middle.
"""
function dotted_plane(scene, framenumber)
    distance = mod(framenumber,50)
    if distance <= 25
        for radius = distance:5:800
            setdash("dot")
            sethue("gray30")
            A, B = [Point(x, 0) for x in [-radius, radius]]
            circle(O, radius, :stroke)
        end
    else
        for radius = (50-distance):5:800
            setdash("dot")
            sethue("gray30")
            A, B = [Point(x, 0) for x in [-radius, radius]]
            circle(O, radius, :stroke)
        end
    end
    return
end

"""
Function for background.
"""
function backdrop(scene, framenumber)
    background("white")
end


"""
Load the gif to the array.

The load function is he FileIO function, but its possibilities are enchanced by
    the ImageMagic package. ImageMagic should not be loaded to project- it is
    handled by FileIO.
"""
function load_gif_to_array(fullpath)
    return img = load(fullpath)
end
