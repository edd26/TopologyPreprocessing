"""
    check_existing_dir(dir_path::String)

Checks if the directory under `dir_path` exists. If not, throws IOError
"""
function check_existing_dir(dir_path::String)
    if !isdir(dir_path)
        @warn "Folder $(data_path) does not exists in current directory."
        @warn "Terminating execution."
        throw(ErrorException("Can nor find folder \""*dir_path*"\"."))
    end
end
