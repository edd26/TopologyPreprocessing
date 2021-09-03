#!/bin/bash

git clone --depth=50 --branch=master https://github.com/edd26/TopologyPreprocessing.jl.git 
# git clone --depth=50 --branch=master https://github.com/JuliaImages/Images.jl.git

cd TopologyPreprocessing.jl
# cd Images.jl

export JULIA_PROJECT=@.

export JL_PKG=Images

git fetch --unshallow

julia --project=.ci/ -t 1 --trace-compile=stdout -e 'using Pkg; Pkg.instantiate(); Pkg.build(verbose=true)'
julia --check-bounds=yes -t 1 --trace-compile=stdout --color=yes -e "ENV[\"JULIA_DEBUG\"] = \"all\"; using Pkg; Pkg.test(coverage=true);"
