#!/bin/bash

git clone --depth=50 --branch=master https://github.com/edd26/TopologyPreprocessing.jl.git 

cd TopologyPreprocessing.jl

export JULIA_PROJECT=@.

julia --color=yes -e "VERSION >= v\"0.7.0-DEV.3630\" && using InteractiveUtils; versioninfo()"

export JL_PKG=TopologyPreprocessing

git fetch --unshallow

julia --color=yes -e "ENV[\"JULIA_DEBUG\"] = \"all\"; using Pkg; Pkg.build(verbose=true);"

julia --check-bounds=yes --color=yes -e "ENV[\"JULIA_DEBUG\"] = \"all\"; using Pkg; Pkg.test(coverage=true);"
