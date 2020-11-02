# TopologyPreprocessing

[![Build Status](https://travis-ci.com/edd26/TopologyPreprocessing.jl.svg?branch=master)](https://travis-ci.com/edd26/TopologyPreprocessing.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/edd26/TopologyPreprocessing.jl?svg=true)](https://ci.appveyor.com/project/edd26/TopologyPreprocessing-jl)
[![Codecov](https://codecov.io/gh/edd26/TopologyPreprocessing.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/edd26/TopologyPreprocessing.jl)
[![Coveralls](https://coveralls.io/repos/github/edd26/TopologyPreprocessing.jl/badge.svg?branch=master)](https://coveralls.io/github/edd26/TopologyPreprocessing.jl?branch=master)


## Installation

At the moment, this package is not added to the Julia Package registry, so the
installation must be done manually. This can be done in 2 ways:

1. Modifying Julia ``LOAD_PATH``:
    1. Go to your destination folder.
    2. Download the repository:

        git pull https://github.com/edd26/TopologyPreprocessing.jl.git

    3. Open Julia startup file, the default would be:

        ~/.julia/config/startup.jl

    4. Add the following line:

    push!(LOAD_PATH, "repo_download_path/TopologyPreprocessing/src/")

    where `repo_download_path` is the path to where repository was pulled.

    5. Restart Julia for the change to work.
2. Cloning the repo into Julia packages.
    1. navigate to Julia packages folder, the default is:

        ~/.julia/packages
    2. Download the repository:

        git pull https://github.com/edd26/TopologyPreprocessing.jl.git

