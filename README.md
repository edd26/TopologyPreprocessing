# TopologyPreprocessing

 ![Continuous Integration (CI) status](https://github.com/edd26/TopologyPreprocessing/actions/workflows/CI_julia.yml/badge.svg)


## Installation

This package registeration is being processed now. After being registered, to install it, run the following.

```julia
julia> using Pkg
julia> Pkg.add("TopologyPreprocessing")
```

## Optional dependencies

**Eirene** and **Plots**/**StatsPlots** are optional (weak) dependencies. The core matrix preprocessing functions work without them. Additional functionality is unlocked when these packages are loaded:

- `using Eirene` — enables `get_bettis`, `get_barcodes`, `get_vectorized_bettis` and related persistent homology functions
- `using Plots, StatsPlots` — enables all plotting functions (`plot_bettis`, `plot_barcodes`, `boxplot_lifetime`, etc.)

### Note on Eirene version

The registered version of Eirene (1.3.x) has outdated dependencies that may cause installation issues. It is recommended to use **Eirene 1.4 from the master branch** instead:

```julia
using Pkg
Pkg.add(url="https://github.com/henselman-petrusek/Eirene.jl")
```

Once Eirene is installed and loaded, all Eirene-dependent functions activate automatically:

```julia
using TopologyPreprocessing
using Eirene

C = eirene(my_matrix, maxdim=3, model="vr")
bettis = get_bettis(C, 3)
```
