# Trixi2Vtk

**Trixi2Vtk.jl** can convert the HDF5-based output files created by
[Trixi.jl](https://github.com/trixi-framework/Trixi.jl) (solution or restart
files) to VTK files. The converted files can then be further processed with
[ParaView](https://www.paraview.org) to generate publication-quality
visualizations. Trixi2Vtk is part of the [Trixi framework](https://github.com/trixi-framework).

For further documentation of Trixi.jl and tutorials explaining how to use Trixi
with Trixi2Vtk, please refer to the [documentation of Trixi](https://trixi-framework.github.io/Trixi.jl/stable/).


## Installation

If you have not yet installed Julia, please follow the instructions for your
operating system found [here](https://julialang.org/downloads/platform/).
Trixi2Vtk works with Julia v1.5 or higher.
Official binaries are available for Windows, macOS, Linux, and FreeBSD.

Trixi2Vtk is a registered Julia package. Thus, you can install it via
```julia
julia> import Pkg

julia> Pkg.add("Trixi2Vtk")
```


## [Authors](@id authors-index-md)

Trixi2Vtk is maintained by the
[Trixi authors](https://github.com/trixi-framework/Trixi.jl/blob/master/AUTHORS.md).
Its principal developers are
[Michael Schlottke-Lakemper](https://www.mi.uni-koeln.de/NumSim/schlottke-lakemper)
(University of Cologne, Germany) and
[Hendrik Ranocha](https://ranocha.de) (KAUST, Saudi Arabia).


## License and contributing

Trixi is licensed under the MIT license (see [License](@ref)).
