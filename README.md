# Trixi2Vtk.jl

<!-- [![Docs-stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://trixi-framework.github.io/Trixi2Vtk.jl/stable) -->
[![Docs-dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://trixi-framework.github.io/Trixi2Vtk.jl/dev)
[![Build Linux & macOS](https://travis-ci.com/trixi-framework/Trixi2Vtk.jl.svg?branch=master)](https://travis-ci.com/trixi-framework/Trixi2Vtk.jl)
[![Build Windows](https://ci.appveyor.com/api/projects/status/3vw2i3iy9n9641wq?svg=true)](https://ci.appveyor.com/project/ranocha/trixi2vtk-jl)
[![Codecov](https://codecov.io/gh/trixi-framework/Trixi2Vtk.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/trixi-framework/Trixi2Vtk.jl)
[![Coveralls](https://coveralls.io/repos/github/trixi-framework/Trixi2Vtk.jl/badge.svg?branch=master)](https://coveralls.io/github/trixi-framework/Trixi2Vtk.jl?branch=master)
[![License: MIT](https://img.shields.io/badge/License-MIT-success.svg)](https://opensource.org/licenses/MIT)
<!-- [![GitHub commits since tagged version](https://img.shields.io/github/commits-since/trixi-framework/Trixi2Vtk.jl/v0.1.0.svg?style=social&logo=github)](https://github.com/trixi-framework/Trixi2Vtk.jl) -->

With **Trixi2Vtk.jl** you can convert the HDF5-based output files created by
[Trixi.jl](https://github.com/trixi-framework/Trixi.jl) (solution or restart
files) to VTK files. The converted files can then be further processed with
[ParaView](https://www.paraview.org) to generate publication-quality
visualizations. Trixi2Vtk is part of the [Trixi
framework](https://github.com/trixi-framework).


## Installation
If you have not yet installed Julia, please follow the instructions for your
operating system found [here](https://julialang.org/downloads/platform/).
Trixi2Vtk works with Julia v1.5.

You can then install Trixi2Vtk, the postprocessing tools, and the respective dependencies by
performing the following steps:

  1. Clone the repository:
     ```bash
     git clone git@github.com:trixi-framework/Trixi2Vtk.jl.git
     ```
  2. Enter the cloned directory and run the following command to install all
     required dependencies:
     ```bash
     julia --project=. -e 'import Pkg; Pkg.instantiate()'
     ```


## Usage
Enter the root directory `Trixi2Vtk.jl/` and execute
```bash
julia --project=@.
```
This will start an interactive Julia session (REPL) using the project setup
of Trixi2Vtk.jl. If you have installed Trixi2Vtk.jl in your default project environment,
you can just start Julia as usual
```bash
julia
```
In the Julia REPL, you need to load the package Trixi2Vtk
```julia
julia> using Trixi2Vtk
```
To process an HDF5 file generate by Trixi.jl, execute
```julia
Trixi2Vtk.convert("out/solution_000000.h5")
```
This will convert `out/solution_000000.h5` to an unstructured VTK file with a
`.vtu` file extension.

Sometimes it can be helpful to use Trixi2Vtk non-interactively in batch mode, e.g.,
when starting a simulation from another script. This is possible by directly passing
the code that shall be executed to Julia
```bash
julia -e 'using Trixi2Vtk; Trixi2Vtk.convert("out/restart_*")'
```


## Authors
Trixi2Vtk is maintained by the
[Trixi authors](https://github.com/trixi-framework/Trixi.jl/blob/master/AUTHORS.md).
Its principal developers are
[Michael Schlottke-Lakemper](https://www.mi.uni-koeln.de/NumSim/schlottke-lakemper)
(University of Cologne, Germany) and [Hendrik Ranocha](https://ranocha.de)
(KAUST, Saudi Arabia).


## License and contributing
Trixi2Vtk is licensed under the MIT license (see [LICENSE.md](LICENSE.md)).
