module Trixi2Vtk

# Include other packages
using EllipsisNotation #..
using Glob: glob
using HDF5: h5open, attributes, haskey
using ProgressMeter: @showprogress, Progress, next!
using StaticArrays: SVector
using TimerOutputs
using Trixi: Trixi, transfinite_mapping, coordinates2mapping, polynomial_interpolation_matrix,
             gauss_lobatto_nodes_weights, TreeMesh, StructuredMesh, UnstructuredMesh2D, P4estMesh, T8codeMesh
using WriteVTK: vtk_grid, MeshCell, VTKCellTypes, vtk_save, paraview_collection

# Include all top-level submodule files
include("interpolate.jl")
include("io.jl")
include("pointlocators.jl")
include("vtktools.jl")

# Include top-level conversion method
include("convert.jl")

# export types/functions that define the public API of Trixi2Vtk
export trixi2vtk

end # module Trixi2Vtk
