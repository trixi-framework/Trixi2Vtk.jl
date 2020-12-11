module TestManual

using Test
using Documenter
using Trixi2Vtk

DocMeta.setdocmeta!(Trixi2Vtk, :DocTestSetup, :(using Trixi2Vtk); recursive=true)
doctest(Trixi2Vtk, manual=false)

end # module
