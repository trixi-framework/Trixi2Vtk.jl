module TestManual

using Test
using Documenter
using Trixi2Vtk

@testset "Manual" begin
  DocMeta.setdocmeta!(Trixi2Vtk, :DocTestSetup, :(using Trixi2Vtk); recursive=true)
  doctest(Trixi2Vtk, manual=false)
end

end # module
