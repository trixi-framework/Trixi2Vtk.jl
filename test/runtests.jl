using Test

@testset "Trixi2Vtk" begin
  include("test_2d.jl")
  include("test_3d.jl")
  include("test_manual.jl")
end
