module Test2D

using Test
using Trixi2Vtk

include("test_trixi2vtk.jl")

# Start with a clean environment: remove Trixi output directory if it exists
outdir = "out"
isdir(outdir) && rm(outdir, recursive=true)


@testset "2D" begin
  run_trixi("parameters.toml", n_steps_max=1)

  @testset "uniform mesh" begin
    test_trixi2vtk_convert("solution_000000.h5", outdir,
        hashes=[("solution_000000.vtu", "1ec2c93c0c9c4f4992dea54afaf2a348ece0160e"),
                ("solution_000000_celldata.vtu", "e396c3ba63276347966d4264cf0f52d592221830")])
  end

  @testset "uniform mesh with vti output" begin
    test_trixi2vtk_convert("restart_000001.h5", outdir,
        hashes=[("restart_000001.vti", "a6f5f7f1c5610ab58658335afba57b3427a6e731"),
                ("restart_000001_celldata.vtu", "e396c3ba63276347966d4264cf0f52d592221830")],
        format=:vti)
  end
end

# Clean up afterwards: delete Trixi output directory
@test_nowarn rm(outdir, recursive=true)

end
