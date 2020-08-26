module Test2D

using Test
using Trixi2Vtk

# Start with a clean environment: remove Trixi output directory if it exists
outdir = "out"
isdir(outdir) && rm(outdir, recursive=true)


@testset "2D" begin
  run_trixi("2d/parameters.toml", n_steps_max=1)

  @testset "uniform mesh" begin
    test_trixi2vtk("solution_000000.h5", outdir,
        hashes=[("solution_000000.vtu", "56ba70356fd6761432cd43d807c91961b3e7832e"),
                ("solution_000000_celldata.vtu", "318de9b537e9724c190c4499e4cf8dca50ffff14")])
  end
end

# Clean up afterwards: delete Trixi output directory
@test_nowarn rm(outdir, recursive=true)

end
