module Test3D

using Test
using Trixi2Vtk

include("test_trixi2vtk.jl")

# Start with a clean environment: remove Trixi output directory if it exists
outdir = "out"
isdir(outdir) && rm(outdir, recursive=true)

# Create artifacts directory where all files that should be preserved will be stored
artifacts_dir = joinpath(pathof(Trixi2Vtk) |> dirname |> dirname, "artifacts")
if !isdir(artifacts_dir)
  mkdir(artifacts_dir)
end


@testset "3D" begin
  @testset "TreeMesh" begin
    isdir(outdir) && rm(outdir, recursive=true)
    run_trixi(joinpath("tree_3d_dgsem", "elixir_advection_extended.jl"), maxiters=1)

    @timed_testset "uniform mesh" begin
      test_trixi2vtk("solution_000000.h5", outdir,
          hashes=[("solution_000000.vtu", "6ab3aa525851187ee0839e1d670a254a66be4ad7"),
                  ("solution_000000_celldata.vtu", "fdfee2d4200ecdad08067b37908412813016f4e7")])
      @test_nowarn trixi2vtk("mesh.h5")

      # Store output files as artifacts to facilitate debugging of failing tests
      outfiles = ("solution_000000.vtu", "solution_000000_celldata.vtu")
      testname = "3d-tree-mesh-uniform-mesh"
      for outfile in outfiles
        println("Copying '", abspath(joinpath(outdir, outfile)),
                "' to '", abspath(joinpath(artifacts_dir, testname * "-" * outfile)),
                "'...")
        cp(joinpath(outdir, outfile), joinpath(artifacts_dir, testname * "-" * outfile), force=true)
      end
    end
  end

  @testset "StructuredMesh" begin
    isdir(outdir) && rm(outdir, recursive=true)
    run_trixi(joinpath("structured_3d_dgsem", "elixir_advection_basic.jl"), maxiters=1)

    @timed_testset "basic" begin
      if Sys.isapple()
        # This file has a different hash on macOS for some reason
        test_trixi2vtk("solution_000000.h5", outdir,
            hashes=[("solution_000000.vtu", "58e07f981fd6c005ea17e47054bd509c2c66d771")])
      else
        test_trixi2vtk("solution_000000.h5", outdir,
            hashes=[("solution_000000.vtu", "58e07f981fd6c005ea17e47054bd509c2c66d771")])
      end
      @test_nowarn trixi2vtk("mesh.h5")

      # Store output files as artifacts to facilitate debugging of failing tests
      outfiles = ("solution_000000.vtu",)
      testname = "3d-curved-mesh-basic"
      for outfile in outfiles
        println("Copying '", abspath(joinpath(outdir, outfile)),
                "' to '", abspath(joinpath(artifacts_dir, testname * "-" * outfile)),
                "'...")
        cp(joinpath(outdir, outfile), joinpath(artifacts_dir, testname * "-" * outfile), force=true)
      end
    end
  end

  @testset "P4estMesh" begin
    isdir(outdir) && rm(outdir, recursive=true)
    run_trixi(joinpath("p4est_3d_dgsem", "elixir_advection_amr_unstructured_curved.jl"), maxiters=1)

    @timed_testset "unstructured curved" begin
      if Sys.isapple()
        # This file has a different hash on macOS for some reason
        test_trixi2vtk("solution_000000.h5", outdir,
          hashes=[("solution_000000.vtu", "fe0f45809ef6f3f0c1b7f5331198585f406923c9")])
      else
        test_trixi2vtk("solution_000000.h5", outdir,
          hashes=[("solution_000000.vtu", "0fa5a099378d153aa3a1bb7dcf3559ea5d6bf9c5")])
      end
      @test_nowarn trixi2vtk("mesh.h5")

      # Store output files as artifacts to facilitate debugging of failing tests
      outfiles = ("solution_000000.vtu",)
      testname = "3d-p4est-mesh-unstructured-curved"
      for outfile in outfiles
        println("Copying '", abspath(joinpath(outdir, outfile)),
                "' to '", abspath(joinpath(artifacts_dir, testname * "-" * outfile)),
                "'...")
        cp(joinpath(outdir, outfile), joinpath(artifacts_dir, testname * "-" * outfile), force=true)
      end
    end
  end
end

# Clean up afterwards: delete Trixi output directory
@test_skip rm(outdir, recursive=true)

end

