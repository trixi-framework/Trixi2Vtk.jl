module Test2D

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


@testset "2D" begin
  @testset "TreeMesh" begin
    isdir(outdir) && rm(outdir, recursive=true)
    run_trixi(joinpath("tree_2d_dgsem", "elixir_advection_extended.jl"), maxiters=1)

    @timed_testset "uniform mesh" begin
      if Sys.iswindows()
        # This test fails on Windows due to globbing not working
        test_trixi2vtk("solution_000000.h5", outdir,
            hashes=[("solution_000000.vtu", "1ec2c93c0c9c4f4992dea54afaf2a348ece0160e"),
                    ("solution_000000_celldata.vtu", "e396c3ba63276347966d4264cf0f52d592221830")])
        outfiles = ("solution_000000.vtu", "solution_000000_celldata.vtu")

      else
        test_trixi2vtk("solution_00000*.h5", outdir,
            hashes=[("solution_000000.vtu", "1ec2c93c0c9c4f4992dea54afaf2a348ece0160e"),
                    ("solution_000000_celldata.vtu", "e396c3ba63276347966d4264cf0f52d592221830"),
                    ("solution_00000.pvd", "b9e6742dc2b397b14d8d3964e90dcfadea5d98cb"),
                    ("solution_00000_celldata.pvd", "c12004714a980581450cd4bad16a2541c5ec8f26")])
        outfiles = ("solution_000000.vtu", "solution_000000_celldata.vtu",
                    "solution_00000.pvd", "solution_00000_celldata.pvd")
      end
      @test_nowarn trixi2vtk("mesh.h5")

      # Store output files as artifacts to facilitate debugging of failing tests
      testname = "2d-tree-mesh-uniform-mesh"
      for outfile in outfiles
        println("Copying '", abspath(joinpath(outdir, outfile)),
                "' to '", abspath(joinpath(artifacts_dir, testname * "-" * outfile)),
                "'...")
        cp(joinpath(outdir, outfile), joinpath(artifacts_dir, testname * "-" * outfile), force=true)
      end
    end

    @timed_testset "uniform mesh with vti output" begin
      if Sys.isapple()
        # This test fails on MacOS due to differing binary VTI files (even though the contents match)
        test_trixi2vtk("restart_000001.h5", outdir,
            format=:vti)
      else
        test_trixi2vtk("restart_000001.h5", outdir,
            hashes=[("restart_000001.vti", "9ade067c71f1f6492242a8aa215bd0d633caf9bc"),
                    ("restart_000001_celldata.vtu", "e396c3ba63276347966d4264cf0f52d592221830")],
            format=:vti)
      end

      # Store output files as artifacts to facilitate debugging of failing tests
      outfiles = ("restart_000001.vti", "restart_000001_celldata.vtu")
      testname = "2d-tree-mesh-uniform-mesh-with-vti-output"
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
    run_trixi(joinpath("structured_2d_dgsem", "elixir_advection_waving_flag.jl"), maxiters=1)

    @timed_testset "waving flag" begin
      if Sys.isapple()
        # This file has a different hash on macOS for some reason
        test_trixi2vtk("solution_000000.h5", outdir,
            hashes=[("solution_000000.vtu", "a93dbd393647627a861d890568e65598be0062f9")])
      else
        test_trixi2vtk("solution_000000.h5", outdir,
            hashes=[("solution_000000.vtu", "564701ed0a9a90230f3a67f8bddd0616c818319b")])
      end
      @test_nowarn trixi2vtk("mesh.h5")

      # Store output files as artifacts to facilitate debugging of failing tests
      outfiles = ("solution_000000.vtu",)
      testname = "2d-curved-mesh-waving-flag"
      for outfile in outfiles
        println("Copying '", abspath(joinpath(outdir, outfile)),
                "' to '", abspath(joinpath(artifacts_dir, testname * "-" * outfile)),
                "'...")
        cp(joinpath(outdir, outfile), joinpath(artifacts_dir, testname * "-" * outfile), force=true)
      end
    end

    @timed_testset "waving flag (supersampling)" begin
      if Sys.isapple()
        # This file has a different hash on macOS for some reason
        test_trixi2vtk("solution_000000.h5", outdir, nvisnodes=6,
            hashes=[("solution_000000.vtu", "c6d74ab831bf4b6de2ba8cf537b6653ad611cfe7")])
      else
        test_trixi2vtk("solution_000000.h5", outdir, nvisnodes=6,
            hashes=[("solution_000000.vtu", "ae8c059c110aaabe2ed7dcfa8516d336c15ba618")])
      end

      # Store output files as artifacts to facilitate debugging of failing tests
      outfiles = ("solution_000000.vtu",)
      testname = "2d-curved-mesh-waving-flag-supersampling"
      for outfile in outfiles
        println("Copying '", abspath(joinpath(outdir, outfile)),
                "' to '", abspath(joinpath(artifacts_dir, testname * "-" * outfile)),
                "'...")
        cp(joinpath(outdir, outfile), joinpath(artifacts_dir, testname * "-" * outfile), force=true)
      end
    end
  end

  @testset "UnstructuredMesh2D" begin
    isdir(outdir) && rm(outdir, recursive=true)
    run_trixi(joinpath("unstructured_2d_dgsem", "elixir_euler_basic.jl"), maxiters=1)

    @timed_testset "basic" begin
      if Sys.isapple()
        # This file has a different hash on macOS for some reason
        test_trixi2vtk("solution_000000.h5", outdir,
            hashes=[("solution_000000.vtu", "0daedeea99d03d53b925ce5691bd7924abe88861")])
      else
        test_trixi2vtk("solution_000000.h5", outdir,
            hashes=[("solution_000000.vtu", "acc03f295d2b7daee2cb0b4e29b90035c5e92bd7")])
      end
      @test_nowarn trixi2vtk("mesh.h5")

      # Store output files as artifacts to facilitate debugging of failing tests
      outfiles = ("solution_000000.vtu",)
      testname = "2d-unstructured-quad-basic"
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
    run_trixi(joinpath("p4est_2d_dgsem", "elixir_euler_source_terms_nonconforming_unstructured_flag.jl"), initial_refinement_level=0, maxiters=1)

    @timed_testset "nonperiodic" begin
      if Sys.isapple()
        # This file has a different hash on macOS for some reason
        test_trixi2vtk("solution_000000.h5", outdir,
          hashes=[("solution_000000.vtu", "476eaaa9c13c3a8dba826b614a68a9d3e7979e6b")])
      else
        test_trixi2vtk("solution_000000.h5", outdir,
          hashes=[("solution_000000.vtu", "a80aadb353ce6ec40baa1b94d278c480f17d0419")])
      end
      @test_nowarn trixi2vtk("mesh.h5")

      # Store output files as artifacts to facilitate debugging of failing tests
      outfiles = ("solution_000000.vtu",)
      testname = "2d-p4est-mesh-nonperiodic"
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
@test_nowarn rm(outdir, recursive=true)

end
