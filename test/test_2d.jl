module Test2D

using Test
using Trixi2Vtk

include("test_trixi2vtk.jl")

# Start with a clean environment: remove Trixi output directory if it exists
outdir = "out"
isdir(outdir) && rm(outdir, recursive=true)

# Create empty artifacts directory where all files that should be preserved will be stored
artifacts_dir = joinpath(pathof(Trixi2Vtk) |> dirname |> dirname, "artifacts")
isdir(artifacts_dir) && rm(artifacts_dir, recursive=true)
mkdir(artifacts_dir)


@testset "2D" begin
  @testset "TreeMesh" begin
    isdir(outdir) && rm(outdir, recursive=true)
    run_trixi(joinpath("2d", "elixir_advection_extended.jl"), maxiters=1)

    @testset "uniform mesh" begin
      test_trixi2vtk("solution_000000.h5", outdir,
          hashes=[("solution_000000.vtu", "1ec2c93c0c9c4f4992dea54afaf2a348ece0160e"),
                  ("solution_000000_celldata.vtu", "9b20ba10df0d2d0fbd15916e5da0ed72ade9890b")])

      # Store output files as artifacts to facilitate debugging of failing tests
      outfiles = ("solution_000000.vtu", "solution_000000_celldata.vtu")
      testname = "tree-mesh-uniform-mesh"
      for outfile in outfiles
        println("Copying '", abspath(joinpath(outdir, outfile)),
                "' to '", abspath(joinpath(artifacts_dir, testname * "-" * outfile)),
                "'...")
        cp(joinpath(outdir, outfile), joinpath(artifacts_dir, testname * "-" * outfile))
      end
    end

    @testset "uniform mesh with vti output" begin
      if Sys.isapple()
        # This test fails on MacOS due to differing binary VTI files (even though the contents match)
        test_trixi2vtk("restart_000001.h5", outdir,
            format=:vti)
      else
        test_trixi2vtk("restart_000001.h5", outdir,
            hashes=[("restart_000001.vti", "664f25ab018a373774b5aad69ad3f2f5a3b21649"),
                    ("restart_000001_celldata.vtu", "9b20ba10df0d2d0fbd15916e5da0ed72ade9890b")],
            format=:vti)
      end

      # Store output files as artifacts to facilitate debugging of failing tests
      outfiles = ("restart_000001.vti", "restart_000001_celldata.vtu")
      testname = "tree-mesh-uniform-mesh-with-vti-output"
      for outfile in outfiles
        println("Copying '", abspath(joinpath(outdir, outfile)),
                "' to '", abspath(joinpath(artifacts_dir, testname * "-" * outfile)),
                "'...")
        cp(joinpath(outdir, outfile), joinpath(artifacts_dir, testname * "-" * outfile))
      end
    end
  end

  @testset "CurvedMesh" begin
    isdir(outdir) && rm(outdir, recursive=true)
    run_trixi(joinpath("2d", "elixir_advection_waving_flag.jl"), maxiters=1)

    @testset "waving flag" begin
      test_trixi2vtk("solution_000000.h5", outdir,
          hashes=[("solution_000000.vtu", "a93dbd393647627a861d890568e65598be0062f9")])

      # Store output files as artifacts to facilitate debugging of failing tests
      outfiles = ("solution_000000.vtu",)
      testname = "curved-mesh-waving-flag"
      for outfile in outfiles
        println("Copying '", abspath(joinpath(outdir, outfile)),
                "' to '", abspath(joinpath(artifacts_dir, testname * "-" * outfile)),
                "'...")
        cp(joinpath(outdir, outfile), joinpath(artifacts_dir, testname * "-" * outfile))
      end
    end

    @testset "waving flag (supersampling)" begin
      test_trixi2vtk("solution_000000.h5", outdir, nvisnodes=6,
          hashes=[("solution_000000.vtu", "c6d74ab831bf4b6de2ba8cf537b6653ad611cfe7")])

      # Store output files as artifacts to facilitate debugging of failing tests
      outfiles = ("solution_000000.vtu",)
      testname = "curved-mesh-waving-flag-supersampling"
      for outfile in outfiles
        println("Copying '", abspath(joinpath(outdir, outfile)),
                "' to '", abspath(joinpath(artifacts_dir, testname * "-" * outfile)),
                "'...")
        cp(joinpath(outdir, outfile), joinpath(artifacts_dir, testname * "-" * outfile))
      end
    end
  end

  @testset "UnstructuredQuadMesh" begin
    isdir(outdir) && rm(outdir, recursive=true)
    run_trixi(joinpath("2d", "elixir_euler_unstructured_quad_basic.jl"), maxiters=1)

    @testset "basic" begin
      test_trixi2vtk("solution_000000.h5", outdir,
          hashes=[("solution_000000.vtu", "0daedeea99d03d53b925ce5691bd7924abe88861")])

      # Store output files as artifacts to facilitate debugging of failing tests
      outfiles = ("solution_000000.vtu",)
      testname = "unstructured-quad-basic"
      for outfile in outfiles
        println("Copying '", abspath(joinpath(outdir, outfile)),
                "' to '", abspath(joinpath(artifacts_dir, testname * "-" * outfile)),
                "'...")
        cp(joinpath(outdir, outfile), joinpath(artifacts_dir, testname * "-" * outfile))
      end
    end
  end
end

# Clean up afterwards: delete Trixi output directory
@test_nowarn rm(outdir, recursive=true)

end
