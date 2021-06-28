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

    @testset "uniform mesh" begin
      if Sys.iswindows()
        # This test fails on Windows due to globbing not working
        test_trixi2vtk("solution_000000.h5", outdir,
            hashes=[("solution_000000.vtu", "1ec2c93c0c9c4f4992dea54afaf2a348ece0160e"),
                    ("solution_000000_celldata.vtu", "9b20ba10df0d2d0fbd15916e5da0ed72ade9890b")])
        outfiles = ("solution_000000.vtu", "solution_000000_celldata.vtu")

      else
        test_trixi2vtk("solution_00000*.h5", outdir,
            hashes=[("solution_000000.vtu", "1ec2c93c0c9c4f4992dea54afaf2a348ece0160e"),
                    ("solution_000000_celldata.vtu", "9b20ba10df0d2d0fbd15916e5da0ed72ade9890b"),
                    ("solution_00000.pvd", "7ba2f8f1927e90ebd4209aab890c58a20acf63f4"),
                    ("solution_00000_celldata.pvd", "448a7130a608ed9f7e4630033b9e1338b1403f7b")])
        outfiles = ("solution_000000.vtu", "solution_000000_celldata.vtu",
                    "solution_00000.pvd", "solution_00000_celldata.pvd")
      end

      # Store output files as artifacts to facilitate debugging of failing tests
      testname = "2d-tree-mesh-uniform-mesh"
      for outfile in outfiles
        println("Copying '", abspath(joinpath(outdir, outfile)),
                "' to '", abspath(joinpath(artifacts_dir, testname * "-" * outfile)),
                "'...")
        cp(joinpath(outdir, outfile), joinpath(artifacts_dir, testname * "-" * outfile), force=true)
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

    @testset "waving flag" begin
      test_trixi2vtk("solution_000000.h5", outdir,
          hashes=[("solution_000000.vtu", "a93dbd393647627a861d890568e65598be0062f9")])

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

    @testset "waving flag (supersampling)" begin
      test_trixi2vtk("solution_000000.h5", outdir, nvisnodes=6,
          hashes=[("solution_000000.vtu", "c6d74ab831bf4b6de2ba8cf537b6653ad611cfe7")])

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

    @testset "basic" begin
      test_trixi2vtk("solution_000000.h5", outdir,
          hashes=[("solution_000000.vtu", "0daedeea99d03d53b925ce5691bd7924abe88861")])

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
    run_trixi(joinpath("p4est_2d_dgsem", "elixir_euler_source_terms_nonperiodic.jl"), initial_refinement_level=1, maxiters=1)

    @testset "nonperiodic" begin
      if Sys.isapple()
        # This file has a different hash on macOS for some reason
        test_trixi2vtk("solution_000000.h5", outdir,
          hashes=[("solution_000000.vtu", "9050e4d141bf4099db907ef2c55768cff98a7291")])
      else
        test_trixi2vtk("solution_000000.h5", outdir,
          hashes=[("solution_000000.vtu", "0aca9d5a7b112d0469fa76b215e6c2d48b038af4")])
      end

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
