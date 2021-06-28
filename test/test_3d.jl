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

    @testset "uniform mesh" begin
      test_trixi2vtk("solution_000000.h5", outdir,
          hashes=[("solution_000000.vtu", "6ab3aa525851187ee0839e1d670a254a66be4ad7"),
                  ("solution_000000_celldata.vtu", "99c782d732a4d1f6764013c3fba2cdaddf3927ab")])

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

    @testset "basic" begin
      test_trixi2vtk("solution_000000.h5", outdir,
          hashes=[("solution_000000.vtu", "57c58480d8f99b7e9d365cb3fa71790db42de750")])

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
    run_trixi(joinpath("p4est_3d_dgsem", "elixir_advection_basic.jl"), maxiters=1)

    @testset "unstructured curved" begin
      test_trixi2vtk("solution_000000.h5", outdir,
          hashes=[("solution_000000.vtu", "2e4fd4eaedfcebacfa0e062305f697f11641c86c")])

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

