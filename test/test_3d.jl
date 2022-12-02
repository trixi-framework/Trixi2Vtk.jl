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

# Point to the directory containing the reference VTK files
# TODO: Not sure how to get these files
refdir = joinpath(pathof(Trixi2Vtk) |> dirname |> dirname, "test", "reference_files")

@testset "3D" begin
  @testset "TreeMesh" begin
    isdir(outdir) && rm(outdir, recursive=true)
    run_trixi(joinpath(examples_dir(), "tree_3d_dgsem", "elixir_euler_blob_amr.jl"), maxiters=4)

    @timed_testset "mesh data" begin
      @test_nowarn trixi2vtk(joinpath(outdir, "mesh_000004.h5"), output_directory=outdir)
      out_file = joinpath(outdir,"mesh_000004_celldata.vtu")
      ref_file = joinpath(refdir, "3d", "treemesh", "dgsem_blob_amr_mesh_04.vtu")
      compare_cell_info(out_file, ref_file)
    end

    @timed_testset "solution celldata" begin
      @test_nowarn trixi2vtk(joinpath(outdir, "solution_000004.h5"), output_directory=outdir)
      out_file = joinpath(outdir,"solution_000004_celldata.vtu")
      ref_file = joinpath(refdir, "3d", "treemesh", "dgsem_blob_amr_celldata_04.vtu")
      compare_cell_info(out_file, ref_file)
    end

    @timed_testset "reinterpolate with nonuniform data" begin
      # Create and test output with reinterpolation (default options: `reinterpolate=true, data_is_uniform=false`)
      @test_nowarn trixi2vtk(joinpath(outdir, "solution_000004.h5"), output_directory=outdir)
      out_file = joinpath(outdir,"solution_000004.vtu")
      ref_file = joinpath(refdir, "3d", "treemesh", "dgsem_blob_amr_reinterp_04.vtu")
      compare_cell_info(out_file, ref_file)
    end

    @timed_testset "do not reinterpolate with nonuniform data" begin
      # Create and test output without reinterpolation on LGL nodes
      @test_nowarn trixi2vtk(joinpath(outdir, "solution_000004.h5"), output_directory=outdir, reinterpolate=false)
      out_file = joinpath(outdir,"solution_000004.vtu")
      ref_file = joinpath(refdir, "3d", "treemesh", "dgsem_blob_amr_no_reinterp_04.vtu")
      compare_point_info(out_file, ref_file)
    end

    @timed_testset "do not reinterpolate with uniform data" begin
      # Create and test output without reinterpolation on uniform nodes
      # OBS! This is a dummy test just to exercise code. The resulting plot will look weird.
      @test_nowarn trixi2vtk(joinpath(outdir, "solution_000004.h5"), output_directory=outdir, reinterpolate=false, data_is_uniform=true)
      out_file = joinpath(outdir,"solution_000004.vtu")
      ref_file = joinpath(refdir, "3d", "treemesh", "dgsem_blob_amr_no_reinterp_uniform_04.vtu")
      compare_point_info(out_file, ref_file)
    end
  end

  @testset "StructuredMesh" begin
    isdir(outdir) && rm(outdir, recursive=true)
    run_trixi(joinpath(examples_dir(), "structured_3d_dgsem", "elixir_advection_nonperiodic_curved.jl"), maxiters=1)

    @timed_testset "mesh data" begin
      @test_nowarn trixi2vtk(joinpath(outdir, "mesh.h5"), output_directory=outdir)
      out_file = joinpath(outdir,"mesh_celldata.vtu")
      ref_file = joinpath(refdir, "3d", "structuredmesh", "dgsem_adv_mesh_01.vtu")
      compare_cell_info(out_file, ref_file)
    end

    @timed_testset "solution celldata" begin
      @test_nowarn trixi2vtk(joinpath(outdir, "solution_000001.h5"), output_directory=outdir)
      out_file = joinpath(outdir,"solution_000001_celldata.vtu")
      ref_file = joinpath(refdir, "3d", "structuredmesh", "dgsem_adv_celldata_01.vtu")
      compare_cell_info(out_file, ref_file)
    end

    @timed_testset "reinterpolate with nonuniform data" begin
      # Create and test output with reinterpolation (default options: `reinterpolate=true, data_is_uniform=false`)
      @test_nowarn trixi2vtk(joinpath(outdir, "solution_000001.h5"), output_directory=outdir)
      out_file = joinpath(outdir,"solution_000001.vtu")
      ref_file = joinpath(refdir, "3d", "structuredmesh", "dgsem_adv_reinterp_01.vtu")
      compare_point_info(out_file, ref_file)
    end

    @timed_testset "do not reinterpolate with nonuniform data" begin
      # Create and test output without reinterpolation on LGL nodes
      @test_nowarn trixi2vtk(joinpath(outdir, "solution_000001.h5"), output_directory=outdir, reinterpolate=false)
      out_file = joinpath(outdir,"solution_000001.vtu")
      ref_file = joinpath(refdir, "3d", "structuredmesh", "dgsem_adv_no_reinterp_01.vtu")
      compare_point_info(out_file, ref_file)
    end

    @timed_testset "do not reinterpolate with uniform data" begin
      # Create and test output without reinterpolation on uniform nodes
      # OBS! This is a dummy test just to exercise code. The resulting plot will look weird.
      @test_nowarn trixi2vtk(joinpath(outdir, "solution_000001.h5"), output_directory=outdir, reinterpolate=false, data_is_uniform=true)
      out_file = joinpath(outdir,"solution_000001.vtu")
      ref_file = joinpath(refdir, "3d", "structuredmesh", "dgsem_adv_no_reinterp_uniform_01.vtu")
      compare_point_info(out_file, ref_file)
    end
  end

  @testset "P4estMesh" begin
    isdir(outdir) && rm(outdir, recursive=true)
    run_trixi(joinpath(examples_dir(), "p4est_3d_dgsem", "elixir_advection_cubed_sphere.jl"), maxiters=2)

    @timed_testset "mesh data" begin
      @test_nowarn trixi2vtk(joinpath(outdir, "mesh.h5"), output_directory=outdir)
      out_file = joinpath(outdir,"mesh_celldata.vtu")
      ref_file = joinpath(refdir, "3d", "p4estmesh", "dgsem_adv_sphere_mesh_02.vtu")
      compare_cell_info(out_file, ref_file)
    end

    @timed_testset "solution celldata" begin
      @test_nowarn trixi2vtk(joinpath(outdir, "solution_000002.h5"), output_directory=outdir)
      out_file = joinpath(outdir,"solution_000002_celldata.vtu")
      ref_file = joinpath(refdir, "3d", "p4estmesh", "dgsem_adv_sphere_celldata_02.vtu")
      compare_cell_info(out_file, ref_file)
    end

    @timed_testset "reinterpolate with nonuniform data" begin
      # Create and test output with reinterpolation (default options: `reinterpolate=true, data_is_uniform=false`)
      @test_nowarn trixi2vtk(joinpath(outdir, "solution_000002.h5"), output_directory=outdir)
      out_file = joinpath(outdir,"solution_000002.vtu")
      ref_file = joinpath(refdir, "3d", "p4estmesh", "dgsem_adv_sphere_reinterp_02.vtu")
      compare_point_info(out_file, ref_file)
    end

    @timed_testset "do not reinterpolate with nonuniform data" begin
      # Create and test output without reinterpolation on LGL nodes
      @test_nowarn trixi2vtk(joinpath(outdir, "solution_000002.h5"), output_directory=outdir, reinterpolate=false)
      out_file = joinpath(outdir,"solution_000002.vtu")
      ref_file = joinpath(refdir, "3d", "p4estmesh", "dgsem_adv_sphere_no_reinterp_02.vtu")
      compare_point_info(out_file, ref_file)
    end

    @timed_testset "do not reinterpolate with uniform data" begin
      # Create and test output without reinterpolation on uniform nodes
      # OBS! This is a dummy test just to exercise code. The resulting plot will look weird.
      @test_nowarn trixi2vtk(joinpath(outdir, "solution_000002.h5"), output_directory=outdir, reinterpolate=false, data_is_uniform=true)
      out_file = joinpath(outdir,"solution_000002.vtu")
      ref_file = joinpath(refdir, "3d", "p4estmesh", "dgsem_adv_sphere_no_reinterp_uniform_02.vtu")
      compare_point_info(out_file, ref_file)
    end
  end
end

# Clean up afterwards: delete Trixi output directory
@test_skip rm(outdir, recursive=true)

end

##
# All old test; very broken
##

# tset "3D" begin
#   @testset "TreeMesh" begin
#     isdir(outdir) && rm(outdir, recursive=true)
#     run_trixi(joinpath(examples_dir(), "tree_3d_dgsem", "elixir_advection_extended.jl"), maxiters=1)

#     @timed_testset "uniform mesh" begin
#       test_trixi2vtk("solution_000000.h5", outdir,
#           hashes=[("solution_000000.vtu", "6ab3aa525851187ee0839e1d670a254a66be4ad7"),
#                   ("solution_000000_celldata.vtu", "fdfee2d4200ecdad08067b37908412813016f4e7")])
#       @test_nowarn trixi2vtk(joinpath(outdir, "mesh.h5"))

#       # Store output files as artifacts to facilitate debugging of failing tests
#       outfiles = ("solution_000000.vtu", "solution_000000_celldata.vtu")
#       testname = "3d-tree-mesh-uniform-mesh"
#       for outfile in outfiles
#         println("Copying '", abspath(joinpath(outdir, outfile)),
#                 "' to '", abspath(joinpath(artifacts_dir, testname * "-" * outfile)),
#                 "'...")
#         cp(joinpath(outdir, outfile), joinpath(artifacts_dir, testname * "-" * outfile), force=true)
#       end
#     end
#   end

#   @testset "StructuredMesh" begin
#     isdir(outdir) && rm(outdir, recursive=true)
#     run_trixi(joinpath(examples_dir(), "structured_3d_dgsem", "elixir_advection_basic.jl"), maxiters=1)

#     @timed_testset "basic" begin
#       if Sys.isapple()
#         # This file has a different hash on macOS for some reason
#         test_trixi2vtk("solution_000000.h5", outdir,
#             hashes=[("solution_000000.vtu", "58e07f981fd6c005ea17e47054bd509c2c66d771")])
#       else
#         test_trixi2vtk("solution_000000.h5", outdir,
#             hashes=[("solution_000000.vtu", "58e07f981fd6c005ea17e47054bd509c2c66d771")])
#       end
#       @test_nowarn trixi2vtk(joinpath(outdir, "mesh.h5"))

#       # Store output files as artifacts to facilitate debugging of failing tests
#       outfiles = ("solution_000000.vtu",)
#       testname = "3d-curved-mesh-basic"
#       for outfile in outfiles
#         println("Copying '", abspath(joinpath(outdir, outfile)),
#                 "' to '", abspath(joinpath(artifacts_dir, testname * "-" * outfile)),
#                 "'...")
#         cp(joinpath(outdir, outfile), joinpath(artifacts_dir, testname * "-" * outfile), force=true)
#       end
#     end
#   end

#   @testset "P4estMesh" begin
#     isdir(outdir) && rm(outdir, recursive=true)
#     run_trixi(joinpath(examples_dir(), "p4est_3d_dgsem", "elixir_advection_amr_unstructured_curved.jl"), maxiters=1)

#     @timed_testset "unstructured curved" begin
#       if Sys.isapple()
#         # This file has a different hash on macOS for some reason
#         test_trixi2vtk("solution_000000.h5", outdir,
#           hashes=[("solution_000000.vtu", "fe0f45809ef6f3f0c1b7f5331198585f406923c9")])
#       else
#         test_trixi2vtk("solution_000000.h5", outdir,
#           hashes=[("solution_000000.vtu", "0fa5a099378d153aa3a1bb7dcf3559ea5d6bf9c5")])
#       end
#       @test_nowarn trixi2vtk(joinpath(outdir, "mesh.h5"))

#       # Store output files as artifacts to facilitate debugging of failing tests
#       outfiles = ("solution_000000.vtu",)
#       testname = "3d-p4est-mesh-unstructured-curved"
#       for outfile in outfiles
#         println("Copying '", abspath(joinpath(outdir, outfile)),
#                 "' to '", abspath(joinpath(artifacts_dir, testname * "-" * outfile)),
#                 "'...")
#         cp(joinpath(outdir, outfile), joinpath(artifacts_dir, testname * "-" * outfile), force=true)
#       end
#     end
#   end
# end
