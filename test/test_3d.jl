module Test3D

using Test
using Trixi2Vtk

include("test_trixi2vtk.jl")

# Start with a clean environment: remove Trixi output directory if it exists
outdir = "out"
isdir(outdir) && rm(outdir, recursive=true)


@testset "3D" begin
  @testset "TreeMesh" begin
    isdir(outdir) && rm(outdir, recursive=true)
    run_trixi(joinpath(examples_dir(), "tree_3d_dgsem", "elixir_euler_blob_amr.jl"), maxiters=4)

    @timed_testset "mesh data" begin
      @test_nowarn trixi2vtk(joinpath(outdir, "mesh_000004.h5"), output_directory=outdir)
      out_file = joinpath(outdir,"mesh_000004_celldata.vtu")
      remote_filename = joinpath("3d", "treemesh", "dgsem_blob_amr_mesh_04.vtu")
      ref_file = get_test_reference_file("dgsem_blob_amr_mesh_04.vtu", remote_filename)
      compare_cell_info(out_file, ref_file)
    end

    @timed_testset "solution celldata" begin
      @test_nowarn trixi2vtk(joinpath(outdir, "solution_000004.h5"), output_directory=outdir)
      out_file = joinpath(outdir,"solution_000004_celldata.vtu")
      remote_filename = joinpath("3d", "treemesh", "dgsem_blob_amr_celldata_04.vtu")
      ref_file = get_test_reference_file("dgsem_blob_amr_celldata_04.vtu", remote_filename)
      compare_cell_info(out_file, ref_file)
    end

    @timed_testset "reinterpolate with nonuniform data" begin
      rm(joinpath(outdir, "solution_000004.vtu")) # delete existing output file
      # Create and test output with reinterpolation (default options: `reinterpolate=true, data_is_uniform=false`)
      @test_nowarn trixi2vtk(joinpath(outdir, "solution_000004.h5"), output_directory=outdir)
      out_file = joinpath(outdir,"solution_000004.vtu")
      remote_filename = joinpath("3d", "treemesh", "dgsem_blob_amr_reinterp_04.vtu")
      ref_file = get_test_reference_file("dgsem_blob_amr_reinterp_04.vtu", remote_filename)
      compare_cell_info(out_file, ref_file)
    end

    @timed_testset "do not reinterpolate with nonuniform data" begin
      rm(joinpath(outdir, "solution_000004.vtu")) # delete existing output file
      # Create and test output without reinterpolation on LGL nodes
      @test_nowarn trixi2vtk(joinpath(outdir, "solution_000004.h5"), output_directory=outdir, reinterpolate=false)
      out_file = joinpath(outdir,"solution_000004.vtu")
      remote_filename = joinpath("3d", "treemesh", "dgsem_blob_amr_no_reinterp_04.vtu")
      ref_file = get_test_reference_file("dgsem_blob_amr_no_reinterp_04.vtu", remote_filename)
      compare_point_info(out_file, ref_file)
    end

    @timed_testset "do not reinterpolate with uniform data" begin
      rm(joinpath(outdir, "solution_000004.vtu")) # delete existing output file
      # Create and test output without reinterpolation on uniform nodes
      # OBS! This is a dummy test just to exercise code. The resulting plot will look weird.
      @test_nowarn trixi2vtk(joinpath(outdir, "solution_000004.h5"), output_directory=outdir, reinterpolate=false, data_is_uniform=true)
      out_file = joinpath(outdir,"solution_000004.vtu")
      remote_filename = joinpath("3d", "treemesh", "dgsem_blob_amr_no_reinterp_uniform_04.vtu")
      ref_file = get_test_reference_file("dgsem_blob_amr_no_reinterp_uniform_04.vtu", remote_filename)
      compare_point_info(out_file, ref_file)
    end
  end

  @testset "StructuredMesh" begin
    isdir(outdir) && rm(outdir, recursive=true)
    run_trixi(joinpath(examples_dir(), "structured_3d_dgsem", "elixir_advection_nonperiodic_curved.jl"), maxiters=1)

    @timed_testset "mesh data" begin
      @test_nowarn trixi2vtk(joinpath(outdir, "mesh.h5"), output_directory=outdir)
      out_file = joinpath(outdir,"mesh_celldata.vtu")
      remote_filename = joinpath("3d", "structuredmesh", "dgsem_adv_mesh_01.vtu")
      ref_file = get_test_reference_file("dgsem_adv_mesh_01.vtu", remote_filename)
      compare_cell_info(out_file, ref_file)
    end

    @timed_testset "solution celldata" begin
      @test_nowarn trixi2vtk(joinpath(outdir, "solution_000001.h5"), output_directory=outdir)
      out_file = joinpath(outdir,"solution_000001_celldata.vtu")
      remote_filename = joinpath("3d", "structuredmesh", "dgsem_adv_celldata_01.vtu")
      ref_file = get_test_reference_file("dgsem_adv_celldata_01.vtu", remote_filename)
      compare_cell_info(out_file, ref_file)
    end

    @timed_testset "reinterpolate with nonuniform data" begin
      rm(joinpath(outdir, "solution_000001.vtu")) # delete existing output file
      # Create and test output with reinterpolation (default options: `reinterpolate=true, data_is_uniform=false`)
      @test_nowarn trixi2vtk(joinpath(outdir, "solution_000001.h5"), output_directory=outdir)
      out_file = joinpath(outdir,"solution_000001.vtu")
      remote_filename = joinpath("3d", "structuredmesh", "dgsem_adv_reinterp_01.vtu")
      ref_file = get_test_reference_file("dgsem_adv_reinterp_01.vtu", remote_filename)
      compare_point_info(out_file, ref_file)
    end

    @timed_testset "do not reinterpolate with nonuniform data" begin
      rm(joinpath(outdir, "solution_000001.vtu")) # delete existing output file
      # Create and test output without reinterpolation on LGL nodes
      @test_nowarn trixi2vtk(joinpath(outdir, "solution_000001.h5"), output_directory=outdir, reinterpolate=false)
      out_file = joinpath(outdir,"solution_000001.vtu")
      remote_filename = joinpath("3d", "structuredmesh", "dgsem_adv_no_reinterp_01.vtu")
      ref_file = get_test_reference_file("dgsem_adv_no_reinterp_01.vtu", remote_filename)
      compare_point_info(out_file, ref_file)
    end

    @timed_testset "do not reinterpolate with uniform data" begin
      rm(joinpath(outdir, "solution_000001.vtu")) # delete existing output file
      # Create and test output without reinterpolation on uniform nodes
      # OBS! This is a dummy test just to exercise code. The resulting plot will look weird.
      @test_nowarn trixi2vtk(joinpath(outdir, "solution_000001.h5"), output_directory=outdir, reinterpolate=false, data_is_uniform=true)
      out_file = joinpath(outdir,"solution_000001.vtu")
      remote_filename = joinpath("3d", "structuredmesh", "dgsem_adv_no_reinterp_uniform_01.vtu")
      ref_file = get_test_reference_file("dgsem_adv_no_reinterp_uniform_01.vtu", remote_filename)
      compare_point_info(out_file, ref_file)
    end
  end

  @testset "P4estMesh" begin
    isdir(outdir) && rm(outdir, recursive=true)
    run_trixi(joinpath(examples_dir(), "p4est_3d_dgsem", "elixir_advection_cubed_sphere.jl"), maxiters=2)

    @timed_testset "mesh data" begin
      @test_nowarn trixi2vtk(joinpath(outdir, "mesh.h5"), output_directory=outdir)
      out_file = joinpath(outdir,"mesh_celldata.vtu")
      remote_filename = joinpath("3d", "p4estmesh", "dgsem_adv_sphere_mesh_02.vtu")
      ref_file = get_test_reference_file("dgsem_adv_sphere_mesh_02.vtu", remote_filename)
      compare_cell_info(out_file, ref_file)
    end

    @timed_testset "solution celldata" begin
      @test_nowarn trixi2vtk(joinpath(outdir, "solution_000002.h5"), output_directory=outdir)
      out_file = joinpath(outdir,"solution_000002_celldata.vtu")
      remote_filename = joinpath("3d", "p4estmesh", "dgsem_adv_sphere_celldata_02.vtu")
      ref_file = get_test_reference_file("dgsem_adv_sphere_celldata_02.vtu", remote_filename)
      compare_cell_info(out_file, ref_file)
    end

    @timed_testset "reinterpolate with nonuniform data" begin
      rm(joinpath(outdir, "solution_000002.vtu")) # delete existing output file
      # Create and test output with reinterpolation (default options: `reinterpolate=true, data_is_uniform=false`)
      @test_nowarn trixi2vtk(joinpath(outdir, "solution_000002.h5"), output_directory=outdir)
      out_file = joinpath(outdir,"solution_000002.vtu")
      remote_filename = joinpath("3d", "p4estmesh", "dgsem_adv_sphere_reinterp_02.vtu")
      ref_file = get_test_reference_file("dgsem_adv_sphere_reinterp_02.vtu", remote_filename)
      compare_point_info(out_file, ref_file)
    end

    @timed_testset "do not reinterpolate with nonuniform data" begin
      rm(joinpath(outdir, "solution_000002.vtu")) # delete existing output file
      # Create and test output without reinterpolation on LGL nodes
      @test_nowarn trixi2vtk(joinpath(outdir, "solution_000002.h5"), output_directory=outdir, reinterpolate=false)
      out_file = joinpath(outdir,"solution_000002.vtu")
      remote_filename = joinpath("3d", "p4estmesh", "dgsem_adv_sphere_no_reinterp_02.vtu")
      ref_file = get_test_reference_file("dgsem_adv_sphere_no_reinterp_02.vtu", remote_filename)
      compare_point_info(out_file, ref_file)
    end

    @timed_testset "do not reinterpolate with uniform data" begin
      rm(joinpath(outdir, "solution_000002.vtu")) # delete existing output file
      # Create and test output without reinterpolation on uniform nodes
      # OBS! This is a dummy test just to exercise code. The resulting plot will look weird.
      @test_nowarn trixi2vtk(joinpath(outdir, "solution_000002.h5"), output_directory=outdir, reinterpolate=false, data_is_uniform=true)
      out_file = joinpath(outdir,"solution_000002.vtu")
      remote_filename = joinpath("3d", "p4estmesh", "dgsem_adv_sphere_no_reinterp_uniform_02.vtu")
      ref_file = get_test_reference_file("dgsem_adv_sphere_no_reinterp_uniform_02.vtu", remote_filename)
      compare_point_info(out_file, ref_file)
    end
  end
end

# Clean up afterwards: delete Trixi output directory and reference file directory
@test_nowarn rm(outdir, recursive=true)
@test_nowarn rm(TEST_REFERENCE_DIR, recursive=true)

end
