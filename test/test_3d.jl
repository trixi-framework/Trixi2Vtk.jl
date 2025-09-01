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
    run_trixi(joinpath(examples_dir(), "tree_3d_dgsem", "elixir_euler_blob_amr.jl"), maxiters=4)

    @timed_testset "mesh data" begin
      # create the output file to be tested
      @test_nowarn trixi2vtk(joinpath(outdir, "mesh_" * LEADING_ZEROS * "000004.h5"), output_directory=outdir)
      outfilename = "mesh_" * LEADING_ZEROS * "000004_celldata.vtu"
      out_file = joinpath(outdir, outfilename)

      # save output file to `artifacts` to facilitate debugging of failing tests
      testname = "3d-treemesh"
      cp(out_file, joinpath(artifacts_dir, testname * "-" * outfilename), force=true)

      # remote file path is actually a URL so it always has the same path structure
      remote_filename = "3d/treemesh/dgsem_blob_amr_mesh_04.vtu"
      ref_file = get_test_reference_file("dgsem_blob_amr_mesh_04.vtu", remote_filename)
      compare_cell_data(out_file, ref_file)
    end

    @timed_testset "solution celldata" begin
      # create the output file to be tested
      @test_nowarn trixi2vtk(joinpath(outdir, "solution_" * LEADING_ZEROS * "000004.h5"), output_directory=outdir)
      outfilename = "solution_" * LEADING_ZEROS * "000004_celldata.vtu"
      out_file = joinpath(outdir, outfilename)

      # save output file to `artifacts` to facilitate debugging of failing tests
      testname = "3d-treemesh"
      cp(out_file, joinpath(artifacts_dir, testname * "-" * outfilename), force=true)

      # remote file path is actually a URL so it always has the same path structure
      remote_filename = "3d/treemesh/dgsem_blob_amr_celldata_04.vtu"
      ref_file = get_test_reference_file("dgsem_blob_amr_celldata_04.vtu", remote_filename)
      compare_cell_data(out_file, ref_file)
    end

    @timed_testset "reinterpolate with nonuniform data" begin
      # Create and test output with reinterpolation (default options: `reinterpolate=true, data_is_uniform=false`)
      @test_nowarn trixi2vtk(joinpath(outdir, "solution_" * LEADING_ZEROS * "000004.h5"), output_directory=outdir)
      outfilename = "solution_" * LEADING_ZEROS * "000004.vtu"
      out_file = joinpath(outdir, outfilename)

      # save output file to `artifacts` to facilitate debugging of failing tests
      testname = "3d-treemesh-reinterp"
      cp(out_file, joinpath(artifacts_dir, testname * "-" * outfilename), force=true)

      # remote file path is actually a URL so it always has the same path structure
      remote_filename = "3d/treemesh/dgsem_blob_amr_reinterp_04.vtu"
      ref_file = get_test_reference_file("dgsem_blob_amr_reinterp_04.vtu", remote_filename)
      compare_cell_data(out_file, ref_file)
    end

    @timed_testset "do not reinterpolate with nonuniform data" begin
      # Create and test output without reinterpolation on LGL nodes
      @test_nowarn trixi2vtk(joinpath(outdir, "solution_" * LEADING_ZEROS * "000004.h5"), output_directory=outdir, reinterpolate=false)
      outfilename = "solution_" * LEADING_ZEROS * "000004.vtu"
      out_file = joinpath(outdir, outfilename)

      # save output file to `artifacts` to facilitate debugging of failing tests
      testname = "3d-treemesh-no-reinterp"
      cp(out_file, joinpath(artifacts_dir, testname * "-" * outfilename), force=true)

      # remote file path is actually a URL so it always has the same path structure
      remote_filename = "3d/treemesh/dgsem_blob_amr_no_reinterp_04.vtu"
      ref_file = get_test_reference_file("dgsem_blob_amr_no_reinterp_04.vtu", remote_filename)
      compare_point_data(out_file, ref_file)
    end

    @timed_testset "do not reinterpolate with uniform data" begin
      # Create and test output without reinterpolation on uniform nodes
      # OBS! This is a dummy test just to exercise code. The resulting plot will look weird.
      @test_nowarn trixi2vtk(joinpath(outdir, "solution_" * LEADING_ZEROS * "000004.h5"), output_directory=outdir, reinterpolate=false, data_is_uniform=true)
      outfilename = "solution_" * LEADING_ZEROS * "000004.vtu"
      out_file = joinpath(outdir, outfilename)

      # save output file to `artifacts` to facilitate debugging of failing tests
      testname = "3d-treemesh-no-reinterp-uniform"
      cp(out_file, joinpath(artifacts_dir, testname * "-" * outfilename), force=true)

      # remote file path is actually a URL so it always has the same path structure
      remote_filename = "3d/treemesh/dgsem_blob_amr_no_reinterp_uniform_04.vtu"
      ref_file = get_test_reference_file("dgsem_blob_amr_no_reinterp_uniform_04.vtu", remote_filename)
      compare_point_data(out_file, ref_file)
    end
  end

  if !Sys.iswindows() && get(ENV, "CI", nothing) == "true"
    # OBS! Only `TreeMesh` results are tested on Windows runners due to memory limits.
    #      All remaining mesh types are tested on Ubuntu and Mac
    @testset "StructuredMesh" begin
      isdir(outdir) && rm(outdir, recursive=true)
      run_trixi(joinpath(examples_dir(), "structured_3d_dgsem", "elixir_advection_nonperiodic_curved.jl"), maxiters=1)

      @timed_testset "mesh data" begin
        # create the output file to be tested
        @test_nowarn trixi2vtk(joinpath(outdir, "mesh.h5"), output_directory=outdir)
        outfilename = "mesh_celldata.vtu"
        out_file = joinpath(outdir, outfilename)

        # save output file to `artifacts` to facilitate debugging of failing tests
        testname = "3d-structuredmesh"
        cp(out_file, joinpath(artifacts_dir, testname * "-" * outfilename), force=true)

        # remote file path is actually a URL so it always has the same path structure
        remote_filename = "3d/structuredmesh/dgsem_adv_mesh_01.vtu"
        ref_file = get_test_reference_file("dgsem_adv_mesh_01.vtu", remote_filename)
        compare_cell_data(out_file, ref_file)
      end

      @timed_testset "solution celldata" begin
        # create the output file to be tested
        @test_nowarn trixi2vtk(joinpath(outdir, "solution_" * LEADING_ZEROS * "000001.h5"), output_directory=outdir)
        outfilename = "solution_" * LEADING_ZEROS * "000001_celldata.vtu"
        out_file = joinpath(outdir, outfilename)

        # save output file to `artifacts` to facilitate debugging of failing tests
        testname = "3d-structuredmesh"
        cp(out_file, joinpath(artifacts_dir, testname * "-" * outfilename), force=true)

        # remote file path is actually a URL so it always has the same path structure
        remote_filename = "3d/structuredmesh/dgsem_adv_celldata_01.vtu"
        ref_file = get_test_reference_file("dgsem_adv_celldata_01.vtu", remote_filename)
        compare_cell_data(out_file, ref_file)
      end

      @timed_testset "reinterpolate with nonuniform data" begin
        # Create and test output with reinterpolation (default options: `reinterpolate=true, data_is_uniform=false`)
        @test_nowarn trixi2vtk(joinpath(outdir, "solution_" * LEADING_ZEROS * "000001.h5"), output_directory=outdir)
        outfilename = "solution_" * LEADING_ZEROS * "000001.vtu"
        out_file = joinpath(outdir, outfilename)

        # save output file to `artifacts` to facilitate debugging of failing tests
        testname = "3d-structuredmesh-reinterp"
        cp(out_file, joinpath(artifacts_dir, testname * "-" * outfilename), force=true)

        # remote file path is actually a URL so it always has the same path structure
        remote_filename = "3d/structuredmesh/dgsem_adv_reinterp_01.vtu"
        ref_file = get_test_reference_file("dgsem_adv_reinterp_01.vtu", remote_filename)
        compare_point_data(out_file, ref_file)
      end

      @timed_testset "do not reinterpolate with nonuniform data" begin
        # Create and test output without reinterpolation on LGL nodes
        @test_nowarn trixi2vtk(joinpath(outdir, "solution_" * LEADING_ZEROS * "000001.h5"), output_directory=outdir, reinterpolate=false)
        outfilename = "solution_" * LEADING_ZEROS * "000001.vtu"
        out_file = joinpath(outdir, outfilename)

        # save output file to `artifacts` to facilitate debugging of failing tests
        testname = "3d-structuredmesh-no-reinterp"
        cp(out_file, joinpath(artifacts_dir, testname * "-" * outfilename), force=true)

        # remote file path is actually a URL so it always has the same path structure
        remote_filename = "3d/structuredmesh/dgsem_adv_no_reinterp_01.vtu"
        ref_file = get_test_reference_file("dgsem_adv_no_reinterp_01.vtu", remote_filename)
        compare_point_data(out_file, ref_file)
      end

      @timed_testset "do not reinterpolate with uniform data" begin
        # Create and test output without reinterpolation on uniform nodes
        # OBS! This is a dummy test just to exercise code. The resulting plot will look weird.
        @test_nowarn trixi2vtk(joinpath(outdir, "solution_" * LEADING_ZEROS * "000001.h5"), output_directory=outdir, reinterpolate=false, data_is_uniform=true)
        outfilename = "solution_" * LEADING_ZEROS * "000001.vtu"
        out_file = joinpath(outdir, outfilename)

        # save output file to `artifacts` to facilitate debugging of failing tests
        testname = "3d-structuredmesh-no-reinterp-uniform"
        cp(out_file, joinpath(artifacts_dir, testname * "-" * outfilename), force=true)

        # remote file path is actually a URL so it always has the same path structure
        remote_filename = "3d/structuredmesh/dgsem_adv_no_reinterp_uniform_01.vtu"
        ref_file = get_test_reference_file("dgsem_adv_no_reinterp_uniform_01.vtu", remote_filename)
        compare_point_data(out_file, ref_file)
      end
    end

    @testset "P4estMesh" begin
      isdir(outdir) && rm(outdir, recursive=true)
      run_trixi(joinpath(examples_dir(), "p4est_3d_dgsem", "elixir_advection_cubed_sphere.jl"), maxiters=2)

      @timed_testset "mesh data" begin
        # create the output file to be tested
        @test_nowarn trixi2vtk(joinpath(outdir, "mesh.h5"), output_directory=outdir)
        outfilename = "mesh_celldata.vtu"
        out_file = joinpath(outdir, outfilename)

        # save output file to `artifacts` to facilitate debugging of failing tests
        testname = "3d-p4estmesh"
        cp(out_file, joinpath(artifacts_dir, testname * "-" * outfilename), force=true)

        # remote file path is actually a URL so it always has the same path structure
        remote_filename = "3d/p4estmesh/dgsem_adv_sphere_mesh_02.vtu"
        ref_file = get_test_reference_file("dgsem_adv_sphere_mesh_02.vtu", remote_filename)
        compare_cell_data(out_file, ref_file)
      end

      @timed_testset "solution celldata" begin
        # create the output file to be tested
        @test_nowarn trixi2vtk(joinpath(outdir, "solution_" * LEADING_ZEROS * "000002.h5"), output_directory=outdir)
        outfilename = "solution_" * LEADING_ZEROS * "000002_celldata.vtu"
        out_file = joinpath(outdir, outfilename)

        # save output file to `artifacts` to facilitate debugging of failing tests
        testname = "3d-p4estmesh"
        cp(out_file, joinpath(artifacts_dir, testname * "-" * outfilename), force=true)

        # remote file path is actually a URL so it always has the same path structure
        remote_filename = "3d/p4estmesh/dgsem_adv_sphere_celldata_02.vtu"
        ref_file = get_test_reference_file("dgsem_adv_sphere_celldata_02.vtu", remote_filename)
        compare_cell_data(out_file, ref_file)
      end

      @timed_testset "reinterpolate with nonuniform data" begin
        # Create and test output with reinterpolation (default options: `reinterpolate=true, data_is_uniform=false`)
        @test_nowarn trixi2vtk(joinpath(outdir, "solution_" * LEADING_ZEROS * "000002.h5"), output_directory=outdir)
        outfilename = "solution_" * LEADING_ZEROS * "000002.vtu"
        out_file = joinpath(outdir, outfilename)

        # save output file to `artifacts` to facilitate debugging of failing tests
        testname = "3d-p4estmesh-reinterp"
        cp(out_file, joinpath(artifacts_dir, testname * "-" * outfilename), force=true)

        # remote file path is actually a URL so it always has the same path structure
        remote_filename = "3d/p4estmesh/dgsem_adv_sphere_reinterp_02.vtu"
        ref_file = get_test_reference_file("dgsem_adv_sphere_reinterp_02.vtu", remote_filename)
        compare_point_data(out_file, ref_file)
      end

      @timed_testset "do not reinterpolate with nonuniform data" begin
        # Create and test output without reinterpolation on LGL nodes
        @test_nowarn trixi2vtk(joinpath(outdir, "solution_" * LEADING_ZEROS * "000002.h5"), output_directory=outdir, reinterpolate=false)
        outfilename = "solution_" * LEADING_ZEROS * "000002.vtu"
        out_file = joinpath(outdir, outfilename)

        # save output file to `artifacts` to facilitate debugging of failing tests
        testname = "3d-p4estmesh-no-reinterp"
        cp(out_file, joinpath(artifacts_dir, testname * "-" * outfilename), force=true)

        # remote file path is actually a URL so it always has the same path structure
        remote_filename = "3d/p4estmesh/dgsem_adv_sphere_no_reinterp_02.vtu"
        ref_file = get_test_reference_file("dgsem_adv_sphere_no_reinterp_02.vtu", remote_filename)
        compare_point_data(out_file, ref_file)
      end

      @timed_testset "do not reinterpolate with uniform data" begin
        # Create and test output without reinterpolation on uniform nodes
        # OBS! This is a dummy test just to exercise code. The resulting plot will look weird.
        @test_nowarn trixi2vtk(joinpath(outdir, "solution_" * LEADING_ZEROS * "000002.h5"), output_directory=outdir, reinterpolate=false, data_is_uniform=true)
        outfilename = "solution_" * LEADING_ZEROS * "000002.vtu"
        out_file = joinpath(outdir, outfilename)

        # save output file to `artifacts` to facilitate debugging of failing tests
        testname = "3d-p4estmesh-no-reinterp-uniform"
        cp(out_file, joinpath(artifacts_dir, testname * "-" * outfilename), force=true)

        # remote file path is actually a URL so it always has the same path structure
        remote_filename = "3d/p4estmesh/dgsem_adv_sphere_no_reinterp_uniform_02.vtu"
        ref_file = get_test_reference_file("dgsem_adv_sphere_no_reinterp_uniform_02.vtu", remote_filename)
        compare_point_data(out_file, ref_file)
      end
    end

    @testset "T8codeMesh" begin
      isdir(outdir) && rm(outdir, recursive=true)
      run_trixi(joinpath(examples_dir(), "t8code_3d_dgsem", "elixir_advection_cubed_sphere.jl"), maxiters=2)

      @timed_testset "mesh data" begin
        # create the output file to be tested
        @test_nowarn trixi2vtk(joinpath(outdir, "mesh.h5"), output_directory=outdir)
        outfilename = "mesh_celldata.vtu"
        out_file = joinpath(outdir, outfilename)

        # save output file to `artifacts` to facilitate debugging of failing tests
        testname = "3d-t8codemesh"
        cp(out_file, joinpath(artifacts_dir, testname * "-" * outfilename), force=true)

        # remote file path is actually a URL so it always has the same path structure
        remote_filename = "3d/t8codemesh/t8_dgsem_adv_sphere_mesh_02.vtu"
        ref_file = get_test_reference_file("t8_dgsem_adv_sphere_mesh_02.vtu", remote_filename)
        compare_cell_data(out_file, ref_file)
      end

      @timed_testset "solution celldata" begin
        # create the output file to be tested
        @test_nowarn trixi2vtk(joinpath(outdir, "solution_" * LEADING_ZEROS * "000002.h5"), output_directory=outdir)
        outfilename = "solution_" * LEADING_ZEROS * "000002_celldata.vtu"
        out_file = joinpath(outdir, outfilename)

        # save output file to `artifacts` to facilitate debugging of failing tests
        testname = "3d-t8codemesh"
        cp(out_file, joinpath(artifacts_dir, testname * "-" * outfilename), force=true)

        # remote file path is actually a URL so it always has the same path structure
        remote_filename = "3d/t8codemesh/t8_dgsem_adv_sphere_celldata_02.vtu"
        ref_file = get_test_reference_file("t8_dgsem_adv_sphere_celldata_02.vtu", remote_filename)
        compare_cell_data(out_file, ref_file)
      end

      @timed_testset "reinterpolate with nonuniform data" begin
        # Create and test output with reinterpolation (default options: `reinterpolate=true, data_is_uniform=false`)
        @test_nowarn trixi2vtk(joinpath(outdir, "solution_" * LEADING_ZEROS * "000002.h5"), output_directory=outdir)
        outfilename = "solution_" * LEADING_ZEROS * "000002.vtu"
        out_file = joinpath(outdir, outfilename)

        # save output file to `artifacts` to facilitate debugging of failing tests
        testname = "3d-t8codemesh-reinterp"
        cp(out_file, joinpath(artifacts_dir, testname * "-" * outfilename), force=true)

        # remote file path is actually a URL so it always has the same path structure
        remote_filename = "3d/t8codemesh/t8_dgsem_adv_sphere_reinterp_02.vtu"
        ref_file = get_test_reference_file("t8_dgsem_adv_sphere_reinterp_02.vtu", remote_filename)
        compare_point_data(out_file, ref_file)
      end

      @timed_testset "do not reinterpolate with nonuniform data" begin
        # Create and test output without reinterpolation on LGL nodes
        @test_nowarn trixi2vtk(joinpath(outdir, "solution_" * LEADING_ZEROS * "000002.h5"), output_directory=outdir, reinterpolate=false)
        outfilename = "solution_" * LEADING_ZEROS * "000002.vtu"
        out_file = joinpath(outdir, outfilename)

        # save output file to `artifacts` to facilitate debugging of failing tests
        testname = "3d-t8codemesh-no-reinterp"
        cp(out_file, joinpath(artifacts_dir, testname * "-" * outfilename), force=true)

        # remote file path is actually a URL so it always has the same path structure
        remote_filename = "3d/t8codemesh/t8_dgsem_adv_sphere_no_reinterp_02.vtu"
        ref_file = get_test_reference_file("t8_dgsem_adv_sphere_no_reinterp_02.vtu", remote_filename)
        compare_point_data(out_file, ref_file)
      end

      @timed_testset "do not reinterpolate with uniform data" begin
        # Create and test output without reinterpolation on uniform nodes
        # OBS! This is a dummy test just to exercise code. The resulting plot will look weird.
        @test_nowarn trixi2vtk(joinpath(outdir, "solution_" * LEADING_ZEROS * "000002.h5"), output_directory=outdir, reinterpolate=false, data_is_uniform=true)
        outfilename = "solution_" * LEADING_ZEROS * "000002.vtu"
        out_file = joinpath(outdir, outfilename)

        # save output file to `artifacts` to facilitate debugging of failing tests
        testname = "3d-t8codemesh-no-reinterp-uniform"
        cp(out_file, joinpath(artifacts_dir, testname * "-" * outfilename), force=true)

        # remote file path is actually a URL so it always has the same path structure
        remote_filename = "3d/t8codemesh/t8_dgsem_adv_sphere_no_reinterp_uniform_02.vtu"
        ref_file = get_test_reference_file("t8_dgsem_adv_sphere_no_reinterp_uniform_02.vtu", remote_filename)
        compare_point_data(out_file, ref_file)
      end
    end
  end
end

# Clean up afterwards: delete Trixi output directory and reference file directory
@test_nowarn rm(outdir, recursive=true)
@test_nowarn rm(TEST_REFERENCE_DIR, recursive=true)

end
