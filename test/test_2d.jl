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
    run_trixi(joinpath(examples_dir(), "tree_2d_dgsem", "elixir_euler_sedov_blast_wave.jl"), maxiters=10)

    @timed_testset "mesh data" begin
      # create the output file to be tested
      @test_nowarn trixi2vtk(joinpath(outdir, "mesh_" * LEADING_ZEROS * "000010.h5"), output_directory=outdir)
      outfilename = "mesh_" * LEADING_ZEROS * "000010_celldata.vtu"
      out_file = joinpath(outdir, outfilename)

      # save output file to `artifacts` to facilitate debugging of failing tests
      testname = "2d-treemesh"
      cp(out_file, joinpath(artifacts_dir, testname * "-" * outfilename), force=true)

      # remote file path is actually a URL so it always has the same path structure
      remote_filename = "2d/treemesh/dgsem_sedov_amr_mesh_10.vtu"
      ref_file = get_test_reference_file("dgsem_sedov_amr_mesh_10.vtu", remote_filename)
      compare_cell_data(out_file, ref_file)
    end

    @timed_testset "solution celldata" begin
      # create the output file to be tested
      @test_nowarn trixi2vtk(joinpath(outdir, "solution_" * LEADING_ZEROS * "000010.h5"), output_directory=outdir)
      outfilename = "solution_" * LEADING_ZEROS * "000010_celldata.vtu"
      out_file = joinpath(outdir, outfilename)

      # save output file to `artifacts` to facilitate debugging of failing tests
      testname = "2d-treemesh"
      cp(out_file, joinpath(artifacts_dir, testname * "-" * outfilename), force=true)

      # remote file path is actually a URL so it always has the same path structure
      remote_filename = "2d/treemesh/dgsem_sedov_amr_celldata_10.vtu"
      ref_file = get_test_reference_file("dgsem_sedov_amr_celldata_10.vtu", remote_filename)
      compare_cell_data(out_file, ref_file)
    end

    @timed_testset "reinterpolate with nonuniform data with VTU format" begin
      # Create and test output with reinterpolation (default options: `reinterpolate=true, data_is_uniform=false`)
      @test_nowarn trixi2vtk(joinpath(outdir, "solution_" * LEADING_ZEROS * "000010.h5"), output_directory=outdir)
      outfilename = "solution_" * LEADING_ZEROS * "000010.vtu"
      out_file = joinpath(outdir, outfilename)

      # save output file to `artifacts` to facilitate debugging of failing tests
      testname = "2d-treemesh-reinterp"
      cp(out_file, joinpath(artifacts_dir, testname * "-" * outfilename), force=true)

      # remote file path is actually a URL so it always has the same path structure
      remote_filename = "2d/treemesh/dgsem_sedov_amr_reinterp_10.vtu"
      ref_file = get_test_reference_file("dgsem_sedov_amr_reinterp_10.vtu", remote_filename)
      compare_cell_data(out_file, ref_file)
    end

    @timed_testset "reinterpolate with nonuniform data with VTI format" begin
      # Create and test output with reinterpolation (default options: `reinterpolate=true, data_is_uniform=false`)
      @test_nowarn trixi2vtk(joinpath(outdir, "solution_" * LEADING_ZEROS * "000010.h5"), output_directory=outdir, format=:vti)
      outfilename = "solution_" * LEADING_ZEROS * "000010.vti"
      out_file = joinpath(outdir, outfilename)

      # save output file to `artifacts` to facilitate debugging of failing tests
      testname = "2d-treemesh-vti-format"
      cp(out_file, joinpath(artifacts_dir, testname * "-" * outfilename), force=true)

      # remote file path is actually a URL so it always has the same path structure
      remote_filename = "2d/treemesh/dgsem_sedov_amr_10.vti"
      ref_file = get_test_reference_file("dgsem_sedov_amr_10.vti", remote_filename)
      compare_cell_data(out_file, ref_file)
    end

    @timed_testset "do not reinterpolate with nonuniform data" begin
      # Create and test output without reinterpolation on LGL nodes
      @test_nowarn trixi2vtk(joinpath(outdir, "solution_" * LEADING_ZEROS * "000010.h5"), output_directory=outdir, reinterpolate=false)
      outfilename = "solution_" * LEADING_ZEROS * "000010.vtu"
      out_file = joinpath(outdir, outfilename)

      # save output file to `artifacts` to facilitate debugging of failing tests
      testname = "2d-treemesh-no-reinterp"
      cp(out_file, joinpath(artifacts_dir, testname * "-" * outfilename), force=true)

      # remote file path is actually a URL so it always has the same path structure
      remote_filename = "2d/treemesh/dgsem_sedov_amr_no_reinterp_10.vtu"
      ref_file = get_test_reference_file("dgsem_sedov_amr_no_reinterp_10.vtu", remote_filename)
      compare_point_data(out_file, ref_file)
    end

    @timed_testset "do not reinterpolate with uniform data" begin
      # Create and test output without reinterpolation on uniform nodes
      # OBS! This is a dummy test just to exercise code. The resulting plot will look weird.
      @test_nowarn trixi2vtk(joinpath(outdir, "solution_" * LEADING_ZEROS * "000010.h5"), output_directory=outdir, reinterpolate=false, data_is_uniform=true)
      outfilename = "solution_" * LEADING_ZEROS * "000010.vtu"
      out_file = joinpath(outdir, outfilename)

      # save output file to `artifacts` to facilitate debugging of failing tests
      testname = "2d-treemesh-no-reinterp-uniform"
      cp(out_file, joinpath(artifacts_dir, testname * "-" * outfilename), force=true)

      # remote file path is actually a URL so it always has the same path structure
      remote_filename = "2d/treemesh/dgsem_sedov_amr_no_reinterp_uniform_10.vtu"
      ref_file = get_test_reference_file("dgsem_sedov_amr_no_reinterp_uniform_10.vtu", remote_filename)
      compare_point_data(out_file, ref_file)
    end

    @timed_testset "attempt reinterpolate with uniform data" begin
      # Purposely request a bad configuration and check that an error message gets thrown
      # OBS! Only needs tested once across all mesh types and dimensions
      @test_throws ArgumentError trixi2vtk(joinpath(outdir, "solution_" * LEADING_ZEROS * "000010.h5"), output_directory=outdir, data_is_uniform=true)
    end
  end

  if !Sys.iswindows() && get(ENV, "CI", nothing) == "true"
    # OBS! Only `TreeMesh` results are tested on Windows runners due to memory limits.
    #      All remaining mesh types are tested on Ubuntu and Mac
    @testset "StructuredMesh" begin
      isdir(outdir) && rm(outdir, recursive=true)
      run_trixi(joinpath(examples_dir(), "structured_2d_dgsem", "elixir_advection_waving_flag.jl"), maxiters=1)

      @timed_testset "mesh data" begin
        # create the output file to be tested
        @test_nowarn trixi2vtk(joinpath(outdir, "mesh.h5"), output_directory=outdir)
        outfilename = "mesh_celldata.vtu"
        out_file = joinpath(outdir, outfilename)

        # save output file to `artifacts` to facilitate debugging of failing tests
        testname = "2d-structuredmesh"
        cp(out_file, joinpath(artifacts_dir, testname * "-" * outfilename), force=true)

        # remote file path is actually a URL so it always has the same path structure
        remote_filename = "2d/structuredmesh/dgsem_adv_mesh_01.vtu"
        ref_file = get_test_reference_file("dgsem_adv_mesh_01.vtu", remote_filename)
        compare_cell_data(out_file, ref_file)
      end

      @timed_testset "solution celldata" begin
        # create the output file to be tested
        # OBS! This exercises passing multiple files (in this case 2 files) to `trixi2vtk`
        #      that only needs to be tested once.
        @test_nowarn trixi2vtk(joinpath(outdir, "solution_00000*.h5"), output_directory=outdir)

        outfilename = "solution_" * LEADING_ZEROS * "000001_celldata.vtu"
        out_file = joinpath(outdir, outfilename)

        # save output file to `artifacts` to facilitate debugging of failing tests
        testname = "2d-structuredmesh"
        cp(out_file, joinpath(artifacts_dir, testname * "-" * outfilename), force=true)

        # remote file path is actually a URL so it always has the same path structure
        remote_filename = "2d/structuredmesh/dgsem_adv_celldata_01.vtu"
        ref_file = get_test_reference_file("dgsem_adv_celldata_01.vtu", remote_filename)
        compare_cell_data(out_file, ref_file)
      end

      @timed_testset "reinterpolate with nonuniform data" begin
        # Create and test output with reinterpolation (default options: `reinterpolate=true, data_is_uniform=false`)
        @test_nowarn trixi2vtk(joinpath(outdir, "solution_" * LEADING_ZEROS * "000001.h5"), output_directory=outdir)
        outfilename = "solution_" * LEADING_ZEROS * "000001.vtu"
        out_file = joinpath(outdir, outfilename)

        # save output file to `artifacts` to facilitate debugging of failing tests
        testname = "2d-structuredmesh-reinterp"
        cp(out_file, joinpath(artifacts_dir, testname * "-" * outfilename), force=true)

        # remote file path is actually a URL so it always has the same path structure
        remote_filename = "2d/structuredmesh/dgsem_adv_reinterp_01.vtu"
        ref_file = get_test_reference_file("dgsem_adv_reinterp_01.vtu", remote_filename)
        compare_point_data(out_file, ref_file)
      end

      @timed_testset "do not reinterpolate with nonuniform data" begin
        # Create and test output without reinterpolation on LGL nodes
        @test_nowarn trixi2vtk(joinpath(outdir, "solution_" * LEADING_ZEROS * "000001.h5"), output_directory=outdir, reinterpolate=false)
        outfilename = "solution_" * LEADING_ZEROS * "000001.vtu"
        out_file = joinpath(outdir, outfilename)

        # save output file to `artifacts` to facilitate debugging of failing tests
        testname = "2d-structuredmesh-no-reinterp"
        cp(out_file, joinpath(artifacts_dir, testname * "-" * outfilename), force=true)

        # remote file path is actually a URL so it always has the same path structure
        remote_filename = "2d/structuredmesh/dgsem_adv_no_reinterp_01.vtu"
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
        testname = "2d-structuredmesh-no-reinterp-uniform"
        cp(out_file, joinpath(artifacts_dir, testname * "-" * outfilename), force=true)

        # remote file path is actually a URL so it always has the same path structure
        remote_filename = "2d/structuredmesh/dgsem_adv_no_reinterp_uniform_01.vtu"
        ref_file = get_test_reference_file("dgsem_adv_no_reinterp_uniform_01.vtu", remote_filename)
        compare_point_data(out_file, ref_file)
      end

      @timed_testset "attempt VTI format on unsupported mesh type" begin
        # Purposely request a bad configuration and check that an error message gets thrown
        # OBS! Only needs tested once across all mesh types and dimensions
        @test_throws ArgumentError trixi2vtk(joinpath(outdir, "solution_" * LEADING_ZEROS * "000001.h5"), output_directory=outdir, format=:vti)
      end
    end

    @testset "UnstructuredMesh2D" begin
      isdir(outdir) && rm(outdir, recursive=true)
      run_trixi(joinpath(examples_dir(), "unstructured_2d_dgsem", "elixir_euler_ec.jl"), maxiters=1)

      @timed_testset "mesh data" begin
        # create the output file to be tested
        @test_nowarn trixi2vtk(joinpath(outdir, "mesh.h5"), output_directory=outdir)
        outfilename = "mesh_celldata.vtu"
        out_file = joinpath(outdir, outfilename)

        # save output file to `artifacts` to facilitate debugging of failing tests
        testname = "2d-unstructuredmesh"
        cp(out_file, joinpath(artifacts_dir, testname * "-" * outfilename), force=true)

        # remote file path is actually a URL so it always has the same path structure
        remote_filename = "2d/unstructuredmesh/dgsem_euler_mesh_01.vtu"
        ref_file = get_test_reference_file("dgsem_euler_mesh_01.vtu", remote_filename)
        compare_cell_data(out_file, ref_file)
      end

      @timed_testset "solution celldata" begin
        # create the output file to be tested
        @test_nowarn trixi2vtk(joinpath(outdir, "solution_" * LEADING_ZEROS * "000001.h5"), output_directory=outdir)
        outfilename = "solution_" * LEADING_ZEROS * "000001_celldata.vtu"
        out_file = joinpath(outdir, outfilename)

        # save output file to `artifacts` to facilitate debugging of failing tests
        testname = "2d-unstructuredmesh"
        cp(out_file, joinpath(artifacts_dir, testname * "-" * outfilename), force=true)

        # remote file path is actually a URL so it always has the same path structure
        remote_filename = "2d/unstructuredmesh/dgsem_euler_celldata_01.vtu"
        ref_file = get_test_reference_file("dgsem_euler_celldata_01.vtu", remote_filename)
        compare_cell_data(out_file, ref_file)
      end

      @timed_testset "reinterpolate with nonuniform data" begin
        # Create and test output with reinterpolation (default options: `reinterpolate=true, data_is_uniform=false`)
        @test_nowarn trixi2vtk(joinpath(outdir, "solution_" * LEADING_ZEROS * "000001.h5"), output_directory=outdir)
        outfilename = "solution_" * LEADING_ZEROS * "000001.vtu"
        out_file = joinpath(outdir, outfilename)

        # save output file to `artifacts` to facilitate debugging of failing tests
        testname = "2d-unstructuredmesh-reinterp"
        cp(out_file, joinpath(artifacts_dir, testname * "-" * outfilename), force=true)

        # remote file path is actually a URL so it always has the same path structure
        remote_filename = "2d/unstructuredmesh/dgsem_euler_reinterp_01.vtu"
        ref_file = get_test_reference_file("dgsem_euler_reinterp_01.vtu", remote_filename)
        compare_point_data(out_file, ref_file)
      end

      @timed_testset "do not reinterpolate with nonuniform data" begin
        # Create and test output without reinterpolation on LGL nodes
        @test_nowarn trixi2vtk(joinpath(outdir, "solution_" * LEADING_ZEROS * "000001.h5"), output_directory=outdir, reinterpolate=false)
        outfilename = "solution_" * LEADING_ZEROS * "000001.vtu"
        out_file = joinpath(outdir, outfilename)

        # save output file to `artifacts` to facilitate debugging of failing tests
        testname = "2d-unstructuredmesh-no-reinterp"
        cp(out_file, joinpath(artifacts_dir, testname * "-" * outfilename), force=true)

        # remote file path is actually a URL so it always has the same path structure
        remote_filename = "2d/unstructuredmesh/dgsem_euler_no_reinterp_01.vtu"
        ref_file = get_test_reference_file("dgsem_euler_no_reinterp_01.vtu", remote_filename)
        compare_point_data(out_file, ref_file)
      end

      @timed_testset "do not reinterpolate with uniform data" begin
        # Create and test output without reinterpolation on uniform nodes
        # OBS! This is a dummy test just to exercise code. The resulting plot will look weird.
        @test_nowarn trixi2vtk(joinpath(outdir, "solution_" * LEADING_ZEROS * "000001.h5"), output_directory=outdir, reinterpolate=false, data_is_uniform=true)
        outfilename = "solution_" * LEADING_ZEROS * "000001.vtu"
        out_file = joinpath(outdir, outfilename)

        # save output file to `artifacts` to facilitate debugging of failing tests
        testname = "2d-unstructuredmesh-no-reinterp-uniform"
        cp(out_file, joinpath(artifacts_dir, testname * "-" * outfilename), force=true)

        # remote file path is actually a URL so it always has the same path structure
        remote_filename = "2d/unstructuredmesh/dgsem_euler_no_reinterp_uniform_01.vtu"
        ref_file = get_test_reference_file("dgsem_euler_no_reinterp_uniform_01.vtu", remote_filename)
        compare_point_data(out_file, ref_file)
      end
    end

    @testset "P4estMesh" begin
      isdir(outdir) && rm(outdir, recursive=true)
      run_trixi(joinpath(examples_dir(), "p4est_2d_dgsem", "elixir_mhd_rotor.jl"), maxiters=5)

      @timed_testset "mesh data" begin
        # create the output file to be tested
        @test_nowarn trixi2vtk(joinpath(outdir, "mesh_" * LEADING_ZEROS * "000005.h5"), output_directory=outdir)
        outfilename = "mesh_" * LEADING_ZEROS * "000005_celldata.vtu"
        out_file = joinpath(outdir, outfilename)

        # save output file to `artifacts` to facilitate debugging of failing tests
        testname = "2d-p4estmesh"
        cp(out_file, joinpath(artifacts_dir, testname * "-" * outfilename), force=true)

        # remote file path is actually a URL so it always has the same path structure
        remote_filename = "2d/p4estmesh/dgsem_rotor_amr_mesh_05.vtu"
        ref_file = get_test_reference_file("dgsem_rotor_amr_mesh_05.vtu", remote_filename)
        compare_cell_data(out_file, ref_file)
      end

      @timed_testset "solution celldata" begin
        # create the output file to be tested
        @test_nowarn trixi2vtk(joinpath(outdir, "solution_" * LEADING_ZEROS * "000005.h5"), output_directory=outdir)
        outfilename = "solution_" * LEADING_ZEROS * "000005_celldata.vtu"
        out_file = joinpath(outdir, outfilename)

        # save output file to `artifacts` to facilitate debugging of failing tests
        testname = "2d-p4estmesh"
        cp(out_file, joinpath(artifacts_dir, testname * "-" * outfilename), force=true)

        # remote file path is actually a URL so it always has the same path structure
        remote_filename = "2d/p4estmesh/dgsem_rotor_amr_celldata_05.vtu"
        ref_file = get_test_reference_file("dgsem_rotor_amr_celldata_05.vtu", remote_filename)
        compare_cell_data(out_file, ref_file)
      end

      @timed_testset "reinterpolate with nonuniform data" begin
        # Create and test output with reinterpolation (default options: `reinterpolate=true, data_is_uniform=false`)
        @test_nowarn trixi2vtk(joinpath(outdir, "solution_" * LEADING_ZEROS * "000005.h5"), output_directory=outdir)
        outfilename = "solution_" * LEADING_ZEROS * "000005.vtu"
        out_file = joinpath(outdir, outfilename)

        # save output file to `artifacts` to facilitate debugging of failing tests
        testname = "2d-p4estmesh-reinterp"
        cp(out_file, joinpath(artifacts_dir, testname * "-" * outfilename), force=true)

        # remote file path is actually a URL so it always has the same path structure
        remote_filename = "2d/p4estmesh/dgsem_rotor_amr_reinterp_05.vtu"
        ref_file = get_test_reference_file("dgsem_rotor_amr_reinterp_05.vtu", remote_filename)
        compare_point_data(out_file, ref_file)
      end

      @timed_testset "do not reinterpolate with nonuniform data" begin
        # Create and test output without reinterpolation on LGL nodes
        @test_nowarn trixi2vtk(joinpath(outdir, "solution_" * LEADING_ZEROS * "000005.h5"), output_directory=outdir, reinterpolate=false)
        outfilename = "solution_" * LEADING_ZEROS * "000005.vtu"
        out_file = joinpath(outdir, outfilename)

        # save output file to `artifacts` to facilitate debugging of failing tests
        testname = "2d-p4estmesh-no-reinterp"
        cp(out_file, joinpath(artifacts_dir, testname * "-" * outfilename), force=true)

        # remote file path is actually a URL so it always has the same path structure
        remote_filename = "2d/p4estmesh/dgsem_rotor_amr_no_reinterp_05.vtu"
        ref_file = get_test_reference_file("dgsem_rotor_amr_no_reinterp_05.vtu", remote_filename)
        compare_point_data(out_file, ref_file)
      end

      @timed_testset "do not reinterpolate with uniform data" begin
        # Create and test output without reinterpolation on uniform nodes
        # OBS! This is a dummy test just to exercise code. The resulting plot will look weird.
        @test_nowarn trixi2vtk(joinpath(outdir, "solution_" * LEADING_ZEROS * "000005.h5"), output_directory=outdir, reinterpolate=false, data_is_uniform=true)
        outfilename = "solution_" * LEADING_ZEROS * "000005.vtu"
        out_file = joinpath(outdir, outfilename)

        # save output file to `artifacts` to facilitate debugging of failing tests
        testname = "2d-p4estmesh-no-reinterp-uniform"
        cp(out_file, joinpath(artifacts_dir, testname * "-" * outfilename), force=true)

        # remote file path is actually a URL so it always has the same path structure
        remote_filename = "2d/p4estmesh/dgsem_rotor_amr_no_reinterp_uniform_05.vtu"
        ref_file = get_test_reference_file("dgsem_rotor_amr_no_reinterp_uniform_05.vtu", remote_filename)
        compare_point_data(out_file, ref_file)
      end
    end

    @testset "T8codeMesh" begin
      isdir(outdir) && rm(outdir, recursive=true)
      run_trixi(joinpath(examples_dir(), "t8code_2d_dgsem", "elixir_mhd_rotor.jl"), maxiters=5)

      @timed_testset "mesh data" begin
        # create the output file to be tested
        @test_nowarn trixi2vtk(joinpath(outdir, "mesh_" * LEADING_ZEROS * "000005.h5"), output_directory=outdir)
        outfilename = "mesh_" * LEADING_ZEROS * "000005_celldata.vtu"
        out_file = joinpath(outdir, outfilename)

        # save output file to `artifacts` to facilitate debugging of failing tests
        testname = "2d-t8codemesh"
        cp(out_file, joinpath(artifacts_dir, testname * "-" * outfilename), force=true)

        # remote file path is actually a URL so it always has the same path structure
        remote_filename = "2d/p4estmesh/dgsem_rotor_amr_mesh_05.vtu"
        ref_file = get_test_reference_file("dgsem_rotor_amr_mesh_05.vtu", remote_filename)
        compare_cell_data(out_file, ref_file)
      end

      @timed_testset "solution celldata" begin
        # create the output file to be tested
        @test_nowarn trixi2vtk(joinpath(outdir, "solution_" * LEADING_ZEROS * "000005.h5"), output_directory=outdir)
        outfilename = "solution_" * LEADING_ZEROS * "000005_celldata.vtu"
        out_file = joinpath(outdir, outfilename)

        # save output file to `artifacts` to facilitate debugging of failing tests
        testname = "2d-t8codemesh"
        cp(out_file, joinpath(artifacts_dir, testname * "-" * outfilename), force=true)

        # remote file path is actually a URL so it always has the same path structure
        remote_filename = "2d/p4estmesh/dgsem_rotor_amr_celldata_05.vtu"
        ref_file = get_test_reference_file("dgsem_rotor_amr_celldata_05.vtu", remote_filename)
        compare_cell_data(out_file, ref_file)
      end

      @timed_testset "reinterpolate with nonuniform data" begin
        # Create and test output with reinterpolation (default options: `reinterpolate=true, data_is_uniform=false`)
        @test_nowarn trixi2vtk(joinpath(outdir, "solution_" * LEADING_ZEROS * "000005.h5"), output_directory=outdir)
        outfilename = "solution_" * LEADING_ZEROS * "000005.vtu"
        out_file = joinpath(outdir, outfilename)

        # save output file to `artifacts` to facilitate debugging of failing tests
        testname = "2d-t8codemesh-reinterp"
        cp(out_file, joinpath(artifacts_dir, testname * "-" * outfilename), force=true)

        # remote file path is actually a URL so it always has the same path structure
        remote_filename = "2d/p4estmesh/dgsem_rotor_amr_reinterp_05.vtu"
        ref_file = get_test_reference_file("dgsem_rotor_amr_reinterp_05.vtu", remote_filename)
        compare_point_data(out_file, ref_file)
      end

      @timed_testset "do not reinterpolate with nonuniform data" begin
        # Create and test output without reinterpolation on LGL nodes
        @test_nowarn trixi2vtk(joinpath(outdir, "solution_" * LEADING_ZEROS * "000005.h5"), output_directory=outdir, reinterpolate=false)
        outfilename = "solution_" * LEADING_ZEROS * "000005.vtu"
        out_file = joinpath(outdir, outfilename)

        # save output file to `artifacts` to facilitate debugging of failing tests
        testname = "2d-t8codemesh-no-reinterp"
        cp(out_file, joinpath(artifacts_dir, testname * "-" * outfilename), force=true)

        # remote file path is actually a URL so it always has the same path structure
        remote_filename = "2d/p4estmesh/dgsem_rotor_amr_no_reinterp_05.vtu"
        ref_file = get_test_reference_file("dgsem_rotor_amr_no_reinterp_05.vtu", remote_filename)
        compare_point_data(out_file, ref_file)
      end

      @timed_testset "do not reinterpolate with uniform data" begin
        # Create and test output without reinterpolation on uniform nodes
        # OBS! This is a dummy test just to exercise code. The resulting plot will look weird.
        @test_nowarn trixi2vtk(joinpath(outdir, "solution_" * LEADING_ZEROS * "000005.h5"), output_directory=outdir, reinterpolate=false, data_is_uniform=true)
        outfilename = "solution_" * LEADING_ZEROS * "000005.vtu"
        out_file = joinpath(outdir, outfilename)

        # save output file to `artifacts` to facilitate debugging of failing tests
        testname = "2d-t8codemesh-no-reinterp-uniform"
        cp(out_file, joinpath(artifacts_dir, testname * "-" * outfilename), force=true)

        # remote file path is actually a URL so it always has the same path structure
        remote_filename = "2d/p4estmesh/dgsem_rotor_amr_no_reinterp_uniform_05.vtu"
        ref_file = get_test_reference_file("dgsem_rotor_amr_no_reinterp_uniform_05.vtu", remote_filename)
        compare_point_data(out_file, ref_file)
      end
    end
  end

  @testset "Subcell limiting coefficients" begin
    isdir(outdir) && rm(outdir, recursive=true)
    run_trixi(joinpath(examples_dir(), "tree_2d_dgsem", "elixir_euler_sedov_blast_wave_sc_subcell.jl"),
              maxiters=10, initial_refinement_level=4)

    @timed_testset "without reinterpolation" begin
      # Create and test output without reinterpolation
      @test_nowarn trixi2vtk(joinpath(outdir, "solution_" * LEADING_ZEROS * "000010.h5"), output_directory=outdir, reinterpolate=false)
      outfilename = "solution_" * LEADING_ZEROS * "000010.vtu"
      out_file = joinpath(outdir, outfilename)

      # save output file to `artifacts` to facilitate debugging of failing tests
      testname = "2d-treemesh-shockcapturing-subcell-no-reinterp"
      cp(out_file, joinpath(artifacts_dir, testname * "-" * outfilename), force=true)

      # remote file path is actually a URL so it always has the same path structure
      remote_filename = "2d/treemesh/dgsem_sedov_subcell_no_interp_10.vtu"
      ref_file = get_test_reference_file("dgsem_sedov_subcell_no_interp_10.vtu", remote_filename)
      compare_point_data(out_file, ref_file)
    end

    @timed_testset "with reinterpolation" begin
      # Create and test output without reinterpolation
      @test_nowarn trixi2vtk(joinpath(outdir, "solution_" * LEADING_ZEROS * "000010.h5"), output_directory=outdir, reinterpolate=true)
      outfilename = "solution_" * LEADING_ZEROS * "000010.vtu"
      out_file = joinpath(outdir, outfilename)

      # save output file to `artifacts` to facilitate debugging of failing tests
      testname = "2d-treemesh-shockcapturing-subcell-reinterp"
      cp(out_file, joinpath(artifacts_dir, testname * "-" * outfilename), force=true)

      # remote file path is actually a URL so it always has the same path structure
      remote_filename = "2d/treemesh/dgsem_sedov_subcell_interp_10.vtu"
      ref_file = get_test_reference_file("dgsem_sedov_subcell_interp_10.vtu", remote_filename)
      compare_cell_data(out_file, ref_file)
    end
  end

  @testset "2D surface in 3D using P4estMesh" begin
    isdir(outdir) && rm(outdir, recursive=true)

    @timed_testset "without reinterpolation" begin
      # Get .h5 files from remote repository (in the future, we will run TrixiAtmo.jl here)
      solution_file = get_test_reference_file("dgsem_isolated_mountain_05.h5",
                                              "2d/p4estmesh_surface/dgsem_isolated_mountain_05.h5")
      mesh_file = get_test_reference_file("mesh.h5", "2d/p4estmesh_surface/mesh.h5")
      p4est_data_file = get_test_reference_file("p4est_data",
                                                "2d/p4estmesh_surface/p4est_data")

      # Create and test output without reinterpolation
      @test_nowarn trixi2vtk(solution_file, output_directory=outdir, reinterpolate=false)
      outfilename = "dgsem_isolated_mountain_05.vtu"
      out_file = joinpath(outdir, outfilename)

      # save output file to `artifacts` to facilitate debugging of failing tests
      testname = "2d-p4estmesh-2d-surface-no-reinterp"
      cp(out_file, joinpath(artifacts_dir, testname * "-" * outfilename), force=true)

      # remote file path is actually a URL so it always has the same path structure
      remote_filename = "2d/p4estmesh_surface/dgsem_isolated_mountain_no_reinterp_05.vtu"
      ref_file = get_test_reference_file("dgsem_isolated_mountain_05_no_reinterp.vtu",
                                         remote_filename)
      compare_point_data(out_file, ref_file)
    end

    @timed_testset "with reinterpolation" begin
      # Get .h5 files from remote repository (in the future, we will run TrixiAtmo.jl here)
      solution_file = get_test_reference_file("dgsem_isolated_mountain_05.h5",
                                              "2d/p4estmesh_surface/dgsem_isolated_mountain_05.h5")
      mesh_file = get_test_reference_file("mesh.h5", "2d/p4estmesh_surface/mesh.h5")
      p4est_data_file = get_test_reference_file("p4est_data",
                                                "2d/p4estmesh_surface/p4est_data")

      # Create and test output without reinterpolation
      @test_nowarn trixi2vtk(solution_file, output_directory=outdir, reinterpolate=true)
      outfilename = "dgsem_isolated_mountain_05.vtu"
      out_file = joinpath(outdir, outfilename)

      # save output file to `artifacts` to facilitate debugging of failing tests
      testname = "2d-p4estmesh-2d-surface-no-reinterp"
      cp(out_file, joinpath(artifacts_dir, testname * "-" * outfilename), force=true)

      # remote file path is actually a URL so it always has the same path structure
      remote_filename = "2d/p4estmesh_surface/dgsem_isolated_mountain_reinterp_05.vtu"
      ref_file = get_test_reference_file("dgsem_isolated_mountain_05_reinterp.vtu",
                                        remote_filename)
      compare_point_data(out_file, ref_file)
    end
  end

end

# Clean up afterwards: delete Trixi output directory and reference file directory
@test_nowarn rm(outdir, recursive=true)
@test_nowarn rm(TEST_REFERENCE_DIR, recursive=true)

end
