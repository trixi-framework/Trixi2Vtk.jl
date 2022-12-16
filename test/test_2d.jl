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
      @test_nowarn trixi2vtk(joinpath(outdir, "mesh_000010.h5"), output_directory=outdir)
      outfilename = "mesh_000010_celldata.vtu"
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
      @test_nowarn trixi2vtk(joinpath(outdir, "solution_000010.h5"), output_directory=outdir)
      outfilename = "solution_000010_celldata.vtu"
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
      @test_nowarn trixi2vtk(joinpath(outdir, "solution_000010.h5"), output_directory=outdir)
      outfilename = "solution_000010.vtu"
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
      @test_nowarn trixi2vtk(joinpath(outdir, "solution_000010.h5"), output_directory=outdir, format=:vti)
      outfilename = "solution_000010.vti"
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
      @test_nowarn trixi2vtk(joinpath(outdir, "solution_000010.h5"), output_directory=outdir, reinterpolate=false)
      outfilename = "solution_000010.vtu"
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
      @test_nowarn trixi2vtk(joinpath(outdir, "solution_000010.h5"), output_directory=outdir, reinterpolate=false, data_is_uniform=true)
      outfilename = "solution_000010.vtu"
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
      @test_throws ArgumentError trixi2vtk(joinpath(outdir, "solution_000010.h5"), output_directory=outdir, data_is_uniform=true)
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

        outfilename = "solution_000001_celldata.vtu"
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
        @test_nowarn trixi2vtk(joinpath(outdir, "solution_000001.h5"), output_directory=outdir)
        outfilename = "solution_000001.vtu"
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
        @test_nowarn trixi2vtk(joinpath(outdir, "solution_000001.h5"), output_directory=outdir, reinterpolate=false)
        outfilename = "solution_000001.vtu"
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
        @test_nowarn trixi2vtk(joinpath(outdir, "solution_000001.h5"), output_directory=outdir, reinterpolate=false, data_is_uniform=true)
        outfilename = "solution_000001.vtu"
        out_file = joinpath(outdir, outfilename)

        # save output file to `artifacts` to facilitate debugging of failing tests
        testname = "2d-structuredmesh-no-reinterp-uniform"
        cp(out_file, joinpath(artifacts_dir, testname * "-" * outfilename), force=true)

        # remote file path is actually a URL so it always has the same path structure
        remote_filename = "2d/structuredmesh/dgsem_adv_no_reinterp_uniform_01.vtu"
        ref_file = get_test_reference_file("dgsem_adv_no_reinterp_uniform_01.vtu", remote_filename)
        compare_point_data(out_file, ref_file)
      end

      @timed_testset "attempt VTI format on unsupportd mesh type" begin
        # Purposely request a bad configuration and check that an error message gets thrown
        # OBS! Only needs tested once across all mesh types and dimensions
        @test_throws ArgumentError trixi2vtk(joinpath(outdir, "solution_000001.h5"), output_directory=outdir, format=:vti)
      end
    end

    @testset "UnstructuredMesh2D" begin
      isdir(outdir) && rm(outdir, recursive=true)
      run_trixi(joinpath(examples_dir(), "unstructured_2d_dgsem", "elixir_shallowwater_ec.jl"), maxiters=1)

      @timed_testset "mesh data" begin
        # create the output file to be tested
        @test_nowarn trixi2vtk(joinpath(outdir, "mesh.h5"), output_directory=outdir)
        outfilename = "mesh_celldata.vtu"
        out_file = joinpath(outdir, outfilename)

        # save output file to `artifacts` to facilitate debugging of failing tests
        testname = "2d-unstructuredmesh"
        cp(out_file, joinpath(artifacts_dir, testname * "-" * outfilename), force=true)

        # remote file path is actually a URL so it always has the same path structure
        remote_filename = "2d/unstructuredmesh/dgsem_swe_mesh_01.vtu"
        ref_file = get_test_reference_file("dgsem_swe_mesh_01.vtu", remote_filename)
        compare_cell_data(out_file, ref_file)
      end

      @timed_testset "solution celldata" begin
        # create the output file to be tested
        @test_nowarn trixi2vtk(joinpath(outdir, "solution_000001.h5"), output_directory=outdir)
        outfilename = "solution_000001_celldata.vtu"
        out_file = joinpath(outdir, outfilename)

        # save output file to `artifacts` to facilitate debugging of failing tests
        testname = "2d-unstructuredmesh"
        cp(out_file, joinpath(artifacts_dir, testname * "-" * outfilename), force=true)

        # remote file path is actually a URL so it always has the same path structure
        remote_filename = "2d/unstructuredmesh/dgsem_swe_celldata_01.vtu"
        ref_file = get_test_reference_file("dgsem_swe_celldata_01.vtu", remote_filename)
        compare_cell_data(out_file, ref_file)
      end

      @timed_testset "reinterpolate with nonuniform data" begin
        # Create and test output with reinterpolation (default options: `reinterpolate=true, data_is_uniform=false`)
        @test_nowarn trixi2vtk(joinpath(outdir, "solution_000001.h5"), output_directory=outdir)
        outfilename = "solution_000001.vtu"
        out_file = joinpath(outdir, outfilename)

        # save output file to `artifacts` to facilitate debugging of failing tests
        testname = "2d-unstructuredmesh-reinterp"
        cp(out_file, joinpath(artifacts_dir, testname * "-" * outfilename), force=true)

        # remote file path is actually a URL so it always has the same path structure
        remote_filename = "2d/unstructuredmesh/dgsem_swe_reinterp_01.vtu"
        ref_file = get_test_reference_file("dgsem_swe_reinterp_01.vtu", remote_filename)
        compare_point_data(out_file, ref_file)
      end

      @timed_testset "do not reinterpolate with nonuniform data" begin
        # Create and test output without reinterpolation on LGL nodes
        @test_nowarn trixi2vtk(joinpath(outdir, "solution_000001.h5"), output_directory=outdir, reinterpolate=false)
        outfilename = "solution_000001.vtu"
        out_file = joinpath(outdir, outfilename)

        # save output file to `artifacts` to facilitate debugging of failing tests
        testname = "2d-unstructuredmesh-no-reinterp"
        cp(out_file, joinpath(artifacts_dir, testname * "-" * outfilename), force=true)

        # remote file path is actually a URL so it always has the same path structure
        remote_filename = "2d/unstructuredmesh/dgsem_swe_no_reinterp_01.vtu"
        ref_file = get_test_reference_file("dgsem_swe_no_reinterp_01.vtu", remote_filename)
        compare_point_data(out_file, ref_file)
      end

      @timed_testset "do not reinterpolate with uniform data" begin
        # Create and test output without reinterpolation on uniform nodes
        # OBS! This is a dummy test just to exercise code. The resulting plot will look weird.
        @test_nowarn trixi2vtk(joinpath(outdir, "solution_000001.h5"), output_directory=outdir, reinterpolate=false, data_is_uniform=true)
        outfilename = "solution_000001.vtu"
        out_file = joinpath(outdir, outfilename)

        # save output file to `artifacts` to facilitate debugging of failing tests
        testname = "2d-unstructuredmesh-no-reinterp-uniform"
        cp(out_file, joinpath(artifacts_dir, testname * "-" * outfilename), force=true)

        # remote file path is actually a URL so it always has the same path structure
        remote_filename = "2d/unstructuredmesh/dgsem_swe_no_reinterp_uniform_01.vtu"
        ref_file = get_test_reference_file("dgsem_swe_no_reinterp_uniform_01.vtu", remote_filename)
        compare_point_data(out_file, ref_file)
      end
    end

    @testset "P4estMesh" begin
      isdir(outdir) && rm(outdir, recursive=true)
      run_trixi(joinpath(examples_dir(), "p4est_2d_dgsem", "elixir_mhd_rotor.jl"), maxiters=5)

      @timed_testset "mesh data" begin
        # create the output file to be tested
        @test_nowarn trixi2vtk(joinpath(outdir, "mesh_000005.h5"), output_directory=outdir)
        outfilename = "mesh_000005_celldata.vtu"
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
        @test_nowarn trixi2vtk(joinpath(outdir, "solution_000005.h5"), output_directory=outdir)
        outfilename = "solution_000005_celldata.vtu"
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
        @test_nowarn trixi2vtk(joinpath(outdir, "solution_000005.h5"), output_directory=outdir)
        outfilename = "solution_000005.vtu"
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
        @test_nowarn trixi2vtk(joinpath(outdir, "solution_000005.h5"), output_directory=outdir, reinterpolate=false)
        outfilename = "solution_000005.vtu"
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
        @test_nowarn trixi2vtk(joinpath(outdir, "solution_000005.h5"), output_directory=outdir, reinterpolate=false, data_is_uniform=true)
        outfilename = "solution_000005.vtu"
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
  end
end

# Clean up afterwards: delete Trixi output directory and reference file directory
@test_nowarn rm(outdir, recursive=true)
@test_nowarn rm(TEST_REFERENCE_DIR, recursive=true)

end
