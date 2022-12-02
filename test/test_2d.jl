module Test2D

using Test
using Trixi2Vtk

include("test_trixi2vtk.jl")

# Start with a clean environment: remove Trixi output directory if it exists
outdir = "out"
isdir(outdir) && rm(outdir, recursive=true)


@testset "2D" begin
  @testset "TreeMesh" begin
    isdir(outdir) && rm(outdir, recursive=true)
    run_trixi(joinpath(examples_dir(), "tree_2d_dgsem", "elixir_euler_sedov_blast_wave.jl"), maxiters=10)

    @timed_testset "mesh data" begin
      @test_nowarn trixi2vtk(joinpath(outdir, "mesh_000010.h5"), output_directory=outdir)
      out_file = joinpath(outdir,"mesh_000010_celldata.vtu")
      ref_file = get_test_reference_file("dgsem_sedov_amr_mesh_10.vtu",
                                         joinpath("2d", "treemesh", "dgsem_sedov_amr_mesh_10.vtu"))
      compare_cell_info(out_file, ref_file)
    end

    @timed_testset "solution celldata" begin
      @test_nowarn trixi2vtk(joinpath(outdir, "solution_000010.h5"), output_directory=outdir)
      out_file = joinpath(outdir,"solution_000010_celldata.vtu")
      ref_file = get_test_reference_file("dgsem_sedov_amr_celldata_10.vtu",
                                         joinpath("2d", "treemesh", "dgsem_sedov_amr_celldata_10.vtu"))
      compare_cell_info(out_file, ref_file)
    end

    @timed_testset "reinterpolate with nonuniform data" begin
      # Create and test output with reinterpolation (default options: `reinterpolate=true, data_is_uniform=false`)
      @test_nowarn trixi2vtk(joinpath(outdir, "solution_000010.h5"), output_directory=outdir)
      out_file = joinpath(outdir,"solution_000010.vtu")
      ref_file = get_test_reference_file("dgsem_sedov_amr_reinterp_10.vtu",
                                         joinpath("2d", "treemesh", "dgsem_sedov_amr_reinterp_10.vtu"))
      compare_cell_info(out_file, ref_file)
    end

    @timed_testset "do not reinterpolate with nonuniform data" begin
      # Create and test output without reinterpolation on LGL nodes
      @test_nowarn trixi2vtk(joinpath(outdir, "solution_000010.h5"), output_directory=outdir, reinterpolate=false)
      out_file = joinpath(outdir,"solution_000010.vtu")
      ref_file = get_test_reference_file("dgsem_sedov_amr_no_reinterp_10.vtu",
                                         joinpath("2d", "treemesh", "dgsem_sedov_amr_no_reinterp_10.vtu"))
      compare_point_info(out_file, ref_file)
    end

    @timed_testset "do not reinterpolate with uniform data" begin
      # Create and test output without reinterpolation on uniform nodes
      # OBS! This is a dummy test just to exercise code. The resulting plot will look weird.
      @test_nowarn trixi2vtk(joinpath(outdir, "solution_000010.h5"), output_directory=outdir, reinterpolate=false, data_is_uniform=true)
      out_file = joinpath(outdir,"solution_000010.vtu")
      ref_file = get_test_reference_file("dgsem_sedov_amr_no_reinterp_uniform_10.vtu",
                                         joinpath("2d", "treemesh", "dgsem_sedov_amr_no_reinterp_uniform_10.vtu"))
      compare_point_info(out_file, ref_file)
    end

    @timed_testset "attempt reinterpolate with uniform data" begin
      # Purposely request a bad configuration and check that an error message gets thrown
      # OBS! Only needs tested once across all mesh types and dimensions
      let err = nothing
        try
          trixi2vtk(joinpath(outdir, "solution_000010.h5"), output_directory=outdir, data_is_uniform=true)
        catch err
        end
        @test err isa Exception
        @test sprint(showerror, err) == "uniform data should not be reinterpolated! Set reinterpolate=false and try again."
      end
    end
  end

  @testset "StructuredMesh" begin
    isdir(outdir) && rm(outdir, recursive=true)
    run_trixi(joinpath(examples_dir(), "structured_2d_dgsem", "elixir_advection_waving_flag.jl"), maxiters=1)

    @timed_testset "mesh data" begin
      @test_nowarn trixi2vtk(joinpath(outdir, "mesh.h5"), output_directory=outdir)
      out_file = joinpath(outdir,"mesh_celldata.vtu")
      ref_file = get_test_reference_file("dgsem_adv_mesh_01.vtu",
                                         joinpath("2d", "structuredmesh", "dgsem_adv_mesh_01.vtu"))
      compare_cell_info(out_file, ref_file)
    end

    @timed_testset "solution celldata" begin
      @test_nowarn trixi2vtk(joinpath(outdir, "solution_000001.h5"), output_directory=outdir)
      out_file = joinpath(outdir,"solution_000001_celldata.vtu")
      ref_file = get_test_reference_file("dgsem_adv_celldata_01.vtu",
                                         joinpath("2d", "structuredmesh", "dgsem_adv_celldata_01.vtu"))
      compare_cell_info(out_file, ref_file)
    end

    @timed_testset "reinterpolate with nonuniform data" begin
      # Create and test output with reinterpolation (default options: `reinterpolate=true, data_is_uniform=false`)
      @test_nowarn trixi2vtk(joinpath(outdir, "solution_000001.h5"), output_directory=outdir)
      out_file = joinpath(outdir,"solution_000001.vtu")
      ref_file = get_test_reference_file("dgsem_adv_reinterp_01.vtu",
                                         joinpath("2d", "structuredmesh", "dgsem_adv_reinterp_01.vtu"))
      compare_point_info(out_file, ref_file)
    end

    @timed_testset "do not reinterpolate with nonuniform data" begin
      # Create and test output without reinterpolation on LGL nodes
      @test_nowarn trixi2vtk(joinpath(outdir, "solution_000001.h5"), output_directory=outdir, reinterpolate=false)
      out_file = joinpath(outdir,"solution_000001.vtu")
      ref_file = get_test_reference_file("dgsem_adv_no_reinterp_01.vtu",
                                         joinpath("2d", "structuredmesh", "dgsem_adv_no_reinterp_01.vtu"))
      compare_point_info(out_file, ref_file)
    end

    @timed_testset "do not reinterpolate with uniform data" begin
      # Create and test output without reinterpolation on uniform nodes
      # OBS! This is a dummy test just to exercise code. The resulting plot will look weird.
      @test_nowarn trixi2vtk(joinpath(outdir, "solution_000001.h5"), output_directory=outdir, reinterpolate=false, data_is_uniform=true)
      out_file = joinpath(outdir,"solution_000001.vtu")
      ref_file = get_test_reference_file("dgsem_adv_no_reinterp_uniform_01.vtu",
                                         joinpath("2d", "structuredmesh", "dgsem_adv_no_reinterp_uniform_01.vtu"))
      compare_point_info(out_file, ref_file)
    end
  end

  @testset "UnstructuredMesh2D" begin
    isdir(outdir) && rm(outdir, recursive=true)
    run_trixi(joinpath(examples_dir(), "unstructured_2d_dgsem", "elixir_shallowwater_ec.jl"), maxiters=1)

    @timed_testset "mesh data" begin
      @test_nowarn trixi2vtk(joinpath(outdir, "mesh.h5"), output_directory=outdir)
      out_file = joinpath(outdir,"mesh_celldata.vtu")
      ref_file = get_test_reference_file("dgsem_swe_mesh_01.vtu",
                                         joinpath("2d", "unstructuredmesh", "dgsem_swe_mesh_01.vtu"))
      compare_cell_info(out_file, ref_file)
    end

    @timed_testset "solution celldata" begin
      @test_nowarn trixi2vtk(joinpath(outdir, "solution_000001.h5"), output_directory=outdir)
      out_file = joinpath(outdir,"solution_000001_celldata.vtu")
      ref_file = get_test_reference_file("dgsem_swe_celldata_01.vtu",
                                         joinpath("2d", "unstructuredmesh", "dgsem_swe_celldata_01.vtu"))
      compare_cell_info(out_file, ref_file)
    end

    @timed_testset "reinterpolate with nonuniform data" begin
      # Create and test output with reinterpolation (default options: `reinterpolate=true, data_is_uniform=false`)
      @test_nowarn trixi2vtk(joinpath(outdir, "solution_000001.h5"), output_directory=outdir)
      out_file = joinpath(outdir,"solution_000001.vtu")
      ref_file = get_test_reference_file("dgsem_swe_reinterp_01.vtu",
                                         joinpath("2d", "unstructuredmesh", "dgsem_swe_reinterp_01.vtu"))
      compare_point_info(out_file, ref_file)
    end

    @timed_testset "do not reinterpolate with nonuniform data" begin
      # Create and test output without reinterpolation on LGL nodes
      @test_nowarn trixi2vtk(joinpath(outdir, "solution_000001.h5"), output_directory=outdir, reinterpolate=false)
      out_file = joinpath(outdir,"solution_000001.vtu")
      ref_file = get_test_reference_file("dgsem_swe_no_reinterp_01.vtu",
                                         joinpath("2d", "unstructuredmesh", "dgsem_swe_no_reinterp_01.vtu"))
      compare_point_info(out_file, ref_file)
    end

    @timed_testset "do not reinterpolate with uniform data" begin
      # Create and test output without reinterpolation on uniform nodes
      # OBS! This is a dummy test just to exercise code. The resulting plot will look weird.
      @test_nowarn trixi2vtk(joinpath(outdir, "solution_000001.h5"), output_directory=outdir, reinterpolate=false, data_is_uniform=true)
      out_file = joinpath(outdir,"solution_000001.vtu")
      ref_file = get_test_reference_file("dgsem_swe_no_reinterp_uniform_01.vtu",
                                         joinpath("2d", "unstructuredmesh", "dgsem_swe_no_reinterp_uniform_01.vtu"))
      compare_point_info(out_file, ref_file)
    end
  end

  @testset "P4estMesh" begin
    isdir(outdir) && rm(outdir, recursive=true)
    run_trixi(joinpath(examples_dir(), "p4est_2d_dgsem", "elixir_mhd_rotor.jl"), maxiters=5)

    @timed_testset "mesh data" begin
      @test_nowarn trixi2vtk(joinpath(outdir, "mesh_000005.h5"), output_directory=outdir)
      out_file = joinpath(outdir,"mesh_000005_celldata.vtu")
      ref_file = get_test_reference_file("dgsem_rotor_amr_mesh_05.vtu",
                                         joinpath("2d", "p4estmesh", "dgsem_rotor_amr_mesh_05.vtu"))
      compare_cell_info(out_file, ref_file)
    end

    @timed_testset "solution celldata" begin
      @test_nowarn trixi2vtk(joinpath(outdir, "solution_000005.h5"), output_directory=outdir)
      out_file = joinpath(outdir,"solution_000005_celldata.vtu")
      ref_file = get_test_reference_file("dgsem_rotor_amr_celldata_05.vtu",
                                         joinpath("2d", "p4estmesh", "dgsem_rotor_amr_celldata_05.vtu"))
      compare_cell_info(out_file, ref_file)
    end

    @timed_testset "reinterpolate with nonuniform data" begin
      # Create and test output with reinterpolation (default options: `reinterpolate=true, data_is_uniform=false`)
      @test_nowarn trixi2vtk(joinpath(outdir, "solution_000005.h5"), output_directory=outdir)
      out_file = joinpath(outdir,"solution_000005.vtu")
      ref_file = get_test_reference_file("dgsem_rotor_amr_reinterp_05.vtu",
                                         joinpath("2d", "p4estmesh", "dgsem_rotor_amr_reinterp_05.vtu"))
      compare_point_info(out_file, ref_file)
    end

    @timed_testset "do not reinterpolate with nonuniform data" begin
      # Create and test output without reinterpolation on LGL nodes
      @test_nowarn trixi2vtk(joinpath(outdir, "solution_000005.h5"), output_directory=outdir, reinterpolate=false)
      out_file = joinpath(outdir,"solution_000005.vtu")
      ref_file = get_test_reference_file("dgsem_rotor_amr_no_reinterp_05.vtu",
                                         joinpath("2d", "p4estmesh", "dgsem_rotor_amr_no_reinterp_05.vtu"))
      compare_point_info(out_file, ref_file)
    end

    @timed_testset "do not reinterpolate with uniform data" begin
      # Create and test output without reinterpolation on uniform nodes
      # OBS! This is a dummy test just to exercise code. The resulting plot will look weird.
      @test_nowarn trixi2vtk(joinpath(outdir, "solution_000005.h5"), output_directory=outdir, reinterpolate=false, data_is_uniform=true)
      out_file = joinpath(outdir,"solution_000005.vtu")
      ref_file = get_test_reference_file("dgsem_rotor_amr_no_reinterp_uniform_05.vtu",
                                         joinpath("2d", "p4estmesh", "dgsem_rotor_amr_no_reinterp_uniform_05.vtu"))
      compare_point_info(out_file, ref_file)
    end
  end
end

# Clean up afterwards: delete Trixi output directory
@test_nowarn rm(outdir, recursive=true)

end
