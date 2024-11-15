module TestManual

using Test
using Trixi2Vtk
using Documenter

include("test_trixi2vtk.jl")

# Start with a clean environment: remove Trixi output directory if it exists
outdir = "out"
isdir(outdir) && rm(outdir, recursive=true)

@testset "Manual" begin
  @timed_testset "Documenter" begin
    DocMeta.setdocmeta!(Trixi2Vtk, :DocTestSetup, :(using Trixi2Vtk); recursive=true)
    doctest(Trixi2Vtk, manual=false)
  end

  @timed_testset "trixi2vtk error triggers" begin
    isdir(outdir) && rm(outdir, recursive=true)
    run_trixi(joinpath(examples_dir(), "tree_2d_dgsem", "elixir_advection_extended.jl"), maxiters=1)

    @testset "no input file" begin
      @test_throws ErrorException trixi2vtk()
    end

    @testset "no such files" begin
      @test_throws ErrorException trixi2vtk("this_does_not_exist")
    end

    @testset "unsupported file format" begin
      @test_throws ErrorException trixi2vtk(joinpath(outdir, "solution_" * LEADING_ZEROS * "000000.h5");
                                            output_directory=outdir,
                                            format=:does_not_exist)
    end
  end

  @timed_testset "trixi2vtk for mesh file" begin
    @test_nowarn trixi2vtk(joinpath(outdir, "mesh.h5"); output_directory=outdir)
  end

  @testset "trixi2vtk set number of output nodes" begin
    @test_nowarn trixi2vtk(joinpath(outdir, "solution_" * LEADING_ZEROS * "000000.h5"); nvisnodes=0)

    @test_nowarn trixi2vtk(joinpath(outdir, "solution_" * LEADING_ZEROS * "000000.h5"); nvisnodes=5)
  end

  @timed_testset "pvd_filenames" begin
    @test Trixi2Vtk.pvd_filenames("", "manual", "out") == (joinpath("out", "manual"), joinpath("out", "manual_celldata"))
    @test_throws ErrorException Trixi2Vtk.pvd_filenames(("a", "b"), nothing, "out")
  end

  @timed_testset "expand_filename_patterns" begin
    mktemp() do path, _
      # Test that absolute filepaths work
      expanded = Trixi2Vtk.expand_filename_patterns([path])
      @test basename(expanded[1]) == basename(path)
    end
  end
end

# Clean up afterwards: delete Trixi output directory
@test_nowarn rm(outdir, recursive=true)

end # module
