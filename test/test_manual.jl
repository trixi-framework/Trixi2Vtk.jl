module TestManual

using Test
using Trixi2Vtk
using Documenter

include("test_trixi2vtk.jl")

# Start with a clean environment: remove Trixi output directory if it exists
outdir = "out"
isdir(outdir) && rm(outdir, recursive=true)

@testset "Manual" begin
  @testset "Documenter" begin
    DocMeta.setdocmeta!(Trixi2Vtk, :DocTestSetup, :(using Trixi2Vtk); recursive=true)
    doctest(Trixi2Vtk, manual=false)
  end

  @testset "trixi2vtk error triggers" begin
    isdir(outdir) && rm(outdir, recursive=true)
    run_trixi(joinpath("2d", "elixir_advection_extended.jl"), maxiters=1)

    @testset "no input file" begin
      @test_throws ErrorException trixi2vtk()
    end

    @testset "no such files" begin
      @test_throws ErrorException trixi2vtk("this_does_not_exist")
    end

    @testset "unsupported file format" begin
      @test_throws ErrorException trixi2vtk(joinpath(outdir, "solution_000000.h5");
                                            output_directory=outdir,
                                            format=:does_not_exist)
    end
  end

  @testset "trixi2vtk for mesh file" begin
    @test_nowarn trixi2vtk(joinpath(outdir, "mesh.h5"); output_directory=outdir)
  end

  @testset "pvd_filenames" begin
    @test Trixi2Vtk.pvd_filenames("", "manual", "out") == (joinpath("out", "manual"), joinpath("out", "manual_celldata"))
    @test_throws ErrorException Trixi2Vtk.pvd_filenames(("a", "b"), nothing, "out")
  end
end

# Clean up afterwards: delete Trixi output directory
@test_nowarn rm(outdir, recursive=true)

end # module
