using Test: @test_nowarn, @test
using SHA
using Trixi
using Trixi2Vtk

# pathof(Trixi) returns /path/to/Trixi/src/Trixi.jl, dirname gives the parent directory
const EXAMPLES_DIR = joinpath(pathof(Trixi) |> dirname |> dirname, "examples", "2d")


function run_trixi(parameters_file; parameters...)
  @test_nowarn Trixi.run(joinpath(EXAMPLES_DIR, parameters_file); parameters...)
end


function sha1file(filename)
  open(filename) do f
    bytes2hex(sha1(f))
  end
end


function test_trixi2vtk(filenames, outdir; hashes=nothing, kwargs...)
  @test_nowarn trixi2vtk(joinpath(outdir, filenames); output_directory=outdir, kwargs...)

  if !isnothing(hashes)
    for (filename, hash_expected) in hashes
      hash_measured = sha1file(joinpath(outdir, filename))
      @test hash_expected == hash_measured
    end
  end
end
