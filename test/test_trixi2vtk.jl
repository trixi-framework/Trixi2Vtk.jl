using Test: @test_nowarn, @test
using SHA
using Trixi
using Trixi2Vtk

# pathof(Trixi) returns /path/to/Trixi/src/Trixi.jl, dirname gives the parent directory
const EXAMPLES_DIR = joinpath(pathof(Trixi) |> dirname |> dirname, "examples")


function run_trixi(parameters_file, parameters...)
  @test_nowarn Trixi.run(joinpath(EXAMPLES_DIR, parameters_file); parameters...)
end


function sha1file(filename)
  open(filename) do f
    hash = bytes2hex(sha1(f))
  end

  return hash
end


function test_trixi2vtk_run(filenames, outdir; hashes=nothing, kwargs...)
  @test_nowarn Trixi2Vtk.run(filenames, output_directory=outdir, kwargs...)

  if !isnothing(hashes)
    for filename, hash in hashes
      @test hash == sha1file(joinpath(outdir, filename))
    end
  end
end
