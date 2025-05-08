using Test: @test_nowarn, @test, @testset
using Downloads: download
using Trixi
using Trixi2Vtk
using ReadVTK


# Leading zeros will always be "000" since we have restricted Trixi.jl to 0.9+
const LEADING_ZEROS = "000"


function run_trixi(elixir; kwargs...)
  # evaluate examples in the scope of the module they're called from
  trixi_include(@__MODULE__, elixir; kwargs...)
end


"""
    get_reference_file(filename, remotename; head="main", output_directory=".", force=false)

Retrieve a reference file from the
[`Trixi2Vtk_reference_files` repository](https://github.com/trixi-framework/Trixi2Vtk_reference_files)
at commit/branch `head` and store it in the `output_directory`. If the file already
exists locally, do not download the file again unless `force` is true.
Provide the `remotename` because the file structure in the `Trixi2Vtk_reference_files` repository
is different than what will be created locally. Return the local path to the downloaded file.
"""
function get_reference_file(filename, remotename; head="main", output_directory=".", force=false)
  filepath = joinpath(output_directory, filename)
  if !isfile(filepath) || force
    url = ("https://github.com/trixi-framework/Trixi2Vtk_reference_files/raw/"
           * head
           * "/reference_files/"
           * remotename)

    download(url, filepath)
  end

  return filepath
end

# Commit in the reference file repository for which the test files will be downloaded
# Note: The purpose of using a specific commit hash (instead of `main`) is to be able to tie a given
#       version of Trixi2Vtk to a specific version of the test file repository. This way, also tests
#       for older Trixi2Vtk releases should continue to work.
const TEST_REFERENCE_COMMIT = "857ed155d1b3aeff20474b284e83bab01748237a"

# Local folder to store downloaded reference files. If you change this, also adapt `../.gitignore`!
const TEST_REFERENCE_DIR = "reference_files"


get_test_reference_file(filename, remotename) = get_reference_file(filename, remotename,
                                                                   head=TEST_REFERENCE_COMMIT,
                                                                   output_directory=TEST_REFERENCE_DIR)


# Start with a clean environment: remove reference file directory if it exists
isdir(TEST_REFERENCE_DIR) && rm(TEST_REFERENCE_DIR, recursive=true)
mkpath(TEST_REFERENCE_DIR)


"""
     compare_cell_data(out_filename, ref_filename; atol=500*eps(), rtol=sqrt(eps()))

Test values from the VTK file header and actual (possibly interpolated) cell data. Uses
`out_filename` created during testing and compares against `ref_filename` that comes
from the
[`Trixi2Vtk_reference_files` repository](https://github.com/trixi-framework/Trixi2Vtk_reference_files).
"""
function compare_cell_data(out_filename, ref_filename; atol=500*eps(), rtol=sqrt(eps()))
  ref_vtk = VTKFile(ref_filename)
  vtk = VTKFile(out_filename)

  # Compare header information
  @test vtk.byte_order == ref_vtk.byte_order
  @test vtk.compressor == ref_vtk.compressor
  @test vtk.file_type == ref_vtk.file_type
  @test vtk.version == ref_vtk.version

  # check that the number of cells and points match
  @test vtk.n_cells == ref_vtk.n_cells
  @test vtk.n_points == ref_vtk.n_points

  # Compare data information

  # extract the data sets
  out_cell_data = get_cell_data(vtk)
  ref_cell_data = get_cell_data(ref_vtk)

  # check that the number of variables and their names match
  @test length(out_cell_data.names) == length(ref_cell_data.names)
  @test out_cell_data.names == ref_cell_data.names

  # Note!!
  # Occasionally, for the last equation variable there is an issue
  # that the size of the data extracted with
  # [`ReadVTK.jl`](https://github.com/trixi-framework/ReadVTK.jl)
  # is one byte larger than expected. Somehow this is related to
  # the zlib format and/or how
  # [`WriteVTK.jl`](https://github.com/jipolanco/WriteVTK.jl) created
  # the file in the first place. Even though there is this access error,
  # the `vtu` file is valid and can be used with ParaView or VisIt without issue.
  # Therefore, we check all the variable data except for the last one,
  # unless the file only has a single variable stored as cell data.
  eqn_check_number = length(ref_cell_data.names) == 1 ? length(ref_cell_data.names) : length(ref_cell_data.names)-1

  # check that the actual plot data is (approximately) the same
  for i in 1:eqn_check_number
    variable_name = ref_cell_data.names[i]
    out_data = get_data(out_cell_data[variable_name])
    ref_data = get_data(ref_cell_data[variable_name])
    @test isapprox(out_data, ref_data, atol=atol, rtol=rtol)
  end
end


"""
    compare_point_data(out_filename, ref_filename; atol=500*eps(), rtol=sqrt(eps()))

Test values from the VTK file header and actual (possibly interpolated) point data. Uses
`out_filename` created during testing and compares against `ref_filename` that comes
from the
[`Trixi2Vtk_reference_files` repository](https://github.com/trixi-framework/Trixi2Vtk_reference_files).
"""
function compare_point_data(out_filename, ref_filename; atol=500*eps(), rtol=sqrt(eps()))
  # Load the data from both files
  ref_vtk = VTKFile(ref_filename)
  vtk = VTKFile(out_filename)

  # Compare header information
  @test vtk.byte_order == ref_vtk.byte_order
  @test vtk.compressor == ref_vtk.compressor
  @test vtk.file_type == ref_vtk.file_type
  @test vtk.version == ref_vtk.version

  # check that the number of cells and points match
  @test vtk.n_cells == ref_vtk.n_cells
  @test vtk.n_points == ref_vtk.n_points

  # Compare data information

  # extract the data sets
  out_point_data = get_point_data(vtk)
  ref_point_data = get_point_data(ref_vtk)

  # check that the number of variables and their names match
  @test length(out_point_data.names) == length(ref_point_data.names)
  @test out_point_data.names == ref_point_data.names

  # Note!!
  # Occasionally, for the last equation variable there is an issue
  # that the size of the data extracted with
  # [`ReadVTK.jl`](https://github.com/trixi-framework/ReadVTK.jl)
  # is one byte larger than expected. Somehow this is related to
  # the zlib format and/or how
  # [`WriteVTK.jl`](https://github.com/jipolanco/WriteVTK.jl) created
  # the file in the first place. Even though there is this access error,
  # the `vtu` file is valid and can be used with ParaView or VisIt without issue.
  # Therefore, we check all the variable data except for the last one.

  eqn_check_number = length(ref_point_data.names) == 1 ? length(ref_point_data.names) : length(ref_point_data.names)-1

  # check that the actual plot data is (approximately) the same
  for i in 1:eqn_check_number
    variable_name = ref_point_data.names[i]
    out_data = get_data(out_point_data[variable_name])
    ref_data = get_data(ref_point_data[variable_name])
    @test isapprox(out_data, ref_data, atol=atol, rtol=rtol)
  end
end


"""
    @timed_testset "name of the testset" #= code to test #=

Similar to `@testset`, but wraps the execution of the testset using `@time`
and prints the name of the testset after execution.
"""
macro timed_testset(name, expr)
  @assert name isa String
  quote
    @time @testset $name $expr
    flush(stdout)
    @info("Testset " * $name * " finished.\n")
    flush(stdout)
  end
end
