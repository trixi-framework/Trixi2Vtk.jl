using Test

# We use the `TRIXI_TEST` environment variable to determine the subset of tests to execute.
const TRIXI_TEST = get(ENV, "TRIXI_TEST", "all")

@time @testset "Trixi2Vtk" begin
  @time if TRIXI_TEST == "all" || TRIXI_TEST == "2d" || TRIXI_TEST == "upstream"
    include("test_2d.jl")
  end

  @time if TRIXI_TEST == "all" || TRIXI_TEST == "3d" || TRIXI_TEST == "upstream"
    include("test_3d.jl")
  end

  @time if TRIXI_TEST == "all" || TRIXI_TEST == "manual"
    include("test_manual.jl")
  end
end
