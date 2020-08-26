using Documenter
using Trixi2Vtk

# Define module-wide setup such that the respective modules are available in doctests
DocMeta.setdocmeta!(Trixi2Vtk,
                    :DocTestSetup,
                    :(push!(LOAD_PATH, ".."); using Trixi2Vtk);
                    recursive=true)

# Make documentation
makedocs(
    # Specify modules for which docstrings should be shown
    modules = [Trixi2Vtk],
    # Set sitename
    sitename="Trixi2Vtk",
    # Provide additional formatting options
    format = Documenter.HTML(
        # Disable pretty URLs during manual testing
        prettyurls = get(ENV, "CI", nothing) == "true",
        # Explicitly add favicon as asset
        assets = ["assets/favicon.ico"],
        # Set canonical URL to GitHub pages URL
        canonical = "https://trixi-framework.github.io/Trixi2Vtk.jl/stable"
    ),
    # Explicitly specify documentation structure
    pages = [
        "Home" => "index.md",
        "Reference" => [
            "Trixi2Vtk" => "reference/trixi2vtk.md",
        ],
        "License" => "license.md"
    ]
)

deploydocs(
    repo = "github.com/trixi-framework/Trixi2Vtk.jl",
)
