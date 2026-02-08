using Documenter
using CellBase

makedocs(
    sitename="CellBase.jl",
    format=Documenter.HTML(),
    modules=[CellBase],
    pages=[
        "Home" => "index.md",
        "Getting Started" => "getting_started.md",
        "Guides" => [
            "Working with Cells" => "guides/working_with_cells.md",
            "File I/O" => "guides/file_io.md",
            "Building Structures" => "guides/building_structures.md",
            "Spacegroup Operations" => "guides/spacegroup_operations.md",
        ],
    ],
    warnonly=true,
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo="github.com/zhubonan/CellBase.git",
    devbranch="master",
    devurl="master",
    target="build",
    branch="gh-pages",
)
