using Documenter
using IceSheetStencils

DocMeta.setdocmeta!(IceSheetStencils, :DocTestSetup,
                    :(using IceSheetStencils); recursive = true)

makedocs(
    sitename = "IceSheetStencils.jl",
    modules  = [IceSheetStencils],
    authors  = "Alexander Robinson",
    format   = Documenter.HTML(
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical  = "https://fesmc.github.io/IceSheetStencils.jl",
        edit_link  = "main",
        repolink   = "https://github.com/fesmc/IceSheetStencils.jl",
    ),
    pages = [
        "Home"   => "index.md",
        "Theory" => "theory.md",
        "Discretization" => [
            "Overview"      => "discretization/index.md",
            "C-grid (SSA)"  => "discretization/c_grid.md",
        ],
        "API reference" => "api.md",
    ],
)

deploydocs(
    repo         = "github.com/fesmc/IceSheetStencils.jl.git",
    devbranch    = "main",
    push_preview = true,
)
