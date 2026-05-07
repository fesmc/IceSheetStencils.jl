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
        "SSA" => [
            "Theory" => "ssa/theory.md",
            "Discretization" => [
                "Overview"                        => "ssa/discretization/index.md",
                "C-grid"                          => "ssa/discretization/c_grid.md",
                "Energy-functional formulation"   => "ssa/discretization/energy_functional.md",
                "Extending to a realistic solver" => "ssa/discretization/extending.md",
            ],
        ],
        "Mass conservation" => [
            "Theory" => "mass_conservation/theory.md",
            "Discretization" => [
                "2-D level-set method" => "mass_conservation/discretization/levelset.md",
            ],
        ],
        "API reference" => "api.md",
    ],
)

deploydocs(
    repo         = "github.com/fesmc/IceSheetStencils.jl.git",
    devbranch    = "main",
    push_preview = true,
)
