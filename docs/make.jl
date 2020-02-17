using Documenter, Juniper

makedocs(
    # See https://github.com/JuliaDocs/Documenter.jl/issues/868
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true", mathengine = Documenter.MathJax()),
    # See https://github.com/JuliaOpt/JuMP.jl/issues/1576
    strict = true,
    sitename = "Juniper",
    pages = [
        "Home" => "index.md",
        "Options" => "options.md",
        "Extras" => "extras.md",
#        "Developer" => [],
#        "Library" => "library.md"
    ]
)

deploydocs(
    repo = "github.com/lanl-ansi/Juniper.jl.git",
)