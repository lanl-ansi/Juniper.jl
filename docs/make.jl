using Documenter, Juniper

makedocs(
    modules = [Juniper],
    format = :html,
    sitename = "Juniper",
    pages = [
        "Home" => "index.md",
        "Options" => "options.md",
        "Extras" => "extras.md",
#        "Developer" => [],
        # "Library" => "library.md"
    ]
)

deploydocs(
    deps = nothing,
    make = nothing,
    target = "build",
    repo = "github.com/lanl-ansi/Juniper.jl.git",
    julia = "0.6"
)