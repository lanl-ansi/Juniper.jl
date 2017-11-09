using Documenter, MINLPBnB

makedocs(
    modules = [MINLPBnB],
    format = :html,
    sitename = "MINLPBnB",
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
    repo = "github.com/Wikunia/MINLPBnB.git",
    julia = "0.6"
)