using Documenter, ReplicateBE

makedocs(sitename="ReplicateBE",
    authors = "Vladimir Arnautov",
    linkcheck = false,
    doctest = false,
    pages = [
        "Home" => "index.md",
        "Examples" => "examples.md",
        "Syntax" => "syntax.md",
        "Details" => "details.md",
        "Api" => "api.md"
    ])
