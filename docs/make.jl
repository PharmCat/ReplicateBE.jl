using Documenter, ReplicateBE

makedocs(modules=[ReplicateBE],
    sitename="ReplicateBE",
    authors = "Vladimir Arnautov",
    linkcheck = false,
    doctest = false,
    pages = [
        "Home" => "index.md",
        "Examples" => "examples.md",
        "Syntax" => "syntax.md",
        "Details" => "details.md",
        "Validation" => "testval.md",
        "Structures" => "struct.md",
        "Api" => "api.md"
    ])

#https://github.com/PharmCat/ReplicateBE.jl/blob/master/docs/
#https://github.com/PharmCat/ReplicateBE.jl/blob/master/docs/
