using Documenter, ReplicateBE

makedocs(sitename="ReplicateBE",
    authors = "Vladimir Arnautov",
    linkcheck = false,
    doctest = false)

deploydocs(
    repo    = "github.com/PharmCat/ReplicateBE.jl.git",
    target  = "build",
    deps    = nothing,
    make    = nothing,
)
