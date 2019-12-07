struct RBEDataStructure
    factors::Vector
    Xv::Vector
    Zv::Vector
    yv::Vector
    p::Int
    N::Int
    n::Int
    remlc::AbstractFloat
    maxobs::Int
    mem::MemCache
end

function rebuildcache(data, type)
        data.mem.svec[1] = zeros(type, data.maxobs)
end

struct RBEResults{T <: AbstractFloat}
    reml::T
    β::Vector{T}
    theta::Vector{T}
    H::Matrix{T}
    A::Matrix{T}
    C::Matrix{T}
    gradc::Vector{Matrix}
end

#Design structure
"""
```julia
struct Design
    obs::Int              # Total observation
    subj::Int             # Subjects
    sqn::Int              # Sequences
    pn::Int               # Periods
    fn::Int               # Formulations
    sbf::Vector{Int}      # Subject by formulation
    rankx::Int            # Rank X
    rankxz::Int           # Rank XZ
    df2::Int              # DF (robust): Subjects - Sequence
    df3::Int              # DF (contain): Observations - Rank(XZ)
    df4::Int              # DF (Subj(Form)): sum(sbf) - sqn - 1                 # Experimental
end
```

Trial design description.
"""
struct Design
    obs::Int
    subj::Int
    sqn::Int
    pn::Int
    fn::Int
    sbf::Vector{Int}
    rankx::Int
    rankxz::Int
    df2::Int
    df3::Int
    df4::Int
    function Design(obs, subj, sqn, pn, fn, sbf, rankx, rankxz)
        new(obs, subj, sqn, pn, fn, sbf, rankx, rankxz, subj - sqn, obs - rankxz, sum(sbf) - sqn - 1)::Design
    end
end

function Base.show(io::IO, d::Design)
    println(io, "Trial design")
    println(io, "  Observations:          $(d.obs)")
    println(io, "  Subjects:              $(d.subj)")
    println(io, "  Sequences:             $(d.sqn)")
    println(io, "  Periods:               $(d.pn)")
    println(io, "  Formulations:          $(d.fn)")
    print(io,   "  Subjects(Formulation): $(d.sbf[1])")
    for i = 2:length(d.sbf)
        print(io, "/$(d.sbf[i])")
    end
    println(io, "")
    println(io, "  Rank X:                $(d.rankx)")
    println(io, "  Rank XZ:               $(d.rankxz)")
    println(io, "  DF (robust):           $(d.df2)")
    println(io, "  DF (contain):          $(d.df3)")
    print(io,   "  DF (Subj(Form)):       $(d.df4)")

end
