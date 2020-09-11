struct RBEDataStructure{T <: AbstractFloat}
    factors::Vector
    Xv::Vector{SubArray{T}}
    Zv::Vector{SubArray{T}}
    yv::Vector{SubArray{T}}
    p::Int
    N::Int
    n::Int
    remlc::T
    maxobs::Int
end

function Base.show(io::IO, rbeds::RBEDataStructure)
    print(io,"""
    .factors: $(rbeds.factors)
    .Xv: length = $(length(rbeds.Xv))
    .Zv: length = $(length(rbeds.Zv))
    .yv: length = $(length(rbeds.yv))
    .yv: length = $(length(rbeds.yv))
    .p = $(rbeds.p)
    .N = $(rbeds.N)
    .n = $(rbeds.n)
    .remlc = $(rbeds.remlc)
    .maxobs = $(rbeds.maxobs)
<<<<<<< HEAD
=======
    .mem: MemCache
>>>>>>> 891996a132783ac52ac2996a5b5092f166adff31
    """)
end

struct RBEResults{T <: AbstractFloat}
    reml::T                         # logREML
    β::Vector{T}                    # β Vector
    theta::Vector{T}                # Variance parameter vector
    fixed::EffectTable
    G::AbstractMatrix{T}            # G matrix
    H::Matrix{T}                    # Hessian matrix
    A::Matrix{T}                    # asymptotic variance-covariance matrix ofb θ
    C::Matrix{T}                    # C var(β) p×p variance-covariance matrix
    gradc::Vector{Matrix}
    preoptim::Union{Optim.MultivariateOptimizationResults, Nothing}             # Pre-optimization result object
    optim::Optim.MultivariateOptimizationResults                                # Optimization result object
end

function Base.show(io::IO, rbeds::RBEResults)
    print(io,"""
    .reml
    .β
    .theta
    .fixed
    .G
    .H
    .A
    .C
    .gradc
    .preoptim
    .optim
    """)
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
