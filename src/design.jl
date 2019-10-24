#Design structure
"""
Describes of trial design.
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
