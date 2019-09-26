#Design structure
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
    print(io,   "  DF (Subj(Form)):       $(sum(d.sbf) - length(d.sbf)*d.sqn)")

end
