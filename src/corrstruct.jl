abstract type CorrStruct end

struct ScaledIdentity <: CorrStruct
    p::Union{Int, Vector{Int}, UnitRange{Int}}
    function ScaledIdentity(p::Union{Int, Vector{Int}, UnitRange{Int}})::ScaledIdentity
        new(p)::ScaledIdentity
    end
    function ScaledIdentity(p::NamedTuple)::ScaledIdentity
        new(p.p)::ScaledIdentity
    end
end
function scaledIdentity(p)::ScaledIdentity
    ScaledIdentity(p)
end


struct VarianceComponents <: CorrStruct #
    p::Union{Int, Vector{Int}, UnitRange{Int}}
    function VarianceComponents(p::Union{Int, Vector{Int}, UnitRange{Int}})::VarianceComponents
        new(p)::VarianceComponents
    end
    function VarianceComponents(p::NamedTuple)::VarianceComponents
        new(p.p)::VarianceComponents
    end
end
function varianceComponents(p)::VarianceComponents
    VarianceComponents(p)
end


struct CompoundSymmetry <: CorrStruct
    p::Union{Int, Vector{Int}, UnitRange{Int}}
    function CompoundSymmetry(p::Union{Int, Vector{Int}, UnitRange{Int}})::CompoundSymmetry
        new(p)::CompoundSymmetry
    end
    function CompoundSymmetry(p::NamedTuple)::CompoundSymmetry
        new(p.p)::CompoundSymmetry
    end
end
function compoundSymmetry(p)::CompoundSymmetry
    CompoundSymmetry(p)
end

struct HeterogeneousCompoundSymmetry <: CorrStruct
    p::Union{Int, Vector{Int}, UnitRange{Int}}
    rho::Int
    function HeterogeneousCompoundSymmetry(p::Union{Int, Vector{Int}, UnitRange{Int}}, rho::Int)::HeterogeneousCompoundSymmetry
        new(p, rho)::HeterogeneousCompoundSymmetry
    end
    function HeterogeneousCompoundSymmetry(p::NamedTuple)::HeterogeneousCompoundSymmetry
        new(p.p, p.rho)::HeterogeneousCompoundSymmetry
    end
end
function heterogeneousCompoundSymmetry(p, rho)::HeterogeneousCompoundSymmetry
    HeterogeneousCompoundSymmetry(p, rho)
end
function heterogeneousCompoundSymmetry(p)::HeterogeneousCompoundSymmetry
    HeterogeneousCompoundSymmetry(p)
end


struct VarianceStruct{T1 <: CorrStruct, T2 <: CorrStruct}
    g::T1
    r::T2
    function VarianceStruct(g, r)::VarianceStruct
        new{typeof(g), typeof(r)}(g, r)::VarianceStruct
    end
end


function getrange(vt::typeof(scaledIdentity), zr::Int)
    (1,)
end
function getrange(vt::typeof(varianceComponents), zr::Int)
    (zr,)
end
function getrange(vt::typeof(compoundSymmetry), zr::Int)
    (2,)
end
function getrange(vt::typeof(heterogeneousCompoundSymmetry), zr::Int)
    (zr, 1)
end

function makevariancestruct(g, r, zr)
    #println(g)
    #println(r)
    glt = getrange(g, zr)
    rlt = getrange(r, zr)
    i    = 1
    rp   = nothing
    rrho = nothing
    gp   = nothing
    grho = nothing

    if rlt[1] == 1
        rp = 1
    else
        rp = i:rlt[1]
    end
    i  += rlt[1]

    if length(rlt) == 2
        if rlt[2] == 1 rrho = i else rrho = i:(i+rlt[2]-1) end
        i  += rlt[2]
    end

    if glt[1] == 1
        gp = i
    else
        gp = i:(i+glt[1]-1)
    end
    i  += glt[1]

    if length(glt) == 2
        if glt[2] == 1 grho = i else grho = i:(i+glt[2]-1) end
    end

    gv = g((;:p => gp, :rho => grho))
    rv = r((;:p => rp, :rho => rrho))
    return VarianceStruct(gv, rv)
end

function gmat(theta, vs::VarianceStruct{HeterogeneousCompoundSymmetry, T}) where T <: CorrStruct
    #g = Matrix(Diagonal(theta[vs.g.p]))
    g = sparse(Diagonal(theta[vs.g.p]))
    for r = 1:size(g, 1)
        for c = size(g, 1):-1:r
            if r == c g[r,c] = g[c,c]*g[r,r] else g[r,c] = g[c,c]*g[r,r]*theta[vs.g.rho] end
        end
    end
    Symmetric(g)
end

function rmat(theta, Z, vs::VarianceStruct{T, VarianceComponents}) where T <: CorrStruct
    #=
    r = Diagonal(zeros(eltype(theta), size(Z, 1)))
    for i = 1:size(Z, 1)
        for c = 1:size(Z,2)
            r[i,i] += theta[vs.r.p][c] * Z[i, c]
        end
        r[i,i] = r[i,i] ^ 2
    end
    r=#
    Diagonal(Z*(theta[vs.r.p].^2))
end


function Base.show(io::IO, vs::VarianceStruct)
    println(io, "Random:")
    println(io, "   $(vs.g)")
    println(io, "Repeated:")
      print(io, "   $(vs.r)")
end
function Base.show(io::IO, vs::ScaledIdentity)
    print(io, "ScaledIdentity: p = $(length(vs.p)) ($(vs.p))")
end
function Base.show(io::IO, vs::VarianceComponents)
    print(io, "VarianceComponents: p = $(length(vs.p)) ($(vs.p))")
end
function Base.show(io::IO, vs::ScaledIdentity)
    print(io, "CompoundSymmetry: p = $(length(vs.p)) ($(vs.p))")
end
function Base.show(io::IO, vs::HeterogeneousCompoundSymmetry)
    print(io, "HeterogeneousCompoundSymmetry: p = $(length(vs.p)) ($(vs.p)); rho = $(length(vs.rho)) ($(vs.rho))")
end
