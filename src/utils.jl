#Model Frame utils
function lvec(mm::ModelMatrix, f::Int)
    l = zeros(length(mm.assign))
    for i = 1:length(l)
        if mm.assign == f l[i] = 1 end
    end
end
# L Matrix for TYPE III
function lmatrix(mf::ModelFrame, f::Union{Symbol, AbstractTerm})
    l   = length(mf.f.rhs.terms)
    id  = findterm(mf, f)
    n   = length(mf.f.rhs.terms[id].contrasts.termnames)
    lm  = zeros(n, length(coefnames(mf)))
    vec = Array{Int, 1}(undef, 0)
    for i = 1:l
        if isa(mf.f.rhs.terms[i], InterceptTerm)
            if f == InterceptTerm
                vec = vcat(vec, ones(1))
            else
                vec = vcat(vec, zeros(1))
            end
        elseif mf.f.rhs.terms[i].sym == f
            vec = vcat(vec, ones(length(mf.f.rhs.terms[i].contrasts.termnames)))
        else
            vec = vcat(vec, zeros(length(mf.f.rhs.terms[i].contrasts.termnames)))
        end
    end
    r     = 1
    for i = 1:size(lm, 2)
        if vec[i] == 1
            lm[r, i] = 1
            r +=1
        end
    end
    return lm
end
#Find by Symbol
function findterm(MF::ModelFrame, f::Union{Symbol, AbstractTerm})::Int
    l = length(MF.f.rhs.terms)
    for i = 1:l
        if isa(MF.f.rhs.terms[i], InterceptTerm)
            if f == InterceptTerm return i else continue end
        end
        if MF.f.rhs.terms[i].sym == f return i end
    end
    return 0
end
#Return term fixed effect length by symbol
function termmodellen(MF::ModelFrame, symbol::Symbol)::Int
    id = findterm(MF, symbol)
    return length(MF.f.rhs.terms[id].contrasts.termnames)
end
#Return term levels count by symbol
function termmodelleveln(MF::ModelFrame, symbol::Symbol)::Int
    id = findterm(MF, symbol)
    return length(MF.f.rhs.terms[id].contrasts.levels)
end

function calcci(x::Float64, se::Float64, df::Float64, alpha::Float64, expci::Bool)::Tuple{Float64, Float64}
    q = quantile(TDist(df), 1.0-alpha/2)
    if !expci
        return x-q*se, x+q*se
    else
        return exp(x-q*se), exp(x+q*se)
    end
end
#return contrast F
function contrast(obj::RBE, L::Matrix{T}) where T <: Real
    lcl  = L*obj.C*L'
    lclr = rank(lcl)
    return obj.fixed.est'*L'*inv(lcl)*L*obj.fixed.est/lclr  #?
end
#
function lsm(obj::RBE, L::Matrix{T}) where T <: Real
    lcl  = L*obj.C*L'
    return L*obj.fixed.est, sqrt.(lcl)
end

function emm(obj::RBE, fm::Matrix, lm::Matrix)
    La = lmean(obj::RBE)
    L  = La .* fm
    L  = L  .+ lm
    return lsm(obj, Matrix(L))
end
#General mean contrast L matrix 1xp
function lmean(obj::RBE)
    L    = zeros(1, length(obj.fixed.est))
    L[1] = 1.0
    it    = 2
    for f in obj.factors
        term = findterm(obj.model, f)
        len  = length(obj.model.f.rhs.terms[term].contrasts.termnames)
        dev  = 1/length(obj.model.f.rhs.terms[term].contrasts.levels)
        for i = 1:len
            L[it] = dev
            it  += 1
        end
    end
    return L
end
#-------------------------------------------------------------------------------
function checkdata(X, Z, Xv, Zv, y)
    if size(Z)[2] != 2 error("Size random effect matrix != 2. Not implemented yet!") end
    if length(Xv) != length(Zv) error("Length Xv != Zv !!!") end
    for i = 1:length(Xv)
        if size(Xv[i])[1]  != size(Zv[i])[1] error("Row num of subject $i Xv != Zv !!!") end
    end
end
#-------------------------------------------------------------------------------
#show EffectTable
function addspace(s::String, n::Int)::String
    for i = 1:n
        s = s * Char(' ')
    end
    return s
end
#-------------------------------------------------------------------------------
function geocv(var)
    return sqrt(exp(var)-1.0)
end
#-------------------------------------------------------------------------------
#Subject number by factor (formulation)
function sbjnbyf(df, subj, fac, f)
    sbj = Array{Any, 1}(undef, 0)
    for i = 1:size(df, 1)
        if df[i, fac] == f
            if !any(x -> x == df[i, subj], sbj) push!(sbj, df[i, subj]) end
        end
    end
    return length(sbj)
end
#-------------------------------------------------------------------------------
function Base.show(io::IO, obj::Tuple{Vararg{Tuple{Float64, Float64}}})
    for i in obj
        println(io, i)
    end
end
