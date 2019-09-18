#Model Frame utils
function lvec(mm::ModelMatrix, f::Int)
    l = zeros(length(mm.assign))
    for i = 1:length(l)
        if mm.assign == f l[i] = 1 end
    end
end

function lvec(MF::ModelFrame, f::Union{Symbol, AbstractTerm})::Int
    l = length(MF.f.rhs.terms)
    vec = Array{Int, 1}(undef, 0)
    for i = 1:l
        if isa(MF.f.rhs.terms[i], InterceptTerm)
            if f == InterceptTerm
                vec = vcat(vec, ones(length(1)))
            end
        elseif MF.f.rhs.terms[i].sym == f
            vec = vcat(vec, ones(length(MF.f.rhs.terms[i].contrasts.termnames)))
        end
    end
    return vec
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
#Return length by Symbol
function termmodellen(MF::ModelFrame, symbol::Symbol)::Int
    id = findterm(MF, symbol)
    return length(MF.f.rhs.terms[id].contrasts.termnames)
end
#Confidence interval
function StatsBase.confint(obj::RBE, alpha::Float64; expci::Bool = false, inv::Bool = false, df = :sat)
    ifv = 1
    if inv ifv = -1 end
    if isa(df, Array{Float64, 1})
        if length(obj.fixed.df) != length(df)
            df = obj.fixed.df
        end
    elseif isa(df, Symbol)
        if df == :df2
            df  = zeros(length(obj.fixed.df))
            df .= obj.df2
        else
            df = obj.fixed.df
        end
    end

    a = Array{Tuple{Float64, Float64},1}(undef, length(obj.fixed.est)-1)
    for i = 2:length(obj.fixed.est)
        a[i-1] = calcci(obj.fixed.est[i]*ifv, obj.fixed.se[i], df[i], alpha, expci)
    end
    return Tuple(a)
end
function calcci(x::Float64, se::Float64, df::Float64, alpha::Float64, expci::Bool)::Tuple{Float64, Float64}
    q = quantile(TDist(df), 1.0-alpha/2)
    if !expci
        return x-q*se, x+q*se
    else
        return exp(x-q*se), exp(x+q*se)
    end
end
function Base.show(io::IO, obj::Tuple{Vararg{Tuple{Float64, Float64}}})
    for i in obj
        println(io, i)
    end
end
function reml2(obj::RBE, θ::Array{Float64, 1})
    return -2*reml(obj.yv, obj.Zv, rank(ModelMatrix(obj.model).m), obj.Xv, θ, obj.fixed.est)
end
function contrast(obj::RBE, L::Matrix{T}) where T <: Real
    lcl  = L*obj.C*L'
    lclr = rank(lcl)
    return (L*obj.fixed.est)'*inv(lcl)*(L*obj.fixed.est)/lclr
end
function lsm(obj::RBE, L::Matrix{T}) where T <: Real
    lcl  = L*obj.C*L'
    return L*obj.fixed.est, sqrt.(lcl)
end
function emm(obj::RBE, fm, lm)
    La = lmean(obj::RBE)'
    L  = La .* fm'
    L  = L  .+ lm'
    return lsm(obj, Matrix(L'))
end
function lmean(obj::RBE)
    #coef  = Array{Float64, 1}(undef, length(obj.factors))
    L    = zeros(length(obj.fixed.est))
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
    return Matrix(L')
end
function fixedeffect()
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
