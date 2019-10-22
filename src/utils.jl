#return contrast table
"""
    contrast(rbe::RBE, L::Matrix; numdf = 1, name = "Contrast", memopt = true)::ContrastTable

Return contrast table for L matrix.

``F = \\frac{\\beta'L'(LCL')^{-1}L\\beta}{rank(LCL')}``

DF for one-dimetion case:

``df = \\frac{2(LCL')^{2}}{g'Ag}``

where ``A = 2H``

where ``g = \\triangledown _{\\theta}(LC^{-1}L')``

DF for multi-dimention case:


"""
function contrast(rbe::RBE, L::Matrix; numdf = 0, name = "Contrast", memopt = true)::ContrastTable
    β       = coef(rbe)
    lcl     = L*rbe.C*L'
    lclr    = rank(lcl)
    F       = β'*L'*inv(lcl)*L*β/lclr
    θ       = theta(rbe)

    if numdf == 0 numdf = rank(L) end

    if rank(L) ≥ 2
        vm      = Array{Float64, 1}(undef, size(L, 1))
        for i = 1:length(vm)
            g       = ForwardDiff.gradient(x -> lclgf(L[i:i,:], L[i:i,:]', rbe.Xv, rbe.Zv, x; memopt = memopt), θ)
            df      = 2*((L[i:i,:]*rbe.C*L[i:i,:]')[1])^2/(g'*(rbe.A)*g)
            vm[i]   = df/(df-2)
        end
        df = 2*sum(vm)/(sum(vm)-rank(L))
    else
        g       = ForwardDiff.gradient(x -> lclgf(L, L', rbe.Xv, rbe.Zv, x; memopt = memopt), θ)
        df      = 2*((lcl)[1])^2/(g'*(rbe.A)*g)
    end
    pval    = ccdf(FDist(numdf, df), F)
    return ContrastTable([name], [F], [numdf], [df], [pval])
end
"""
    estimate(rbe::RBE, L::Matrix; df = :sat, name = "Estimate", memopt = true, alpha = 0.05)

Return estimate table for L 1xp matrix.

``estimate = L\\beta``

``se = \\sqrt{LCL'}``

``t = estimate/se``

For ```df = :sat```:

``df = \\frac{2(LCL')^{2}}{g'Ag}``

where ``A = 2H``

where ``g = \\triangledown _{\\theta}(LC^{-1}L')``

For ```df = :cont``` (contain):

``df = N - rank(ZX)``

CI estimate is ``CI = stimate ± t(alpha, df)*se ``
"""
function estimate(rbe::RBE, L::Matrix; df = :sat, name = "Estimate", memopt = true, alpha = 0.05)
    lcl     = L*rbe.C*L'
    β       = coef(rbe)
    est     = (L*β)[1]
    lclr    = rank(lcl)
    se      = sqrt((lcl)[1])
    #F       = β'*L'*inv(lcl)*L*β/lclr
    if df == :sat
        θ       = theta(rbe)
        g       = ForwardDiff.gradient(x -> lclgf(L, L', rbe.Xv, rbe.Zv, x; memopt = memopt), θ)
        df      = 2*((lcl)[1])^2/(g'*(rbe.A)*g)
    elseif df == :cont
        df      = obj.design.df3
    else
        throw(ArgumentError("df unknown!"))
    end
    t       = ((est)/se)
    pval    = ccdf(TDist(df), abs(t))*2
    return EstimateTable([name], [est], [se], [df], [t], [pval], [est - se*quantile(TDist(df), 1-alpha/2)], [est + se*quantile(TDist(df), 1-alpha/2)], alpha)
end

# L Matrix for TYPE III
"""
    lmatrix(mf::ModelFrame, f::Union{Symbol, AbstractTerm})

L matrix for factor.
"""
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
#
"""
    lsm(rbe::RBE, L::Matrix)

Deprecated.
"""
#Deprecated
function lsm(rbe::RBE, L::Matrix)
    lcl  = L*rbe.C*L'
    return L*coef(rbe), sqrt.(lcl)
end
#
"""
    emm(obj::RBE, fm::Matrix, lm::Matrix)

Matrix mask.
"""
function emm(obj::RBE, fm::Matrix, lm::Matrix)
    La = lmean(obj::RBE)
    L  = La .* fm
    L  = L  .+ lm
    return lsm(obj, Matrix(L))
end
#General mean contrast L matrix 1xp
"""
    lmean(obj::RBE)

Return L-matrix for general mean.
"""
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
function addspace(s::String, n::Int; first = false)::String
    if n > 0
        for i = 1:n
            if first s = Char(' ') * s else s = s * Char(' ') end
        end
    end
    return s
end

function printmatrix(io::IO, m::Matrix)
    sm = string.(m)
    lv = maximum(length.(m), dims = 1)
    for r = 1:size(sm, 1)
        for c = 1:size(sm, 1)
            print(io, addspace(sm[r,c], lv[c] - length(sm[r,c]))*"   ")
        end
        println(io, "")
    end
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
#=
function Base.show(io::IO, obj::Tuple{Vararg{Tuple{Float64, Float64}}})
    for i in obj
        println(io, i)
    end
end
=#
