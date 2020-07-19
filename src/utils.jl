#return contrast table
"""
    contrast(rbe::RBE, L::Matrix; name = "Contrast", memopt = true)::ContrastTable

Return contrast table for L matrix. Table include:
* F
* NumDF
* DF
* P|f|

```math
F = \\frac{\\beta'L'(LCL')^{-1}L\\beta}{rank(LCL')}
```

where

```math
C = (\\sum_{i=1}^{n} X_i'V_i^{-1}X_i)^{-1}
```

DF for one-dimetion case:

```math
df = \\frac{2(LCL')^{2}}{g'Ag}
```

where ``A = 2H^{-1}``

where ``g = \\triangledown _{\\theta}(LC_{\\theta}L')``

DF for multi-dimention case see Schaalje et al 2002.

p value calculated with:

```julia
pval    = ccdf(FDist(numdf, df), F)
```
"""
function contrast(rbe::RBE, L::Matrix;  name = "Contrast", memopt = true)::ContrastTable
    F, ndf, df, pval = contrastvec(rbe.data, rbe.result, L)
    return ContrastTable([name], [F], [ndf], [df], [pval])

end

"""
    estimate(rbe::RBE, L::Matrix; df = :sat, name = "Estimate", memopt = true, alpha = 0.05)

Return estimate table for L 1xp matrix. Table include:
* estimate value
* SE
* DF
* t
* P|t|
* CI Upper
* CI Lower

```math
estimate = L\\beta
```

```math
SE = \\sqrt{LCL'}
```

```math
t = estimate / SE
```

For ```df = :sat```:

```math
df = \\frac{2(LCL')^{2}}{g'Ag}
```

where

```math
C = (\\sum_{i=1}^{n} X_i'V_i^{-1}X_i)^{-1}
```

where ``A = 2H^{-1}``

where ``g = \\triangledown _{\\theta}(LC_{\\theta}L')``

For ```df = :cont``` (contain):

```math
df = N - rank(ZX)
```

CI estimate is:

```math
CI = estimate ± t(alpha, df)*SE
```

Example of L matrix if length of fixed effect vector is 6, estimate for 4-th value:

```julia
L = [0 0 0 1 0 0]
```
"""
function estimate(rbe::RBE, L::Matrix; df = :sat, name = "Estimate", memopt = true, alpha = 0.05)::EstimateTable
    est, se, t = estimatevec(rbe.data, rbe.result, L)
    if df == :sat
        df      = sattdf(rbe.data, rbe.result.gradc, rbe.result.A, rbe.result.C, L, L*rbe.result.C*L')
    elseif df == :cont
        df      = rbe.design.df3
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
    vec = Vector{Int}(undef, 0)
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
#
#-------------------------------------------------------------------------------
function checkdata(X::Matrix, Z::Matrix, Xv::Vector, Zv::Vector, y::Vector)
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
    lv = maximum(length.(sm), dims = 1)
    for r = 1:size(sm, 1)
        for c = 1:size(sm, 2)
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
function intravar(rbe::RBE)
    terms = rbe.rmodel.f.rhs.terms[2].contrasts.termnames
    θ     = theta(rbe)
    return Dict(terms .=> θ[1:2])
end
