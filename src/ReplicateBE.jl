# ReplicateBE
# Author: Vladimir Arnautov aka PharmCat
# Copyright © 2019 Vladimir Arnautov aka PharmCat <mail@pharmcat.net>
# Licence: GNU General Public License v3.0

module ReplicateBE

using DataFrames, Distributions, StatsModels, StatsBase, ForwardDiff, LinearAlgebra, Optim, Memoize, TimerOutputs

    export RBE, rbe, reml2, show, confint, contrast, lsm, emm, lmean
    import Base.show
    import StatsBase.confint
    import Statistics.var

const LOG2PI = log(2π)
const MEMOPT = true

struct RBE
    model::ModelFrame               #Model frame
    rmodel::ModelFrame              #Random effect model
    factors::Array{Symbol, 1}       #Factor list
    β::Array{Float64, 1}            #β coefficients (fixed effect)
    θ0::Array{Float64, 1}           #Initial variance paramethers
    θ::Array{Float64, 1}            #Final variance paramethers
    reml::Float64                   #-2REML
    se::Array{Float64, 1}           #SE for each β level
    f::Array{Float64, 1}            #F for each β level
    df::Array{Float64, 1}           #DF (degree of freedom) for each β level (Satterthwaite)
    df2::Float64                    #DF N / pn - sn
    R::Array{Matrix{Float64},1}     #R matrices for each subject
    V::Array{Matrix{Float64},1}     #V matrices for each subject
    G::Matrix{Float64}              #G matrix
    C::Matrix{Float64}              #C var(β) p×p variance-covariance matrix
    A::Matrix{Float64}              #asymptotic variance-covariance matrix ofb θ
    H::Matrix{Float64}              #Hessian matrix
    X::Matrix                       #Matrix for fixed effects
    Z::Matrix                       #Matrix for random effects
    Xv::Array{Matrix{Float64},1}    #X matrices for each subject
    Zv::Array{Matrix{Float64},1}    #Z matrices for each subject
    yv::Array{Array{Float64, 1},1}  #responce vectors for each subject
    detH::Float64                   #Hessian determinant
    preoptim::Optim.MultivariateOptimizationResults        #Pre-optimization result object
    optim::Optim.MultivariateOptimizationResults           #Optimization result object
end

#=
struct MemAlloc
    ZG
    β
    θ2
    r
    tr
    tX
    tXiV
    Xvβ
    triVv
end
=#

#Find by Symbol
function findterm(MF::ModelFrame, symbol::Symbol)::Int
    l = length(MF.f.rhs.terms)
    for i = 1:l
        if isa(MF.f.rhs.terms[i], InterceptTerm) continue end
        if MF.f.rhs.terms[i].sym == symbol return i end
    end
    return 0
end
#Return length by Symbol
function termmodellen(MF::ModelFrame, symbol::Symbol)::Int
    id = findterm(MF, symbol)
    return length(MF.f.rhs.terms[id].contrasts.termnames)
end

include("show.jl")
include("utils.jl")
include("memalloc.jl")
include("deprecated.jl")
#-------------------------------------------------------------------------------
"""
    Mixed model fitting function
"""
function rbe(df; dvar::Symbol,
    subject::Symbol,
    formulation::Symbol,
    period::Symbol,
    sequence::Symbol,
    g_tol::Float64 = 1e-8, x_tol::Float64 = 0.0, f_tol::Float64 = 0.0, iterations::Int = 100,
    store_trace = false, extended_trace = false, show_trace = false, memopt = true)

    to = TimerOutput()

    categorical!(df, subject);
    categorical!(df, formulation);
    categorical!(df, period);
    categorical!(df, sequence);
    @timeit to "sort" sort!(df, [subject, formulation, period])
    Xf = @eval(@formula($dvar ~ $sequence + $period + $formulation))
    Zf = @eval(@formula($dvar ~ 0 + $formulation))
    MF = ModelFrame(Xf, df)
    RMF = ModelFrame(Zf, df, contrasts = Dict(formulation => StatsModels.FullDummyCoding()))
    X  = ModelMatrix(MF).m
    Z  = ModelMatrix(RMF).m
    p  = rank(X)
    y  = df[:, dvar]                                                            #Dependent variable

    #Make pre located arrays with matrices for each subject
    @timeit to "sortsubj" Xv, Zv, yv = sortsubjects(df, subject, X, Z, y)
    n  = length(Xv)
    N  = sum(length.(yv))
    pn = termmodellen(MF, period)
    sn = termmodellen(MF, sequence)
    #pn = length(MF.contrasts[period].levels)
    #sn = length(MF.contrasts[sequence].levels)

    #Memory pre-allocation arrays for matrix computations
    memc, memc2, memc3, memc4 = memcalloc(p, 2, yv)
    #Check data
    @timeit to "check" checkdata(X, Z, Xv, Zv, y)

    #Calculate initial fixed parameters
    qro   = qr(X)
    @timeit to "lm"  β     = inv(qro.R)*qro.Q'*y

    #Calculate initial variance
    iv = initvar(df, dvar, formulation, subject)
    if iv[1] < iv[3] || iv[2] < iv[3] iv[1] = iv[2] = 2*iv[3] end
    θvec0 = [iv[3], iv[3], iv[1]-iv[3], iv[2]-iv[3], 0.501]

    #Prelocatiom for G, R, V, V⁻¹ matrices
    G     = zeros(2, 2)
    Rv    = Array{Array{Float64,2}, 1}(undef, n)
    Vv    = Array{Array{Float64,2}, 1}(undef, n)
    iVv   = Array{Array{Float64,2}, 1}(undef, n)
    matvecz!(Rv, Zv)
    matvecz!(Vv, Zv)
    matvecz!(iVv, Zv)

    #First step optimization (pre-optimization)
    od = OnceDifferentiable(x -> -2*reml(yv, Zv, p, Xv, x, β), θvec0; autodiff = :forward)
    #remlf(x)   = -reml2!(yv, Zv, p, n, N, Xv, G, Rv, Vv, iVv, x, β, memc, memc2, memc3, memc4)
    #remlf(x)  = -reml2b!(yv, Zv, p, Xv, G, Rv, Vv, iVv, x, β, memc, memc2, memc3, memc4)
    method = LBFGS()
    #method = ConjugateGradient()
    #method = NelderMead()

    limeps=eps()
    @timeit to "o1" pO = optimize(od, [limeps, limeps, limeps, limeps, limeps], [Inf, Inf, Inf, Inf, 1.0], θvec0,  Fminbox(method), Optim.Options(g_tol = 1e-1))
    #pO = optimize(remlf,  [limeps, limeps, limeps, limeps, limeps], [Inf, Inf, Inf, Inf, 1.0], θvec0,  Fminbox(method), Optim.Options(g_tol = 1e-3))
    θ  = Optim.minimizer(pO)

    #Final optimization
    #Provide gradient function for Optim
    #Not used yet
    #g!(storage, θx) = copyto!(storage, ForwardDiff.gradient(x -> -2*reml(yv, Zv, p, Xv, x, β), θx))
    #REML function for optimization
    td = TwiceDifferentiable(x -> -2*remlb(yv, Zv, p, Xv, x, β), θvec0; autodiff = :forward)
    #remlfb(x) = -reml2b!(yv, Zv, p, n, N, Xv, G, Rv, Vv, iVv, x, β, memc, memc2, memc3, memc4)
    #@timeit to "o2" O  = optimize(remlfb, θ, method=Newton(),  g_tol=g_tol, x_tol=x_tol, f_tol=f_tol, allow_f_increases = true, store_trace = store_trace, extended_trace = extended_trace, show_trace = show_trace)
    @timeit to "o2" O  = optimize(td, θ, method=Newton(),  g_tol=g_tol, x_tol=x_tol, f_tol=f_tol, allow_f_increases = true, store_trace = store_trace, extended_trace = extended_trace, show_trace = show_trace)
    θ  = Optim.minimizer(O)

    #Get reml
    #remlv = -2*reml(yv, Zv, p, Xv, θ, β)
    @timeit to "remlv"  remlv = -reml2b!(yv, Zv, p, n, N, Xv, G, Rv, Vv, iVv, θ, β, memc, memc2, memc3, memc4)

    #θ[5] can not be more than 1.0
    if θ[5] > 1 θ[5] = 1 end

    #Get Hessian matrix (H) with ForwardDiff
    @timeit to "H" H         = ForwardDiff.hessian(x -> -2*reml(yv, Zv, p, Xv, x, β), θ)
    dH        = det(H)
    H[5,:] .= 0
    H[:,5] .= 0

    #Secondary parameters calculation
    A            = 2*pinv(H)
    @timeit to "etc" se, F, df, C = ctrst(p, Xv, Zv, iVv, θ, β, A)
    df2          = N / pn - sn
    #println(to)                               #!!!should be checked!!!
    return RBE(MF, RMF, [sequence, period, formulation], β, θvec0, θ, remlv, se, F, df, df2, Rv, Vv, G, C, A, H, X, Z, Xv, Zv, yv, dH, pO, O)
end #END OF rbe()
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
"""
    Make X, Z mtrices and vector y for each subject;
"""
function sortsubjects(df::DataFrame, sbj::Symbol, X::Matrix, Z::Matrix, y::Vector)
    u = unique(df[:, sbj])
    Xa = Array{Array{Float64,2}, 1}(undef, length(u))
    Za = Array{Array{Float64,2}, 1}(undef, length(u))
    ya = Array{Array{Float64,1}, 1}(undef, length(u))
    for i = 1:length(u)
        v = findall(x->x==u[i], df[:, sbj])
        Xs = Array{Float64, 1}(undef, 0)
        Zs = Array{Float64, 1}(undef, 0)
        ys = Array{Float64, 1}(undef, 0)
        for r in v
            append!(Xs, X[r, :])
            append!(Zs, Z[r, :])
            push!(ys, y[r])
        end
        Xa[i] = Matrix(reshape(Xs, size(X)[2], :)')
        Za[i] = Matrix(reshape(Zs, size(Z)[2], :)')
        ya[i] = ys
    end
    for i = 1:length(u)
        for c = 1:length(u)
            if Za[i] == Za[c] && Za[i] !== Za[c] Za[i] = Za[c] end
            if Xa[i] == Xa[c] && Xa[i] !== Xa[c] Xa[i] = Xa[c] end
        end
    end
    return Xa, Za, ya
end
"""
    G matrix
"""
@inline function gmat(σ₁::S, σ₂::T, ρ::U)::Matrix where S <: Real where T <: Real where U <: Real
    if ρ > 1.0 ρ = 1.0 end
    if ρ < 0.0 ρ = 0.0 end
    if σ₁ < 0.0 σ₁ = 1.0e-6 end
    if σ₂ < 0.0 σ₂ = 1.0e-6 end
    cov = sqrt(σ₁ * σ₂) * ρ
    return [σ₁ cov; cov σ₂]
end
@inline function gmat!(G::Matrix{Float64}, σ₁::Float64, σ₂::Float64, ρ::Float64)
    if ρ > 1.0 ρ = 1.0 end
    if ρ < 0.0 ρ = 0.0 end
    if σ₁ < 0.0 σ₁ = 1.0e-6 end
    if σ₂ < 0.0 σ₂ = 1.0e-6 end
    G[1, 1] = σ₁
    G[2, 2] = σ₂
    G[1, 2] = G[2, 1] = sqrt(σ₁ * σ₂) * ρ
    return
end

"""
    R matrix (ForwardDiff+)
"""
@inline function rmat(σ::Vector{S}, Z::Matrix{T})::Matrix where S <: Real where T <: Real
    if σ[1] < 0.0 σ[1] = 1.0e-6 end
    if σ[2] < 0.0 σ[2] = 1.0e-6 end
    return Matrix(Diagonal((Z*σ)[:,1]))
end
@inline function rmat!(R::Matrix{Float64}, σ::Array{Float64, 1}, Z::Matrix{Float64})
    if σ[1] < 0.0 σ[1] = 1.0e-6 end
    if σ[2] < 0.0 σ[2] = 1.0e-6 end
    copyto!(R, Matrix(Diagonal((Z*σ)[:,1])))
    return
end
@memoize function memrmat(σ::Vector, Z::Matrix)::Matrix
    return rmat(σ, Z)
end
@memoize function memzgz(G::Matrix, Z::Matrix)::Matrix
    return Z*G*Z'
end
@memoize function memvmat(ZGZ::Matrix, R::Matrix)::Matrix
    return ZGZ + R
end
@memoize function meminv(m::Matrix)::Matrix
    return inv(m)
end
"""
    Return variance-covariance matrix V
"""
@inline function vmat(G, R, Z)
    V  = Z*G*Z' + R
    return V
end
@inline function vmat!(V::Matrix{Float64}, G::Matrix{Float64}, R::Matrix{Float64}, Z::Matrix{Float64}, memc)
    #copyto!(V, Z*G*Z')
    mul!(memc[size(Z)[1]], Z, G)
    mul!(V, memc[size(Z)[1]], Z')
    V .+= R
    return
end
"""
    Return set of zero feeled matrices
"""
function matvecz!(M, Zv)
    n  = length(Zv)
    for i = 1:n
        M[i] = zeros(size(Zv[i])[1], size(Zv[i])[1])
    end
    return
end
"""
    Return C matrix
    var(β) p×p variance-covariance matrix
"""
@inline function cmat(Xv, Zv, iVv, θ)::Array{Float64, 2}
    p = size(Xv[1])[2]
    C = zeros(p,p)
    for i=1:length(Xv)
        C .+= Xv[i]'*iVv[i]*Xv[i]
    end
    return inv(C)
end
#println("θ₁: ", θ1, " θ₂: ",  θ2,  " θ₃: ", θ3)


"""
    REML function for ForwardDiff
"""
function reml(yv, Zv, p, Xv, θvec, β; memopt::Bool = true)
    maxobs    = maximum(length.(yv))
    #some memory optimizations to reduse allocations
    mXviV     = Array{Array{eltype(θvec), 2}, 1}(undef, maxobs)
    mXviVXv   = zeros(promote_type(Float64, eltype(θvec)), p, p)
    for i = 1:maxobs
        mXviV[i] =  zeros(promote_type(Float64, eltype(θvec)), p, i)
    end
    #---------------------------------------------------------------------------
    n         = length(yv)
    N         = sum(length.(yv))
    G         = gmat(θvec[3], θvec[4], θvec[5])
    c         = (N-p)*LOG2PI
    θ1        = 0
    θ2        = zeros(promote_type(Float64, eltype(θvec)), p, p)
    θ3        = 0
    iV        = nothing
    θr   = [θvec[1], θvec[2]]
    for i = 1:n
        if MEMOPT && memopt
            R    = memrmat(θr, Zv[i])
            V    = memvmat(memzgz(G, Zv[i]), R)
            iV   = meminv(V)
        else
            R    = rmat(θr, Zv[i])
            V    = vmat(G, R, Zv[i])
            iV   = inv(V)
        end
        θ1  += logdet(V)
        #-----------------------------------------------------------------------
        #θ2 += Xv[i]'*iV*Xv[i]
        mul!(mXviV[size(Xv[i])[1]], Xv[i]', iV)
        mul!(mXviVXv, mXviV[size(Xv[i])[1]], Xv[i])
        θ2  += mXviVXv
        #-----------------------------------------------------------------------
        r    = yv[i]-Xv[i]*β
        θ3  += r'*iV*r
    end
    return   -(θ1 + logdet(θ2) + θ3 + c)/2
end

function remlb(yv, Zv, p, Xv, θvec, β; memopt::Bool = true)
    maxobs    = maximum(length.(yv))
    #some memory optimizations to reduse allocations
    mXviV     = Array{Array{eltype(θvec), 2}, 1}(undef, maxobs)
    mXviVXv   = zeros(promote_type(Float64, eltype(θvec)), p, p)
    for i = 1:maxobs
        mXviV[i] =  zeros(promote_type(Float64, eltype(θvec)), p, i)
    end
    #---------------------------------------------------------------------------
    n         = length(yv)
    N         = sum(length.(yv))
    G         = gmat(θvec[3], θvec[4], θvec[5])
    iVv       = Array{Array{eltype(θvec), 2}, 1}(undef, n)
    c         = (N-p)*LOG2PI
    θ1        = 0
    θ2        = zeros(promote_type(Float64, eltype(θvec)), p, p)
    θ3        = 0
    iV        = nothing
    βm        = zeros(promote_type(Float64, eltype(θvec)), p)
    βt        = zeros(promote_type(Float64, eltype(θvec)), p)
    θr        = [θvec[1], θvec[2]]
    for i = 1:n
        if MEMOPT && memopt
            R    = memrmat(θr, Zv[i])
            V    = memvmat(memzgz(G, Zv[i]), R)
            iVv[i]   = meminv(V)
        else
            R    = rmat(θr, Zv[i])
            V    = vmat(G, R, Zv[i])
            iVv[i]   = inv(V)
        end
        θ1  += logdet(V)
        #-----------------------------------------------------------------------
        #θ2 += Xv[i]'*iV*Xv[i]
        mul!(mXviV[size(Xv[i])[1]], Xv[i]', iVv[i])
        mul!(mXviVXv, mXviV[size(Xv[i])[1]], Xv[i])
        θ2  += mXviVXv
        βm  .+= mXviV[size(Xv[i])[1]]*yv[i]
        #-----------------------------------------------------------------------
        #tm   = Xv[i]'*iVv[i]    #Temp matrix for Xv[i]'*iV*Xv[i] and Xv[i]'*iV*yv[i] calc
        #θ2m .+= tm*Xv[i]
        #βm  .+= tm*yv[i]
    end

    mul!(βt, inv(θ2), βm)

    for i = 1:n
        r    = yv[i] - Xv[i]*βt
        θ3  += r'*iVv[i]*r
    end

    return   -(θ1 + logdet(θ2) + θ3 + c)/2
end
"""
    Satterthwaite DF gradient function
"""
function lclgf(L, Lt, Xv, Zv, θ)
    p   = size(Xv[1])[2]
    G   = gmat(θ[3], θ[4], θ[5])
    C   = zeros(promote_type(Float64, eltype(θ)), p, p)
    θr  = [θ[1], θ[2]]
    for i=1:length(Xv)
        iV   = meminv(memvmat(memzgz(G, Zv[i]), memrmat(θr, Zv[i])))
        #R   = rmat(θr, Zv[i])
        #iV  = inv(vmat(G, R, Zv[i]))
        C  += Xv[i]'*iV*Xv[i]
    end
    return (L*inv(C)*Lt)[1]
end
"""
    Secondary param estimation:
    SE
    F
    DF
    C
"""
function ctrst(p, Xv, Zv, iVv, θ, β, A)
    C     = cmat(Xv, Zv, iVv, θ)
    se    = Array{Float64, 1}(undef, p)
    F     = Array{Float64, 1}(undef, p)
    df    = Array{Float64, 1}(undef, p)
    for i = 1:p
        L    = zeros(1, p)
        L[i]    = 1
        Lt      = L'
        lcl     = L*C*Lt                         #lcl     = L*C*L'
        lclr    = rank(lcl)
        se[i]   = sqrt((lcl)[1])
        Lβ      = L*β
        F[i]    = Lβ'*inv(lcl)*Lβ/lclr           #F[i]    = (L*β)'*inv(L*C*L')*(L*β)
        lclg(x) = lclgf(L, Lt, Xv, Zv, x)
        g       = ForwardDiff.gradient(lclg, θ)
        df[i]   = 2*((lcl)[1])^2/(g'*(A)*g)
        #LinearAlgebra.eigen(L*C*L')
    end
    return se, F, df, C
end
#-------------------------------------------------------------------------------
#             REML FOR OPT ALGORITHM
#-------------------------------------------------------------------------------
"""
    REML with β final update
"""
function reml2b!(yv::S, Zv::T, p::Int, n::Int, N::Int, Xv::T, G::Array{Float64, 2}, Rv::T, Vv::T, iVv::T, θvec::Array{Float64, 1}, β::Array{Float64, 1}, memc, memc2, memc3, memc4)::Float64 where T <: Array{Array{Float64, 2}, 1} where S <: Array{Array{Float64, 1}, 1}

    gmat!(G, θvec[3], θvec[4], θvec[5])
    c  = (N-p)*LOG2PI #log(2π)
    θ1 = 0
    θ2  = zeros(p, p)
    θ3 = 0
    iV   = nothing
    fill!(memc4, 0)
    #θ2m  = zeros(p,p)
    βm   = zeros(p)
    θr   = [θvec[1], θvec[2]]
    @inbounds for i = 1:n
        rmat!(Rv[i], θr, Zv[i])
        vmat!(Vv[i], G, Rv[i], Zv[i], memc)
        copyto!(iVv[i], inv(Vv[i]))

        #Rv[i] = memrmat(θr, Zv[i])
        #zgz   = lmemzgz(G,  Zv[i])
        #Vv[i] = memvmat(zgz, Rv[i])
        #iVv[i]= meminv(Vv[i])

        θ1  += logdet(Vv[i])
        mul!(memc2[size(Xv[i])[1]], Xv[i]', iVv[i])
        θ2    .+= memc2[size(Xv[i])[1]]*Xv[i]
        βm    .+= memc2[size(Xv[i])[1]]*yv[i]

        #ToDo
        #mul!(memcX[length(yv[i])], memc2[size(Xv[i])[1]], yv[i])
        #βm  .+= memcX[length(yv[i])]
        #tm   = Xv[i]'*iVv[i]    #Temp matrix for Xv[i]'*iV*Xv[i] and Xv[i]'*iV*yv[i] calc
        #θ2m .+= tm*Xv[i]
        #βm  .+= tm*yv[i]
    end
    mul!(β, inv(θ2), βm)
    for i = 1:n
        copyto!(memc3[length(yv[i])], yv[i])
        memc3[length(yv[i])] .-= Xv[i]*β
        θ3  += memc3[length(yv[i])]'*iVv[i]*memc3[length(yv[i])]
        #Same:
        #r    = yv[i] - Xv[i]*β
        #θ3  += r'*iVv[i]*r
    end
    #θ2       = logdet(θ2m)
    return   -(θ1 + logdet(θ2) + θ3 + c)
end
#-------------------------------------------------------------------------------
"""
    Initial variance computation
"""
function initvar(df, dv, fac, sbj)
    u  = unique(df[:, sbj])
    f  = unique(df[:, fac])
    fv = Array{Float64, 1}(undef, 0)
    for i in f
        fvd = Array{Float64, 1}(undef, 0)
        v = findall(x->x==i, df[:, fac])
        for r in v
            push!(fvd, df[r, dv])
        end
        push!(fv, var(fvd))
    end
    sv = Array{Float64, 1}(undef, 0)
    for i in u
        fvd = Array{Float64, 1}(undef, 0)
        v   = findall(x->x==i, df[:, sbj])
        for r in v
            push!(fvd, df[r, dv])
        end
        push!(sv, var(fvd))
    end
    push!(fv, mean(sv))
    return fv
end

#-------------------------------------------------------------------------------
end # module
