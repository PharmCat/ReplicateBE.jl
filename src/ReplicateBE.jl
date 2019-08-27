# ReplicateBE
# Author: Vladimir Arnautov aka PharmCat
# Copyright © 2019 Vladimir Arnautov aka PharmCat <mail@pharmcat.net>
# Licence: GNU General Public License v3.0

module ReplicateBE

using DataFrames, Distributions, StatsModels, StatsBase, ForwardDiff, LinearAlgebra, Optim

    export RBE, rbe, reml2, show, confint, contrast, lsm, emm, lmean
    import Base.show
    import StatsBase.confint
    import Statistics.var

LOG2PI = log(2π)


struct RBE
    model::ModelFrame
    factors::Array{Symbol, 1}
    β::Array{Float64, 1}
    θ0::Array{Float64, 1}
    θ::Array{Float64, 1}
    reml::Float64
    se::Array{Float64, 1}
    f::Array{Float64, 1}
    df::Array{Float64, 1}
    df2::Float64
    R::Array{Matrix{Float64},1}
    V::Array{Matrix{Float64},1}
    G::Matrix{Float64}
    C::Matrix{Float64}
    A::Matrix{Float64}              #asymptotic variance-covariance matrix ofb θ
    H::Matrix{Float64}
    X::Matrix
    Z::Matrix
    Xv::Array{Matrix{Float64},1}
    Zv::Array{Matrix{Float64},1}
    yv::Array{Array{Float64, 1},1}
    detH::Float64
    preoptim::Optim.MultivariateOptimizationResults
    optim::Optim.MultivariateOptimizationResults
end

include("show.jl")
include("utils.jl")

"""
    Mixed model fitting function
"""
function rbe(df; dvar::Symbol,
    subject::Symbol,
    formulation::Symbol,
    period::Symbol,
    sequence::Symbol,
    g_tol::Float64 = 1e-8, x_tol::Float64 = 1e-8, f_tol::Float64 = 1e-8, iterations::Int = 100)

    categorical!(df, subject);
    categorical!(df, formulation);
    categorical!(df, period);
    categorical!(df, sequence);
    sort!(df, [subject, formulation, period])
    Xf = @eval(@formula($dvar ~ $sequence + $period + $formulation))
    Zf = @eval(@formula($dvar ~ 0 + $formulation))
    MF = ModelFrame(Xf, df)
    X  = ModelMatrix(MF).m
    Z  = ModelMatrix(ModelFrame(Zf, df, contrasts = Dict(formulation => StatsModels.FullDummyCoding()))).m
    p  = rank(X)
    y  = df[:, dvar]
    Xv, Zv, yv = sortsubjects(df, subject, X, Z, y)
    n  = length(Xv)
    N  = sum(length.(yv))
    pn = length(MF.contrasts[period].levels)
    sn = length(MF.contrasts[sequence].levels)

    memc = Array{Array{Float64, 2}, 1}(undef, 4)
    memc[1] = zeros(1,2)
    memc[2] = zeros(2,2)
    memc[3] = zeros(3,2)
    memc[4] = zeros(4,2)
    memc2 = Array{Array{Float64, 2}, 1}(undef, 4)
    memc2[1] = zeros(p,1)
    memc2[2] = zeros(p,2)
    memc2[3] = zeros(p,3)
    memc2[4] = zeros(p,4)

    if size(Z)[2] != 2 error("Size random effect matrix != 2. Not implemented yet!") end
    checkdata(X, Z, Xv, Zv, y)

    qro   = qr(X)
    β     = inv(qro.R)*qro.Q'*y

    iv = initvar(df, dvar, formulation, subject)
    if iv[1] < iv[3] || iv[2] < iv[3] iv[1] = iv[2] = 2*iv[3] end
    θvec0 = [iv[3], iv[3], iv[1]-iv[3], iv[2]-iv[3], 0.001]

    G     = zeros(2, 2)
    Rv    = Array{Array{Float64,2}, 1}(undef, n)
    Vv    = Array{Array{Float64,2}, 1}(undef, n)
    iVv   = Array{Array{Float64,2}, 1}(undef, n)
    matvecz!(Rv, Zv)
    matvecz!(Vv, Zv)
    matvecz!(iVv, Zv)

    remlf(x) = -reml2!(yv, Zv, p, n, N, Xv, G, Rv, Vv, iVv, x, β, memc)

    method =LBFGS()
    #method=ConjugateGradient()
    #method = NelderMead()
    limeps=eps()

    pO = optimize(remlf, [limeps, limeps, limeps, limeps, limeps], [Inf, Inf, Inf, Inf, 1.0], θvec0, Fminbox(method), Optim.Options(g_tol = 1e-2))
    θ  = Optim.minimizer(pO)
    remlf(x) = -reml2b!(yv, Zv, p, Xv, G, Rv, Vv, iVv, x, β, memc, memc2)
    #O  = optimize(remlf, θ, method=Newton(),  g_tol=g_tol, x_tol=x_tol, f_tol=f_tol, callback = βcoef!(p, n, yv, Xv, iVv, β), allow_f_increases = true, store_trace = true, extended_trace = true, show_trace = false)
    O  = optimize(remlf, θ, method=Newton(),  g_tol=g_tol, x_tol=x_tol, f_tol=f_tol, allow_f_increases = true, store_trace = true, extended_trace = true, show_trace = false)
    θ  = Optim.minimizer(O)

    #H  = Optim.trace(O)[end].metadata["h(x)"]
    remlv = remlf(θ)
    #G     = gmat(θ[3], θ[4], θ[5])
    #Rv    = rmatvec(θ[1], θ[2], Zv)
    #Vv    = vmatvec(Zv, G, Rv)
    #iVv   = inv.(Vv)
    #βcoef!(p, n, yv, Xv, iVv, β)
    #β     = βcoef(yv, X, Xv, iVv)
    #
    if θ[5] > 1 θ[5] = 1 end
    remlf2(x) = -2*reml(yv, Zv, p, Xv, x, β)
    H         = ForwardDiff.hessian(remlf2, θ)
    dH        = det(H)
    H[5,:] .= 0
    H[:,5] .= 0

    A            = 2*pinv(H)
    se, F, df, C = ctrst(p, Xv, Zv, iVv, θ, β, A)
    df2          = N / pn - sn
    return RBE(MF, [sequence, period, formulation], β, θvec0, θ, remlv, se, F, df, df2, Rv, Vv, G, C, A, H, X, Z, Xv, Zv, yv, dH, pO, O)
end

"""
    Make X, Z mtrices and vector y for each subject;
"""
function sortsubjects(df::DataFrame, sbj::Symbol, X::Matrix, Z::Matrix, y::Vector)
    u = unique(df[:, sbj])
    Xa = Array{Array{Float64,2}, 1}(undef, 0)
    Za = Array{Array{Float64,2}, 1}(undef, 0)
    ya = Array{Array{Float64,1}, 1}(undef, 0)
    for i in u
        v = findall(x->x==i, df[:, sbj])
        Xs = Array{Float64, 1}(undef, 0)
        Zs = Array{Float64, 1}(undef, 0)
        ys = Array{Float64, 1}(undef, 0)
        for r in v
            append!(Xs, X[r, :])
            append!(Zs, Z[r, :])
            push!(ys, y[r])
        end
        push!(Xa, Matrix(reshape(Xs, size(X)[2], :)'))
        push!(Za, Matrix(reshape(Zs, size(Z)[2], :)'))
        push!(ya, ys)
    end
    return Xa, Za, ya
end
"""
    G matrix
"""
@inline function gmat(σ₁, σ₂, ρ)
    if ρ > 1.0 ρ = 1.0 end
    if σ₁ < 0.0 σ₁ = 1.0e-6 end
    if σ₂ < 0.0 σ₂ = 1.0e-6 end
    cov = sqrt(σ₁ * σ₂) * ρ
    return [σ₁ cov; cov σ₂]
end
@inline function gmat!(G::Matrix{Float64}, σ₁::Float64, σ₂::Float64, ρ::Float64)
    if ρ > 1.0 ρ = 1.0 end
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
@inline function rmat(σ, Z)
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

"""
    Return variance-covariance matrix V
"""
@inline function vmat(G, R, Z)
    V  = Z*G*Z' + R
end
@inline function vmat!(V::Matrix{Float64}, G::Matrix{Float64}, R::Matrix{Float64}, Z::Matrix{Float64}, memc)
    #copyto!(V, Z*G*Z')
    mul!(memc[size(Z)[1]], Z, G)
    mul!(V, memc[size(Z)[1]], Z')
    V .+= R
    return
end
"""
    Return set of R matrices
"""
function rmatvec(σ₁, σ₂, Zvec)
    n  = length(Zvec)
    Ra = Array{Array{Float64,2}, 1}(undef, n)
    for i = 1:n
        Ra[i] = rmat([σ₁, σ₂], Zvec[i])
    end
    return Ra
end
"""
    Return set of V matrices
"""
function vmatvec(Zvec, G, Rvec)
    n  = length(Zvec)
    Va = Array{Array{Float64,2}, 1}(undef, n)
    for i = 1:length(Zvec)
        Va[i] = vmat(G, Rvec[i], Zvec[i])
    end
    return Va
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
function cmat(Xv, Zv, iVv, θ)::Array{Float64, 2}
    p = size(Xv[1])[2]
    C = zeros(p,p)
    for i=1:length(Xv)
        C .+= Xv[i]'*iVv[i]*Xv[i]
    end
    return inv(C)
end
function βcoef(yv, X, Xv, iVv)
    p = rank(X)
    n = length(yv)
    A = zeros(p,p)
    β = zeros(p)
    for i = 1:n
        A = A + (Xv[i]'*iVv[i]*Xv[i])
        β = β .+ Xv[i]'*iVv[i]*yv[i]
    end
    return inv(A)*β
end
function βcoef!(p::Int, n::Int, yv::Array{Array{Float64, 1}, 1}, Xv::Array{Array{Float64, 2}, 1}, iVv::Array{Array{Float64, 2}, 1}, β::Array{Float64, 1})
    A = zeros(p,p)
    β0 = zeros(p)
    for i = 1:n
        A  .+= Xv[i]'*iVv[i]*Xv[i]
        β0 .+=  Xv[i]'*iVv[i]*yv[i]
    end
    copyto!(β, inv(A)*β0)
    return
end
#println("θ₁: ", θ1, " θ₂: ",  θ2,  " θ₃: ", θ3)
"""
    REML function for ForwardDiff
"""
function reml(yv, Zv, p, Xv, θvec, β)
    n = length(yv)
    N = sum(length.(yv))
    G = gmat(θvec[3], θvec[4], θvec[5])
    c  = (N-p)/2*LOG2PI
    θ1 = 0
    θ2 = 0
    θ3 = 0
    iV   = nothing
    θ2m  = zeros(p,p)
    for i = 1:n
        R   = rmat([θvec[1], θvec[2]], Zv[i])
        V   = vmat(G, R, Zv[i])
        iV  = inv(V)
        θ1  += logdet(V)
        θ2m += Xv[i]'*iV*Xv[i]
        r    = yv[i]-Xv[i]*β
        θ3  += r'*iV*r
    end
    θ2       = logdet(θ2m)
    return   -(θ1/2 + θ2/2 + θ3/2 + c)
end
function lcgf(L, Xv, Zv, θ)
    p   = size(Xv[1])[2]
    C   = zeros(p,p)
    G   = gmat(θ[3], θ[4], θ[5])
    for i=1:length(Xv)
        R   = rmat([θ[1], θ[2]], Zv[i])
        iV  = inv(vmat(G, R, Zv[i]))
        C   = C + Xv[i]'*iV*Xv[i]
    end
    return (L*inv(C)*L')[1]
end
"""
    Secondary param estimation:
    SE
    F
    DF
"""
function ctrst(p, Xv, Zv, iVv, θ, β, A)
    C     = cmat(Xv, Zv, iVv, θ)
    se    = Array{Float64, 1}(undef, p)
    F     = Array{Float64, 1}(undef, p)
    df    = Array{Float64, 1}(undef, p)
    for i = 2:p
        L    = zeros(p)
        L[i] = 1
        L    = L'
        lcl  = L*C*L'
        lclr = rank(lcl)
        se[i]   = sqrt((lcl)[1])
        #F[i]    = (L*β)'*inv(L*C*L')*(L*β)
        F[i]    = (L*β)'*inv(lcl)*(L*β)/lclr
        lclg(x) = lcgf(L, Xv, Zv, x)
        g       = ForwardDiff.gradient(lclg, θ)
        df[i]   = 2*((lcl)[1])^2/(g'*(A)*g)
    end
    return se, F, df, C
end
#-------------------------------------------------------------------------------
function checkdata(X, Z, Xv, Zv, y)
    if length(Xv) != length(Zv) error("Length Xv != Zv !!!") end
    for i = 1:length(Xv)
        if size(Xv[i])[1]  != size(Zv[i])[1] error("Row num of subject $i Xv != Zv !!!") end
        #if size(Xv[i])[1]  != 4 error("Subject observation of subject $i != 4, other designs not implemented yet!!!") end
        #if sum(Zv[i][:,1]) != 2 error("Subject $i, formulation 1, not have 2 observation, other solutions not implemented yet!!!") end
        #if sum(Zv[i][:,2]) != 2 error("Subject $i, formulation 2, not have 2 observation, other solutions not implemented yet!!!") end
    end
end


#-------------------------------------------------------------------------------
"""
    Optim reml with β
"""
function reml2b!(yv, Zv, p, Xv, G, Rv, Vv, iVv, θvec, β, memc, memc2)
    n = length(yv)
    N = sum(length.(yv))
    gmat!(G, θvec[3], θvec[4], θvec[5])
    c  = (N-p)*LOG2PI #log(2π)
    θ1 = 0
    #θ2 = 0
    θ3 = 0
    iV   = nothing
    θ2m  = zeros(p,p)
    βm   = zeros(p)
    @inbounds for i = 1:n
        rmat!(Rv[i], [θvec[1], θvec[2]], Zv[i])
        vmat!(Vv[i], G, Rv[i], Zv[i], memc)
        copyto!(iVv[i], inv(Vv[i]))
        θ1  += logdet(Vv[i])
        mul!(memc2[size(Xv[i])[1]], Xv[i]', iVv[i])
        θ2m .+= memc2[size(Xv[i])[1]]*Xv[i]
        βm  .+= memc2[size(Xv[i])[1]]*yv[i]
        #tm   = Xv[i]'*iVv[i]    #Temp matrix for Xv[i]'*iV*Xv[i] and Xv[i]'*iV*yv[i] calc
        #θ2m .+= tm*Xv[i]
        #βm  .+= tm*yv[i]
    end
    mul!(β, inv(θ2m), βm)
    for i = 1:n
        r    = yv[i] - Xv[i]*β
        θ3  += r'*iVv[i]*r
    end
    #θ2       = logdet(θ2m)
    return   -(θ1 + logdet(θ2m) + θ3 + c)
end
"""
    Optim reml without β
"""
function reml2!(yv::S, Zv::T, p::Int, n::Int, N::Int, Xv::T, G::Array{Float64, 2}, Rv::T, Vv::T, iVv::T, θvec::Array{Float64, 1}, β::Array{Float64, 1}, memc)::Float64 where T <: Array{Array{Float64, 2}, 1} where S <: Array{Array{Float64, 1}, 1}

    gmat!(G, θvec[3], θvec[4], θvec[5])
    c  = (N-p)*LOG2PI
    θ1 = 0
    #θ2 = 0
    θ3 = 0
    θ2m  = zeros(p,p)
    for i = 1:n
        rmat!(Rv[i], [θvec[1], θvec[2]], Zv[i])
        vmat!(Vv[i], G, Rv[i], Zv[i], memc)
        copyto!(iVv[i], inv(Vv[i]))
        θ1  += logdet(Vv[i])
        θ2m .+= Xv[i]'*iVv[i]*Xv[i]
        r    = yv[i]-Xv[i]*β
        θ3  .+= r'*iVv[i]*r
    end
    #θ2       = logdet(θ2m)
    return   -(θ1 + logdet(θ2m) + θ3 + c)
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

end # module
