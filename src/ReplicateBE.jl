# ReplicateBE
# Author: Vladimir Arnautov aka PharmCat
# Copyright © 2019 Vladimir Arnautov aka PharmCat <mail@pharmcat.net>
# Licence: GNU General Public License v3.0

module ReplicateBE

using DataFrames, Distributions, StatsModels, ForwardDiff, LinearAlgebra, Optim

    export rbe
    import Base.show
    import Statistics.var

struct RBE
    model::ModelFrame
    factors::Array{Symbol, 1}
    β::Array{Float64, 1}
    θ::Array{Float64, 1}
    reml::Float64
    se::Array{Float64, 1}
    F::Array{Float64, 1}
    DF::Array{Float64, 1}
    R::Array{Matrix{Float64},1}
    V::Array{Matrix{Float64},1}
    G::Matrix{Float64}
    A::Matrix{Float64}
    H::Matrix{Float64}
    detH::Float64
    optim::Optim.MultivariateOptimizationResults
end

include("show.jl")

"""
    Mixed model fitting function
"""
function rbe(df; dvar::Symbol,
    subject::Symbol,
    formulation::Symbol,
    period::Symbol,
    sequence::Symbol,
    g_tol::Float64 = 1e-8, x_tol::Float64 = 0.0, f_tol::Float64 = 0.0)

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
    if size(Z)[2] != 2 error("Size random effect matrix != 2. Not implemented yet!") end

    y = df[:, dvar]
    Xv, Zv, yv = sortsubjects(df, subject, X, Z, y)
    checkdata(X, Z, Xv, Zv, y)

    iv = initvar(df, dvar, formulation, subject)
    if iv[1] < iv[3] || iv[2] < iv[3] iv[1] = iv[2] = 2*iv[3] end
    θvec0 = [iv[3], iv[3], iv[1]-iv[3], iv[2]-iv[3], 0.001]

    remlf(x) = -reml2(yv, Zv, p, Xv, x)

    O  = optimize(remlf, θvec0, method=Newton(),  g_tol=g_tol, x_tol=x_tol, f_tol=f_tol, allow_f_increases = true, store_trace = true, extended_trace = true, show_trace = false)
    θ  = Optim.minimizer(O)
    H  = Optim.trace(O)[end].metadata["h(x)"]
    remlv = remlf(θ)
    G     = gmat(θ[3], θ[4], θ[5])
    Rv    = rmatvec(θ[1], θ[2], Zv)
    Vv    = vmatvec(Zv, G, Rv)
    iVv   = inv.(Vv)
    β     = βcoef(yv, X, Xv, iVv)

    if θ[5] > 1 θ[5] = 1 end
    remlf2(x) = -2*reml(yv, Zv, p, Xv, x, β)
    H         = ForwardDiff.hessian(remlf2, θ)
    dH        = det(H)
    H[5,:] .= 0
    H[:,5] .= 0

    A         = 2*pinv(H)
    se, F, df = ctrst(p, Xv, Zv, θ, β, A)
    return RBE(MF, [sequence, period, formulation], β, θ, remlv, se, F, df, Rv, Vv, G, A, H, dH, O)
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
function gmat(σ₁, σ₂, ρ)
    if ρ > 1.0 ρ = 1.0 end
    if σ₁ < 0.0 σ₁ = 1.0e-6 end
    if σ₂ < 0.0 σ₂ = 1.0e-6 end
    cov = sqrt(σ₁ * σ₂) * ρ
    return [σ₁ cov; cov σ₂]
end

"""
    R matrix (ForwardDiff+)
"""
@inline function rmat(σ, Z)
    if σ[1] < 0.0 σ[1] = 1.0e-6 end
    if σ[2] < 0.0 σ[2] = 1.0e-6 end
    return Matrix(Diagonal((Z*σ)[:,1]))
end

#=
function rmat(σ₁, σ₂, Z)
    if σ₁ < 0.0 σ₁ = 1.0e-6 end
    if σ₂ < 0.0 σ₂ = 1.0e-6 end
    M = zeros(4, 4)
    for i = 1:4
        M[i, i] = Z[i,:]'*[σ₁, σ₂]
    end
    return M
end
=#
#=
function rmat2(σ, Z)
    if σ[1] < 0.0 σ[1] = 1.0e-6 end
    if σ[2] < 0.0 σ[2] = 1.0e-6 end
    return Matrix(Diagonal((Z*σ)[:,1]))
end
=#
"""
    Return variance-covariance matrix V
"""
function cov(G, R, Z)
    V  = Z*G*Z' + R
end
function rmatvec(σ₁, σ₂, Zvec)
    n  = length(Zvec)
    Ra = Array{Array{Float64,2}, 1}(undef, n)
    for i = 1:n
        Ra[i] = rmat([σ₁, σ₂], Zvec[i])
    end
    return Ra
end
function vmatvec(Zvec, G, Rvec)
    Va = Array{Array{Float64,2}, 1}(undef, 0)
    for i = 1:length(Zvec)
        push!(Va, cov(G, Rvec[i], Zvec[i]))
    end
    return Va
end
function cmat(Xv, Zv, θ)
    p = size(Xv[1])[2]
    C = zeros(p,p)
    G   = gmat(θ[3], θ[4], θ[5])
    Rv  = rmatvec(θ[1], θ[2], Zv)
    Vv  = vmatvec(Zv, G, Rv)
    iVv = inv.(Vv)
    for i=1:length(Xv)
        C = C + Xv[i]'*iVv[i]*Xv[i]
    end
    return pinv(C)
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
#=
function reml(yv, Zv, X, Xv, θvec)
    n = length(yv)
    N = sum(length.(yv))
    G = gmat(θvec[3], θvec[4], θvec[5])
    p = rank(X)
    c  = (N-p)/2*log(2π)
    θ1 = 0
    θ2 = 0
    θ3 = 0
    θ2m  = zeros(p,p)
    Rv   = rmatvec(θvec[1], θvec[2], Zv)
    Vv   = vmatvec(Zv, G, Rv)
    iVv  = inv.(Vv)
    β    = βcoef(yv, X, Xv, iVv)
    for i = 1:n
        θ1  += log(det(Vv[i]))
        θ2m += Xv[i]'*iVv[i]*Xv[i]
        r    = yv[i]-Xv[i]*β
        θ3  += r'*iVv[i]*r
    end
    θ2       = log(det(θ2m))
    #println("θ₁: ", θ1, " θ₂: ",  θ2,  " θ₃: ", θ3)
    return   -(θ1/2 + θ2/2 + θ3/2 + c)
end
=#
function reml(yv, Zv, p, Xv, θvec, β)
    n = length(yv)
    N = sum(length.(yv))
    G = gmat(θvec[3], θvec[4], θvec[5])
    c  = (N-p)/2*log(2π)
    θ1 = 0
    θ2 = 0
    θ3 = 0
    iV   = nothing
    θ2m  = zeros(p,p)
    for i = 1:n
        R   = rmat([θvec[1], θvec[2]], Zv[i])
        V   = cov(G, R, Zv[i])
        iV  = inv(V)
        θ1  += logdet(V)
        θ2m += Xv[i]'*iV *Xv[i]
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
        iVv = inv(cov(G, R, Zv[i]))
        C   = C + Xv[i]'*iVv*Xv[i]
    end
    return (L*inv(C)*L')[1]
end
function ctrst(p, Xv, Zv, θ, β, A)
    C     = cmat(Xv, Zv, θ)
    se    = Array{Float64, 1}(undef, 0)
    F     = Array{Float64, 1}(undef, 0)
    df    = Array{Float64, 1}(undef, 0)
    for i = 2:p
        L    = zeros(p)
        L[i] = 1
        L    = L'
        push!(se, sqrt((L*C*L')[1]))
        push!(F, (L*β)'*inv(L*C*L')*(L*β))
        lclg(x) = lcgf(L, Xv, Zv, x)
        g       = ForwardDiff.gradient(lclg, θ)
        push!(df, 2*((L*C*L')[1])^2/(g'*(A)*g))
    end
    return se, F, df
end

function checkdata(X, Z, Xv, Zv, y)
    if length(Xv) != length(Zv) error("Length Xv != Zv !!!") end
    for i = 1:length(Xv)
        if size(Xv[i])[1]  != size(Zv[i])[1] error("Row num of subject $i Xv != Zv !!!") end
        if size(Xv[i])[1]  != 4 error("Subject observation of subject $i != 4, other designs not implemented yet!!!") end
        if sum(Zv[i][:,1]) != 2 error("Subject $i, formulation 1, not have 2 observation, other solutions not implemented yet!!!") end
        if sum(Zv[i][:,2]) != 2 error("Subject $i, formulation 2, not have 2 observation, other solutions not implemented yet!!!") end
    end
end


#-------------------------------------------------------------------------------
function reml2(yv, Zv, p, Xv, θvec)
    n = length(yv)
    N = sum(length.(yv))
    G = gmat(θvec[3], θvec[4], θvec[5])
    c  = (N-p)*log(2π)
    θ1 = 0
    θ2 = 0
    θ3 = 0
    iV   = nothing
    θ2m  = zeros(p,p)
    βm   = zeros(p)
    for i = 1:n
        R   = rmat([θvec[1], θvec[2]], Zv[i])
        V   = cov(G, R, Zv[i])
        iV  = inv(V)
        θ1  += logdet(V)
        tm   = Xv[i]'*iV    #Temp matrix for Xv[i]'*iV*Xv[i] and Xv[i]'*iV*yv[i] calc
        θ2m += tm*Xv[i]
        βm  += tm*yv[i]
    end
    β        = inv(θ2m)*βm
    for i = 1:n
        r    = yv[i]-Xv[i]*β
        θ3  += r'*iV*r
    end
    θ2       = logdet(θ2m)
    return   -(θ1 + θ2 + θ3 + c)
end
#=
function βcoef2(yv, p, Xv, Zv, θvec)
    n   = length(yv)
    A   = zeros(p,p)
    β   = zeros(p)
    G   = gmat(θvec[3], θvec[4], θvec[5])
    for i = 1:n
        R   = rmat([θvec[1], θvec[2]], Zv[i])
        iV = inv(cov(G, R, Zv[i]))
        A = A + (Xv[i]'*iV*Xv[i])
        β = β .+ Xv[i]'*iV*yv[i]
    end
    return inv(A)*β
end
=#
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
