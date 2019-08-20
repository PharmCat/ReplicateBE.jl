# ReplicateBE
# Author: Vladimir Arnautov aka PharmCat
# Copyright © 2019 Vladimir Arnautov aka PharmCat <mail@pharmcat.net>
# Licence: GNU General Public License v3.0

module ReplicateBE

using DataFrames, Distributions, StatsModels, ForwardDiff, LinearAlgebra, Optim

    export rbe
    import Base.show

struct RBE
    model
    factors
    β
    θ
    reml
    se
    F
    DF
    R
    V
    G
    A
    H
    detH
end

#=
function Base.show(obj::RBE)
    obj.model.contrasts
    for f in obj.factors
    end
end
=#
function rbe(df; var::Symbol, subject::Symbol, formulation::Symbol, period::Symbol, sequence::Symbol)
    categorical!(df, subject);
    categorical!(df, formulation);
    categorical!(df, period);
    categorical!(df, sequence);
    sort!(df, [subject, formulation, period])
    Xf = @eval(@formula($var ~ $sequence + $period + $formulation))
    Zf = @eval(@formula($var ~ 0 + $formulation))
    MF = ModelFrame(Xf, df)
    X  = ModelMatrix(MF).m
    Z  = ModelMatrix(ModelFrame(Zf, df, contrasts = Dict(formulation => StatsModels.FullDummyCoding()))).m

    if size(Z)[2] != 2 error("Size random effect matrix != 2. Not implemented yet!") end

    y = df[:, var]
    Xv, Zv, yv = sortsubjects(df, subject, X, Z, y)
    checkdata(X, Z, Xv, Zv, y)

    θvec0 = [0.02, 0.02, 0.2, 0.2, 0.001]

    remlf(x) = -2*reml(yv, Zv, X, Xv, x)

    O  = optimize(remlf, θvec0, method=Newton(),  g_tol= 1e-12, x_tol = 1e-12, store_trace = true, extended_trace = true, show_trace = false)
    θ  = Optim.minimizer(O)
    H  = Optim.trace(O)[end].metadata["h(x)"]
    dH = det(H)
    H[5,:] .= 0
    H[:,5] .= 0
    remlv = remlf(θ)
    G     = gmat(θ[3], θ[4], θ[5])
    Rv    = rmatvec(θ[1], θ[2], Zv)
    Vv    = vmatvec(Zv, G, Rv)
    iVv   = inv.(Vv)

    β     = βcoef(yv, X, Xv, iVv)
    #=
    L    = [0 0 0 0 0 1]
    C    = cmat(Xv, Zv, θ)
    se   = sqrt((L*C*L')[1])
    F    = (L*β)'*inv(L*C*L')*(L*β)
    lclg(x) = lcgf(L, Xv, Zv, x)
    g       = ForwardDiff.gradient(lclg, θ)
    df      = 2*((L*C*L')[1])^2/(g'*(2*pinv(H))*g)
    =#
    A         = 2*pinv(H)
    se, F, df = ctrst(X, Xv, Zv, θ, β, A)
    return RBE(MF, [sequence, period, formulation], β, θ, remlv, se, F, df, Rv, Vv, G, A, H, dH)
end

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

function gmat(σ₁, σ₂, ρ)
    if ρ > 1.0 ρ = 1.0 end
    if σ₁ < 0.0 σ₁ = 1.0e-6 end
    if σ₂ < 0.0 σ₂ = 1.0e-6 end
    cov = sqrt(σ₁ * σ₂) * ρ
    return [σ₁ cov; cov σ₂]
end
#=
function rmat(σ, Z)
    if σ[1] < 0.0 σ[1] = 1.0e-6 end
    if σ[2] < 0.0 σ[2] = 1.0e-6 end
    return Matrix(Diagonal((Z*σ)[:,1]))
end
=#
function rmat(σ₁, σ₂, Z)
    if σ₁ < 0.0 σ₁ = 1.0e-6 end
    if σ₂ < 0.0 σ₂ = 1.0e-6 end
    M = zeros(4, 4)
    for i = 1:4
        M[i, i] = Z[i,:]'*[σ₁, σ₂]
    end
    return M
end
function rmat2(σ, Z)
    if σ[1] < 0.0 σ[1] = 1.0e-6 end
    if σ[2] < 0.0 σ[2] = 1.0e-6 end
    return Matrix(Diagonal((Z*σ)[:,1]))
end
function cov(G, R, Z)
    V  = Z*G*Z' + R
end
function rmatvec(σ₁, σ₂, Zvec)
    n  = length(Zvec)
    Ra = Array{Array{Float64,2}, 1}(undef, n)
    for i = 1:n
        Ra[i] = rmat(σ₁, σ₂, Zvec[i])
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
function lcgf(L, Xv, Zv, θ)
    p   = size(Xv[1])[2]
    C   = zeros(p,p)
    G   = gmat(θ[3], θ[4], θ[5])
    for i=1:length(Xv)
        R   = rmat2([θ[1], θ[2]], Zv[i])
        iVv = inv(cov(G, R, Zv[i]))
        C   = C + Xv[i]'*iVv*Xv[i]
    end
    return (L*inv(C)*L')[1]
end
function ctrst(X, Xv, Zv, θ, β, A)
    p     = rank(X)
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

end # module
