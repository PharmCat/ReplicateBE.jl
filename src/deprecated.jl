"""
    Secondary param estimation:
    SE
    F
    DF
    C
"""
function ctrst(p, Xv, Zv, iVv, θ, β, A, C; memopt::Bool = true)
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
        lclg(x) = lclgf(L, Lt, Xv, Zv, x; memopt = memopt)
        g       = ForwardDiff.gradient(lclg, θ)
        df[i]   = 2*((lcl)[1])^2/(g'*(A)*g)
        #LinearAlgebra.eigen(L*C*L')
    end
    return se, F, df
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
"""
    Optim reml without β
    For Pre-opt
"""
function reml2!(yv::S, Zv::T, p::Int, n::Int, N::Int, Xv::T, G::Array{Float64, 2}, Rv::T, Vv::T, iVv::T, θvec::Array{Float64, 1}, β::Array{Float64, 1}, memc, memc2, memc3, memc4)::Float64 where T <: Array{Array{Float64, 2}, 1} where S <: Array{Array{Float64, 1}, 1}
    gmat!(G, θvec[3], θvec[4], θvec[5])
    c  = (N-p)*LOG2PI
    θ1 = 0
    #θ2 = 0
    θ3 = 0
    fill!(memc4, 0)
    #θ2m  = zeros(p,p)
    θr    = [θvec[1], θvec[2]]
    for i = 1:n
        rmat!(Rv[i], θr, Zv[i])
        vmat!(Vv[i], G, Rv[i], Zv[i], memc)
        copyto!(iVv[i], inv(Vv[i]))
        θ1  += logdet(Vv[i])
        mul!(memc2[size(Xv[i])[1]], Xv[i]', iVv[i])
        memc4 .+= memc2[size(Xv[i])[1]]*Xv[i]
        #θ2m .+= Xv[i]'*iVv[i]*Xv[i]
        copyto!(memc3[length(yv[i])], yv[i])
        memc3[length(yv[i])] .-= Xv[i]*β
        θ3  += memc3[length(yv[i])]'*iVv[i]*memc3[length(yv[i])]
        #r    = yv[i]-Xv[i]*β
        #θ3  .+= r'*iVv[i]*r
    end
    #θ2       = logdet(θ2m)
    return   -(θ1 + logdet(memc4) + θ3 + c)
end

function memcalloc(p, zs, yv)
    maxobs  = maximum(length.(yv))
    memc = Array{Array{Float64, 2}, 1}(undef, maxobs)
    memc2 = Array{Array{Float64, 2}, 1}(undef, maxobs)
    memc3 = Array{Array{Float64, 1}, 1}(undef, maxobs)
    for i = 1:maxobs
        memc[i] = zeros(i, zs)
        memc2[i] = zeros(p, i)
        memc3[i] = zeros(i)
    end
    memc4 = zeros(p, p)
    return memc, memc2, memc3, memc4
end
