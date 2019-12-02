#-------------------------------------------------------------------------------
#                               GENERAL FUNCTIONS
#
#-------------------------------------------------------------------------------
"""
    Make X, Z matrices and vector y for each subject;
"""
function sortsubjects(df::DataFrame, sbj::Symbol, X::Matrix, Z::Matrix, y::Vector)
    u = unique(df[:, sbj])
    Xa = Array{Array{eltype(X),2}, 1}(undef, length(u))
    Za = Array{Array{eltype(Z),2}, 1}(undef, length(u))
    ya = Array{Array{eltype(y),1}, 1}(undef, length(u))
    for i = 1:length(u)
        v = findall(x->x==u[i], df[:, sbj])
        Xs = Array{eltype(X), 1}(undef, 0)
        Zs = Array{eltype(Z), 1}(undef, 0)
        ys = Array{eltype(y), 1}(undef, 0)
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
@inline function gmat(σ::Vector)::Matrix
    #m = Matrix(Diagonal(σ[1:2]))
    #m[1, 2] = m[2, 1] = sqrt(σ[1] * σ[2]) * σ[3]
    #return m
    cov = sqrt(σ[1] * σ[2]) * σ[3]
    return [σ[1] cov; cov σ[2]]
end
"""
    G matrix  (memory pre-allocation)
"""
@inline function gmat!(G::Matrix{T}, σ::Vector) where T <: AbstractFloat
    G[1, 1] = σ[1]
    G[2, 2] = σ[2]
    G[1, 2] = G[2, 1] = sqrt(σ[1] * σ[2]) * σ[3]
    return
end

"""
    R matrix (ForwardDiff+)
"""
@inline function rmat(σ::Vector, Z::Matrix)::Matrix
    return Matrix(Diagonal((Z*σ)))
end
"""
    R matrix  (memory pre-allocation)
"""
@inline function rmat!(R::Matrix{T}, σ::Vector{T}, Z::Matrix{T}) where T <: AbstractFloat
    copyto!(R, Matrix(Diagonal((Z*σ))))
    return
end
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
"""
    Return variance-covariance matrix V
"""
@inline function vmat(G::Matrix, R::Matrix, Z::Matrix)::Matrix
    V  = Z * G * Z' + R
    return V
end
@inline function vmat!(V::Matrix{T}, G::Matrix{T}, R::Matrix{T}, Z::Matrix{T}, memc) where T <: AbstractFloat
    #copyto!(V, Z*G*Z')
    mul!(memc[size(Z)[1]], Z, G)
    mul!(V, memc[size(Z)[1]], Z')
    V .+= R
    return
end
function mvmat(G::Matrix, σ::Vector, Z::Matrix, cache)::Matrix
    h = hash(tuple(σ, Z))
    if h in keys(cache)
        return cache[h]
    else
        V  = Z * G * Z' + Matrix(Diagonal((Z*σ)))
        cache[h] = V
        return V
    end
end
"""
    Feel M with zero feeled matrices
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
@inline function cmat(Xv::Vector{Matrix{T}}, Zv::Vector, iVv::Vector, θ::Vector)::Matrix where T <: AbstractFloat
    p = size(Xv[1])[2]
    C = zeros(p, p)
    for i=1:length(Xv)
        @inbounds C .+= Xv[i]' * iVv[i] * Xv[i]
    end
    return inv(C)
end
#println("θ₁: ", θ1, " θ₂: ",  θ2,  " θ₃: ", θ3)

function minv(M::Matrix, cache::Dict)::Matrix
    h = hash(M)
    if h in keys(cache)
        return cache[h]
    else
        iM = inv(M)
        cache[h] = iM
        return iM
    end
end
function mlogdet(M::Matrix, cache::Dict)
    h = hash(M)
    if h in keys(cache)
        return cache[h]
    else
        iM = logdet(M)
        cache[h] = iM
        return iM
    end
end

"""
    REML function for ForwardDiff
"""
function reml(yv::Vector, Zv::Vector, p::Int, Xv::Vector, θvec::Vector, β::Vector; memopt::Bool = true)
    maxobs    = maximum(length.(yv))
    #some memory optimizations to reduse allocations
    mXviV     = Array{Array{eltype(θvec), 2}, 1}(undef, maxobs)
    mXviVXv   = zeros(promote_type(eltype(yv[1]), eltype(θvec)), p, p)
    for i = 1:maxobs
        mXviV[i] =  zeros(promote_type(eltype(yv[1]), eltype(θvec)), p, i)
    end
    cache     = Dict()
    cachel    = Dict()
    cachem    = Dict()
    #---------------------------------------------------------------------------
    n         = length(yv)
    N         = sum(length.(yv))
    G         = gmat(θvec[3:5])
    c         = (N - p) * LOG2PI
    θ1        = 0
    θ2        = zeros(promote_type(eltype(yv[1]), eltype(θvec)), p, p)
    θ3        = 0
    iV        = nothing
    for i = 1:n
        if MEMOPT && memopt
            #@inbounds R    = mrmat(θr, Zv[i], cacher)
            #@inbounds V    = memvmat(memzgz(G, Zv[i]), R)
            @inbounds V    = mvmat(G, θvec[1:2], Zv[i], cachem)
            iV   = minv(V, cache)
            θ1  += mlogdet(V, cachel)
        else
            @inbounds V    = vmat(G, rmat(θvec[1:2], Zv[i]), Zv[i])
            iV   = inv(V)
            θ1  += logdet(V)
        end

        #-----------------------------------------------------------------------
        #θ2 += Xv[i]'*iV*Xv[i]
        @inbounds mul!(mXviV[size(Xv[i])[1]], Xv[i]', iV)
        @inbounds mul!(mXviVXv, mXviV[size(Xv[i])[1]], Xv[i])
        θ2  += mXviVXv
        #-----------------------------------------------------------------------
        @inbounds r    = yv[i] - Xv[i] * β
        θ3  += r' * iV * r
    end
    return   -(θ1 + logdet(θ2) + θ3 + c)/2
end
"""
    REML estimation with β recalculation
"""
function remlb(yv::Vector, Zv::Vector, p::Int, Xv::Vector, θvec::Vector, β::Vector; memopt::Bool = true)
    maxobs    = maximum(length.(yv))
    #some memory optimizations to reduse allocations
    mXviV     = Array{Array{eltype(θvec), 2}, 1}(undef, maxobs)
    mXviVXv   = zeros(promote_type(eltype(yv[1]), eltype(θvec)), p, p)
    for i = 1:maxobs
        mXviV[i] =  zeros(promote_type(eltype(yv[1]), eltype(θvec)), p, i)
    end
    cache     = Dict()
    cachel    = Dict()
    cachem    = Dict()
    #---------------------------------------------------------------------------
    n         = length(yv)
    N         = sum(length.(yv))
    G         = gmat(θvec[3:5])
    iVv       = Array{Array{eltype(θvec), 2}, 1}(undef, n)
    c         = (N-p)*LOG2PI
    θ1        = 0
    θ2        = zeros(promote_type(eltype(yv[1]), eltype(θvec)), p, p)
    θ3        = 0
    iV        = nothing
    βm        = zeros(promote_type(eltype(yv[1]), eltype(θvec)), p)
    βt        = zeros(promote_type(eltype(yv[1]), eltype(θvec)), p)
    for i = 1:n
        if MEMOPT && memopt
            #@inbounds R        = memrmat(θr, Zv[i])
            #@inbounds V        = memvmat(memzgz(G, Zv[i]), R)
            @inbounds V        = mvmat(G, θvec[1:2], Zv[i], cachem)
            @inbounds iVv[i]   = minv(V, cache)
            θ1                += mlogdet(V, cachel)
        else
            @inbounds R        = rmat(θvec[1:2], Zv[i])
            @inbounds V        = vmat(G, R, Zv[i])
            @inbounds iVv[i]   = inv(V)
            θ1                += logdet(V)
        end

        #-----------------------------------------------------------------------
        #θ2 += Xv[i]'*iV*Xv[i]
        @inbounds mul!(mXviV[size(Xv[i])[1]], Xv[i]', iVv[i])
        @inbounds mul!(mXviVXv, mXviV[size(Xv[i])[1]], Xv[i])
        θ2  += mXviVXv
        @inbounds βm  .+= mXviV[size(Xv[i])[1]] * yv[i]
        #-----------------------------------------------------------------------
        #tm   = Xv[i]'*iVv[i]    #Temp matrix for Xv[i]'*iV*Xv[i] and Xv[i]'*iV*yv[i] calc
        #θ2m .+= tm*Xv[i]
        #βm  .+= tm*yv[i]
    end
    mul!(βt, inv(θ2), βm)
    for i = 1:n
        @inbounds r    = yv[i] - Xv[i] * βt
        @inbounds θ3  += r' * iVv[i] * r
    end

    return   -(θ1 + logdet(θ2) + θ3 + c)/2
end
"""
Satterthwaite DF gradient function.
"""
function lclgf(L, Lt, Xv::Vector, Zv::Vector, θ::Vector; memopt::Bool = true)
    p     = size(Xv[1])[2]
    G     = gmat(θ[3:5])
    C     = zeros(promote_type(eltype(Zv[1]), eltype(θ)), p, p)
    cache     = Dict()
    cachem    = Dict()
    for i = 1:length(Xv)
        if MEMOPT && memopt
            iV   = minv(mvmat(G, θ[1:2], Zv[i], cachem), cache)
        else
            R   = rmat(θ[1:2], Zv[i])
            iV  = inv(vmat(G, R, Zv[i]))
        end
        C  += Xv[i]' * iV * Xv[i]
    end
    return (L * inv(C) * Lt)[1]
end
#-------------------------------------------------------------------------------
#             REML FOR OPT ALGORITHM
#-------------------------------------------------------------------------------
"""
    REML with β final update
"""
function reml2b!(yv::Vector, Zv::Vector, p::Int, n::Int, N::Int,
        Xv::Vector, G::Matrix{T}, Rv::Vector, Vv::Vector, iVv::Vector,
        θvec::Vector{T}, β::Vector{T}, mem::MemAlloc) where T <: AbstractFloat
    gmat!(G, θvec[3:5])
    c  = (N-p)*LOG2PI #log(2π)
    θ1 = 0
    θ2  = zeros(p, p)
    θ3 = 0
    iV   = nothing
    #fill!(mem.mem4, 0)
    βm   = zeros(p)
    cache     = Dict()
    cachel    = Dict()
    @inbounds for i = 1:n
        rmat!(Rv[i], θvec[1:2], Zv[i])
        vmat!(Vv[i], G, Rv[i], Zv[i], mem.mem1)
        #Memopt!
        copyto!(iVv[i], minv(Vv[i], cache))
        θ1  += mlogdet(Vv[i], cachel)
        #-
        mul!(mem.mem2[size(Xv[i])[1]], Xv[i]', iVv[i])
        θ2    .+= mem.mem2[size(Xv[i])[1]] * Xv[i]
        βm    .+= mem.mem2[size(Xv[i])[1]] * yv[i]
    end
    mul!(β, inv(θ2), βm)
    for i = 1:n
        copyto!(mem.mem3[length(yv[i])], yv[i])
        mem.mem3[length(yv[i])] .-= Xv[i] * β
        θ3  += mem.mem3[length(yv[i])]' * iVv[i] * mem.mem3[length(yv[i])]
        #Same:
        #r    = yv[i] - Xv[i]*β
        #θ3  += r'*iVv[i]*r
    end
    return   -(θ1 + logdet(θ2) + θ3 + c)
end
#-------------------------------------------------------------------------------
"""
    Initial variance computation
"""
function initvar(df::DataFrame, dv::Symbol, fac::Symbol, sbj::Symbol)::Vector
    u  = unique(df[:, sbj])
    f  = unique(df[:, fac])

    fv = Array{eltype(df[!, dv]), 1}(undef, 0)
    sb = Array{eltype(df[!, dv]), 1}(undef, 0)
    for i in f
        push!(fv, var(df[df[!, fac] .== i, dv]))
    end
    for i in u
        sv = var(df[df[!, sbj] .== i, dv])
        if sv > 0 push!(sb, sv) end
    end
    push!(fv, mean(sb))
    return fv
end
#-------------------------------------------------------------------------------
function optimcallback(x)
    false
end
#-------------------------------------------------------------------------------
function vlink(σ)
    exp(σ)
end
function vlinkr(σ)
    log(σ)
end

function rholinkpsigmoid(ρ, m)
    return 1/(1 + exp(ρ * m))
end
function rholinkpsigmoidr(ρ, m)
    return log(1/ρ - 1)/m
end

function rholinksigmoid(ρ, m)
    return ρ/sqrt(1 + ρ^2)
end
function rholinksigmoidr(ρ, m)
    return sign(ρ)*sqrt(ρ^2/(1 - ρ^2))
end

function rholinksigmoid2(ρ, m)
    return atan(ρ)
end
function rholinksigmoidr2(ρ, m)
    return tan(ρ)
end

function varlinkmap(θ, r1, r2, f1, f2)
    θl      = similar(θ)
    θl[r1]  = f1.(θ[r1])
    θl[r2]  = f2.(θ[r2])
    return θl
end
