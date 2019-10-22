#-------------------------------------------------------------------------------
#                               GENERAL FUNCTIONS
#
#-------------------------------------------------------------------------------
"""
    Make X, Z matrices and vector y for each subject;
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
@inline function gmat(σ::Vector)::Matrix
    if σ[3] > 1.0 σ[3] = 1.0 end
    if σ[3] < 0.0 σ[3] = 0.0 end
    if σ[1] < 0.0 σ[1] = 1.0e-6 end
    if σ[2] < 0.0 σ[2] = 1.0e-6 end
    cov = sqrt(σ[1] * σ[2]) * σ[3]
    return [σ[1] cov; cov σ[2]]
end
"""
    G matrix  (memory pre-allocation)
"""
@inline function gmat!(G::Matrix{Float64}, σ::Vector)
    if σ[3] > 1.0 σ[3] = 1.0 end
    if σ[3] < 0.0 σ[3] = 0.0 end
    if σ[1] < 0.0 σ[1] = 1.0e-6 end
    if σ[2] < 0.0 σ[2] = 1.0e-6 end
    G[1, 1] = σ[1]
    G[2, 2] = σ[2]
    G[1, 2] = G[2, 1] = sqrt(σ[1] * σ[2]) * σ[3]
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
"""
    R matrix  (memory pre-allocation)
"""
@inline function rmat!(R::Matrix{Float64}, σ::Array{Float64, 1}, Z::Matrix{Float64})
    if σ[1] < 0.0 σ[1] = 1.0e-6 end
    if σ[2] < 0.0 σ[2] = 1.0e-6 end
    copyto!(R, Matrix(Diagonal((Z*σ)[:,1])))
    return
end
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
"""
    Return variance-covariance matrix V
"""
@inline function vmat(G::Matrix{S}, R::Matrix{T}, Z::Matrix{U}) where S <: Real where T <: Real where U <: Real
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
function mvmat(G::Matrix{S}, σ::Vector{T}, Z::Matrix{U}, cache) where S <: Real where T <: Real where U <: Real
    h = hash(tuple(σ, Z))
    if h in keys(cache)
        return cache[h]
    else
        if σ[1] < 0.0 σ[1] = 1.0e-6 end
        if σ[2] < 0.0 σ[2] = 1.0e-6 end
        V  = Z*G*Z' + Matrix(Diagonal((Z*σ)[:,1]))
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
@inline function cmat(Xv::Vector{Matrix{T}}, Zv::Vector, iVv::Vector, θ::Vector)::Array{Float64, 2} where T <: Real
    p = size(Xv[1])[2]
    C = zeros(p,p)
    for i=1:length(Xv)
        @inbounds C .+= Xv[i]'*iVv[i]*Xv[i]
    end
    return inv(C)
end
#println("θ₁: ", θ1, " θ₂: ",  θ2,  " θ₃: ", θ3)

function minv(M, cache)
    h = hash(M)
    if h in keys(cache)
        return cache[h]
    else
        iM = inv(M)
        cache[h] = iM
        return iM
    end
end
function mlogdet(M, cache)
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
    mXviVXv   = zeros(promote_type(Float64, eltype(θvec)), p, p)
    for i = 1:maxobs
        mXviV[i] =  zeros(promote_type(Float64, eltype(θvec)), p, i)
    end
    cache     = Dict()
    cachel    = Dict()
    cachem    = Dict()
    #---------------------------------------------------------------------------
    n         = length(yv)
    N         = sum(length.(yv))
    G         = gmat(θvec[3:5])
    c         = (N-p)*LOG2PI
    θ1        = 0
    θ2        = zeros(promote_type(Float64, eltype(θvec)), p, p)
    θ3        = 0
    iV        = nothing
    θr        = θvec[1:2]
    for i = 1:n
        if MEMOPT && memopt
            #@inbounds R    = mrmat(θr, Zv[i], cacher)
            #@inbounds V    = memvmat(memzgz(G, Zv[i]), R)
            @inbounds V    = mvmat(G, θr, Zv[i], cachem)
            iV   = minv(V, cache)
            θ1  += mlogdet(V, cachel)
        else
            @inbounds V    = vmat(G, rmat(θr, Zv[i]), Zv[i])
            iV   = inv(V)
            θ1  += logdet(V)
        end

        #-----------------------------------------------------------------------
        #θ2 += Xv[i]'*iV*Xv[i]
        @inbounds mul!(mXviV[size(Xv[i])[1]], Xv[i]', iV)
        @inbounds mul!(mXviVXv, mXviV[size(Xv[i])[1]], Xv[i])
        θ2  += mXviVXv
        #-----------------------------------------------------------------------
        @inbounds r    = yv[i]-Xv[i]*β
        θ3  += r'*iV*r
    end
    return   -(θ1 + logdet(θ2) + θ3 + c)/2
end
"""
    REML estimation with β recalculation
"""
function remlb(yv, Zv, p, Xv, θvec, β; memopt::Bool = true)
    maxobs    = maximum(length.(yv))
    #some memory optimizations to reduse allocations
    mXviV     = Array{Array{eltype(θvec), 2}, 1}(undef, maxobs)
    mXviVXv   = zeros(promote_type(Float64, eltype(θvec)), p, p)
    for i = 1:maxobs
        mXviV[i] =  zeros(promote_type(Float64, eltype(θvec)), p, i)
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
    θ2        = zeros(promote_type(Float64, eltype(θvec)), p, p)
    θ3        = 0
    iV        = nothing
    βm        = zeros(promote_type(Float64, eltype(θvec)), p)
    βt        = zeros(promote_type(Float64, eltype(θvec)), p)
    θr        = θvec[1:2]
    for i = 1:n
        if MEMOPT && memopt
            #@inbounds R        = memrmat(θr, Zv[i])
            #@inbounds V        = memvmat(memzgz(G, Zv[i]), R)
            @inbounds V        = mvmat(G, θr, Zv[i], cachem)
            @inbounds iVv[i]   = minv(V, cache)
            θ1                += mlogdet(V, cachel)
        else
            @inbounds R        = rmat(θr, Zv[i])
            @inbounds V        = vmat(G, R, Zv[i])
            @inbounds iVv[i]   = inv(V)
            θ1                += logdet(V)
        end

        #-----------------------------------------------------------------------
        #θ2 += Xv[i]'*iV*Xv[i]
        @inbounds mul!(mXviV[size(Xv[i])[1]], Xv[i]', iVv[i])
        @inbounds mul!(mXviVXv, mXviV[size(Xv[i])[1]], Xv[i])
        θ2  += mXviVXv
        @inbounds βm  .+= mXviV[size(Xv[i])[1]]*yv[i]
        #-----------------------------------------------------------------------
        #tm   = Xv[i]'*iVv[i]    #Temp matrix for Xv[i]'*iV*Xv[i] and Xv[i]'*iV*yv[i] calc
        #θ2m .+= tm*Xv[i]
        #βm  .+= tm*yv[i]
    end
    mul!(βt, inv(θ2), βm)
    for i = 1:n
        @inbounds r    = yv[i] - Xv[i]*βt
        @inbounds θ3  += r'*iVv[i]*r
    end
    return   -(θ1 + logdet(θ2) + θ3 + c)/2
end
"""
Satterthwaite DF gradient function.
"""
function lclgf(L, Lt, Xv, Zv, θ; memopt::Bool = true)
    p     = size(Xv[1])[2]
    G     = gmat(θ[3:5])
    C     = zeros(promote_type(Float64, eltype(θ)), p, p)
    θr    = θ[1:2]
    cache     = Dict()
    cachem    = Dict()
    for i = 1:length(Xv)
        if MEMOPT && memopt
            iV   = minv(mvmat(G, θr, Zv[i], cachem), cache)
        else
            R   = rmat(θr, Zv[i])
            iV  = inv(vmat(G, R, Zv[i]))
        end
        C  += Xv[i]'*iV*Xv[i]
    end
    return (L*inv(C)*Lt)[1]
end
#-------------------------------------------------------------------------------
#             REML FOR OPT ALGORITHM
#-------------------------------------------------------------------------------
"""
    REML with β final update
"""
function reml2b!(yv::S, Zv::T, p::Int, n::Int, N::Int, Xv::T, G::Array{Float64, 2}, Rv::T, Vv::T, iVv::T, θvec::Array{Float64, 1}, β::Array{Float64, 1}, mem::MemAlloc)::Float64 where T <: Array{Array{Float64, 2}, 1} where S <: Array{Array{Float64, 1}, 1}
    gmat!(G, θvec[3:5])
    c  = (N-p)*LOG2PI #log(2π)
    θ1 = 0
    θ2  = zeros(p, p)
    θ3 = 0
    iV   = nothing
    #fill!(mem.mem4, 0)
    βm   = zeros(p)
    θr   = [θvec[1], θvec[2]]
    cache     = Dict()
    cachel    = Dict()
    @inbounds for i = 1:n
        rmat!(Rv[i], θr, Zv[i])
        vmat!(Vv[i], G, Rv[i], Zv[i], mem.mem1)
        #Memopt!
        copyto!(iVv[i], minv(Vv[i], cache))
        θ1  += mlogdet(Vv[i], cachel)
        #-
        mul!(mem.mem2[size(Xv[i])[1]], Xv[i]', iVv[i])
        θ2    .+= mem.mem2[size(Xv[i])[1]]*Xv[i]
        βm    .+= mem.mem2[size(Xv[i])[1]]*yv[i]
    end
    mul!(β, inv(θ2), βm)
    for i = 1:n
        copyto!(mem.mem3[length(yv[i])], yv[i])
        mem.mem3[length(yv[i])] .-= Xv[i]*β
        θ3  += mem.mem3[length(yv[i])]'*iVv[i]*mem.mem3[length(yv[i])]
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
        if length(fvd) > 1 push!(fv, var(fvd)) end
    end
    sv = Array{Float64, 1}(undef, 0)
    for i in u
        fvd = Array{Float64, 1}(undef, 0)
        v   = findall(x->x==i, df[:, sbj])
        for r in v
            push!(fvd, df[r, dv])
        end
        if length(fvd) > 1 push!(sv, var(fvd)) end
    end
    push!(fv, mean(filter!(x -> x !== NaN, sv)))
    return fv
end
