# 6.582 6.568 6.539
#-------------------------------------------------------------------------------
#                               GENERAL FUNCTIONS
#
#-------------------------------------------------------------------------------
"""
    Make X, Z matrices and vector y for each subject;
"""
function sortsubjects(df::DataFrame, sbj::Symbol, X::Matrix, Z::Matrix, y::Vector)
    u = unique(df[:, sbj])
    Xa = Vector{Matrix{eltype(X)}}(undef, length(u))
    Za = Vector{Matrix{eltype(Z)}}(undef, length(u))
    ya = Vector{Vector{eltype(y)}}(undef, length(u))
    for i = 1:length(u)
        v = findall(x->x==u[i], df[:, sbj])
        Xs = Vector{eltype(X)}(undef, 0)
        Zs = Vector{eltype(Z)}(undef, 0)
        ys = Vector{eltype(y)}(undef, 0)
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
@inline function gmat(σ::Vector)::AbstractMatrix
    #m = Matrix(Diagonal(σ[1:2]))
    #m[1, 2] = m[2, 1] = sqrt(σ[1] * σ[2]) * σ[3]
    #return m
    cov = sqrt(σ[1] * σ[2]) * σ[3]
    return Symmetric([σ[1] cov; cov σ[2]])
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
    return Diagonal(Z*σ)
end
"""
    R matrix  (memory pre-allocation)
"""
@inline function rmat!(R::AbstractMatrix{T}, σ::Vector{T}, Z::Matrix{T}) where T <: AbstractFloat
    copyto!(R, Diagonal(Z*σ))
    return
end
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
"""
    Return variance-covariance matrix V
"""
@inline function vmat(G::AbstractMatrix, R::AbstractMatrix, Z::Matrix)::AbstractMatrix
    return  mulαβαtc(Z, G, R)
end
@inline function vmat!(V::Matrix{T}, G::AbstractMatrix{T}, R::AbstractMatrix{T}, Z::Matrix{T}, memc) where T <: AbstractFloat
    #copyto!(V, Z*G*Z')
    mul!(memc[size(Z)[1]], Z, G)
    mul!(V, memc[size(Z)[1]], Z')
    V .+= R
    return
end
function mvmat(G::AbstractMatrix, σ::Vector, Z::Matrix, cache)::Matrix
    h = hash(tuple(σ, Z))
    if h in keys(cache)
        return cache[h]
    else
        #V  = mulαβαtc(Z, G, Diagonal(Z*σ), mem)
        V   = Z * G * Z' + Diagonal(Z*σ)
        cache[h] = V
        return V
    end
end
function mvmat(G::AbstractMatrix, σ::Vector, Z::Matrix, mem, cache)::Matrix
    h = hash(tuple(σ, Z))
    if h in keys(cache)
        return cache[h]
    else
        V  = mulαβαtc(Z, G, Diagonal(Z*σ), mem)
        #V   = Z * G * Z' + Diagonal(Z*σ)
        cache[h] = V
        return V
    end
end
function mvmatall(G::AbstractMatrix, σ::Vector, Z::Matrix, mem, cache)
    #h = hash(tuple(σ, Z))
    h = hash(Z)
    if h in keys(cache)
        return cache[h]
    else
        V   = mulαβαtc(Z, G, Diagonal(Z*σ), mem)
        #V   = Z * G * Z' + Diagonal(Z*σ)
        iV  = inv(V)
        ldV = logdet(V)
        cache[h] = (V, iV, ldV)
        return V, iV, ldV
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
#-------------------------------------------------------------------------------
#             REML FOR OPT ALGORITHM
#-------------------------------------------------------------------------------
"""
    -2 REML with β final update
"""
function reml2b!(data::RBEDataStructure, G::Matrix{T}, Rv::Vector, Vv::Vector, iVv::Vector,
        θvec::Vector{T}, β::Vector{T}, mem::MemAlloc) where T <: AbstractFloat

    rebuildcache(data, promote_type(eltype(first(data.yv)), eltype(θvec)))
    gmat!(G, θvec[3:5])

    θ1 = 0
    θ2 = zeros(data.p, data.p)
    θ3 = 0
    #iV   = nothing
    #fill!(mem.mem4, 0)
    βm   = zeros(data.p)
    cache     = Dict()
    cachel    = Dict()
    @inbounds for i = 1:data.n
        rmat!(Rv[i], θvec[1:2], data.Zv[i])
        #vmat!(Vv[i], G, Rv[i], Zv[i], mem.mem1)

        #Memopt!
        Vv[i], iVv[i], ldV         = mvmatall(G, θvec[1:2], data.Zv[i], first(data.mem.svec), cache)

        #copyto!(iVv[i], minv(Vv[i], cache))
        #θ1  += mlogdet(Vv[i], cachel)
        θ1  += ldV
        #-
        #θ2 += Xv[i]'*iVv[i]*Xv[i]
        mulall!(θ2, βm, data.Xv[i], iVv[i], data.yv[i], first(data.mem.svec))
    end
    mul!(β, inv(θ2), βm)
    for i = 1:data.n

        #copyto!(mem.mem3[length(yv[i])], yv[i])
        #mem.mem3[length(yv[i])] .-= Xv[i] * β
        #θ3  += mem.mem3[length(yv[i])]' * iVv[i] * mem.mem3[length(yv[i])]

        #Same:
        #r    = yv[i] - Xv[i]*β
        #θ3  += r'*iVv[i]*r

        @inbounds θ3  += mulall(data.yv[i], data.Xv[i], β, iVv[i], first(data.mem.svec))
    end
    return   θ1 + logdet(θ2) + θ3 + data.remlc
end
"""
    -2 REML function for ForwardDiff
"""
function reml2(data::RBEDataStructure, θvec::Vector, β::Vector; memopt::Bool = true)

    #memory optimizations to reduse allocations (cache rebuild)
    rebuildcache(data, promote_type(eltype(data.yv[1]), eltype(θvec)))
    cache     = Dict()
    #cachel    = Dict()
    #cachem    = Dict()
    #---------------------------------------------------------------------------
    G         = gmat(θvec[3:5])
    θ1        = 0
    θ2        = zeros(promote_type(eltype(data.yv[1]), eltype(θvec)), data.p, data.p)
    θ3        = 0
    iV        = nothing
    for i = 1:data.n
        if MEMOPT && memopt

            ##@inbounds V    = mvmat(G, θvec[1:2], data.Zv[i], first(data.mem.svec), cachem)
            ##iV             = minv(V, cache)
            ##θ1            += mlogdet(V, cachel)

            V, iV, ldV         = mvmatall(G, θvec[1:2], data.Zv[i], first(data.mem.svec), cache)
            θ1                += ldV
        else
            @inbounds V    = vmat(G, rmat(θvec[1:2], data.Zv[i]), data.Zv[i])
            iV             = inv(V)
            θ1            += logdet(V)
        end
        #-----------------------------------------------------------------------
        #θ2 += Xv[i]'*iV*Xv[i]
        @inbounds mulall!(θ2, data.Xv[i], iV, first(data.mem.svec))
        #-----------------------------------------------------------------------
        #r    = yv[i] - Xv[i] * β
        #θ3  += r' * iV * r
        @inbounds θ3  += mulall(data.yv[i], data.Xv[i], β, iV, first(data.mem.svec))
    end
    return   θ1 + logdet(θ2) + θ3 + data.remlc
end
"""
    -2 REML estimation with β recalculation for ForwardDiff
"""
function reml2bfd(data::RBEDataStructure, θvec::Vector; memopt::Bool = true)
    return reml2b(data, θvec; memopt = memopt)[1]
end
function reml2b(data::RBEDataStructure, θvec::Vector; memopt::Bool = true)

    rebuildcache(data, promote_type(eltype(data.yv[1]), eltype(θvec)))
    #maxobs    = data.maxobs
    #mem       = MemCache(eltype(θvec), maxobs)
    #some memory optimizations to reduse allocations
    #mXviV     = Vector{Matrix{eltype(θvec)}}(undef, maxobs)
    #mc        = Vector{Vector{eltype(θvec)}}(undef, data.maxobs)
    #mXviVXv   = zeros(promote_type(eltype(yv[1]), eltype(θvec)), p, p)
    #for i = 1:data.maxobs
        #mXviV[i] =  zeros(promote_type(eltype(yv[1]), eltype(θvec)), p, i)
        #mc[i]    =  zeros(promote_type(eltype(data.yv[1]), eltype(θvec)), i)
    #end
    cache     = Dict()
    #cachel    = Dict()
    #cachem    = Dict()
    #---------------------------------------------------------------------------
    #n         = length(yv)
    #N         = sum(length.(yv))
    G         = gmat(θvec[3:5])
    iVv       = Array{Array{eltype(θvec), 2}, 1}(undef, data.n)
    V         = nothing
    ldV       = nothing
    #c         = data.remlc
    θ1        = 0
    θ2        = zeros(promote_type(eltype(first(data.yv)), eltype(θvec)), data.p, data.p)
    θ3        = 0
    #iV        = nothing
    βm        = zeros(promote_type(eltype(first(data.yv)), eltype(θvec)), data.p)
    β         = zeros(promote_type(eltype(first(data.yv)), eltype(θvec)), data.p)

    for i = 1:data.n
        if MEMOPT && memopt
            #@inbounds R        = memrmat(θr, Zv[i])
            #@inbounds V        = memvmat(memzgz(G, Zv[i]), R)
            ##@inbounds V        = mvmat(G, θvec[1:2], data.Zv[i], cachem)
            ##@inbounds iVv[i]   = minv(V, cache)
            ##θ1                += mlogdet(V, cachel)

            @inbounds V, iVv[i], ldV = mvmatall(G, θvec[1:2], data.Zv[i], first(data.mem.svec), cache)
            
            θ1                += ldV
        else
            @inbounds R        = rmat(θvec[1:2], data.Zv[i])
            @inbounds V        = vmat(G, R, data.Zv[i])
            @inbounds iVv[i]   = inv(V)
            θ1                += logdet(V)
        end

        #-----------------------------------------------------------------------
        #θ2 += Xv[i]'*iVv[i]*Xv[i]
        #-
        #@inbounds mul!(mXviV[size(Xv[i])[1]], Xv[i]', iVv[i])
        #@inbounds mul!(mXviVXv, mXviV[size(Xv[i])[1]], Xv[i])
        #θ2  += mXviVXv
        #@inbounds βm  .+= mXviV[size(Xv[i])[1]] * yv[i]
        mulall!(θ2, βm, data.Xv[i], iVv[i], data.yv[i], first(data.mem.svec))
        #-----------------------------------------------------------------------
        #tm   = Xv[i]'*iVv[i]    #Temp matrix for Xv[i]'*iV*Xv[i] and Xv[i]'*iV*yv[i] calc
        #θ2m .+= tm*Xv[i]
        #βm  .+= tm*yv[i]
    end
    mul!(β, inv(θ2), βm)
    for i = 1:data.n
        #@inbounds r    = yv[i] - Xv[i] * βt
        #@inbounds θ3  += r' * iVv[i] * r

        @inbounds θ3  += mulall(data.yv[i], data.Xv[i], β, iVv[i], first(data.mem.svec))
    end

    return   θ1 + logdet(θ2) + θ3 + data.remlc,  β, θ2
end
#-------------------------------------------------------------------------------
"""
    Return C matrix
    var(β) p×p variance-covariance matrix
"""
@inline function cmat(Xv::Vector{Matrix{T}}, Zv::Vector, iVv::Vector, θ::Vector)::Matrix where T <: AbstractFloat
    p = size(Xv[1])[2]
    C = zeros(p, p)
    for i=1:length(Xv)
        #@inbounds C .+= Xv[i]' * iVv[i] * Xv[i]
        mulall!(C, Xv[i], iVv[i])
    end
    return pinv(C)
end
#-------------------------------------------------------------------------------
# C matrix derivation
#-------------------------------------------------------------------------------
"""
non inverted C matrix gradient function
"""
function cmatgf(Xv::Vector, Zv::Vector, θ::Vector; memopt::Bool = true)
    p      = size(Xv[1], 2)
    jC     = ForwardDiff.jacobian(x -> cmatvec(Xv, Zv, x; memopt = memopt), θ)
    result = Vector{Matrix}(undef, 0)
    for i in 1:length(θ)
        push!(result, reshape(jC[:,i], p, p))
    end
    return result
end
"""
non inverted C matrix in vector form for gradient
"""
function cmatvec(Xv::Vector, Zv::Vector, θ::Vector; memopt::Bool = true)
    p     = size(Xv[1], 2)
    G     = gmat(θ[3:5])
    #mem   = MemCache(eltype(θ), 4)
    C     = zeros(promote_type(eltype(Zv[1]), eltype(θ)), p, p)
    cache     = Dict()
    cachem    = Dict()
    for i = 1:length(Xv)
        if memopt
            iV   = minv(mvmat(G, θ[1:2], Zv[i], cachem), cache)
        else
            R   = rmat(θ[1:2], Zv[i])
            iV  = inv(vmat(G, R, Zv[i]))
        end
        #C  += Xv[i]' * iV * Xv[i]
        mulall!(C, Xv[i], iV)
    end
    return C[:]
end
"""
C matrix gradients
"""
function cmatg(Xv::Vector, Zv::Vector, θ::Vector, C::Matrix; memopt::Bool = true)
    g  = Vector{Matrix}(undef, length(θ))
    jC = cmatgf(Xv, Zv, θ; memopt = memopt)
    for i = 1:length(θ)
        g[i] = (- C * jC[i] * C)
    end
    return g
end
"""
L * C * L' for all C marices
"""
function lclg(gradc, L)
    g  = Vector{eltype(gradc[1])}(undef, length(gradc))
    for i = 1:length(gradc)
        g[i] = (L * gradc[i] * L')[1]
    end
    return g
end
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
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
function rholinksigmoid2r(ρ, m)
    return tan(ρ)
end

function varlinkmap(θ, r1, r2, f1, f2)
    θl      = similar(θ)
    θl[r1]  = f1.(θ[r1])
    θl[r2]  = f2.(θ[r2])
    return θl
end
