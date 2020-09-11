#-------------------------------------------------------------------------------
#                               GENERAL FUNCTIONS
#
#-------------------------------------------------------------------------------
"""
    Make X, Z matrices and vector y for each subject;
"""
function sortsubjects(df::DataFrame, sbj::Symbol, X::Matrix{T}, Z::Matrix{T}, y::Vector{T}) where T
    u = unique(df[!, sbj])
    Xa = Vector{SubArray{T}}(undef, length(u))
    Za = Vector{SubArray{T}}(undef, length(u))
    ya = Vector{SubArray{T}}(undef, length(u))
    @simd for i = 1:length(u)
        @inbounds v = findall(x->x==u[i], df[!, sbj])
        @inbounds Xa[i] = view(X, v, :)
        @inbounds Za[i] = view(Z, v, :)
        @inbounds ya[i] = view(y, v)
    end
    return Xa, Za, ya
end
"""
    G matrix
"""
@inline function gmat(σ::AbstractVector)::AbstractMatrix
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
@inline function rmat(σ::AbstractVector, Z::AbstractMatrix)::Matrix
    return Diagonal(Z*σ)
end
"""
    R matrix  (memory pre-allocation)
"""
@inline function rmat!(R::AbstractMatrix{T}, σ::Vector{T}, Z::AbstractMatrix{T}) where T <: AbstractFloat
    copyto!(R, Diagonal(Z*σ))
    return
end
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
"""
    Return variance-covariance matrix V
"""
@inline function vmat(G::AbstractMatrix, R::AbstractMatrix, Z::AbstractMatrix)::AbstractMatrix
    return  mulαβαtc(Z, G, R)
end
@inline function vmat!(V::Matrix{T}, G::AbstractMatrix{T}, R::AbstractMatrix{T}, Z::AbstractMatrix{T}, memc) where T <: AbstractFloat
    #copyto!(V, Z*G*Z')
    mul!(memc[size(Z)[1]], Z, G)
    mul!(V, memc[size(Z)[1]], Z')
    V .+= R
    return
end
function mvmat(G::AbstractMatrix, σ::Vector, Z::AbstractMatrix, cache)::Matrix
    #h = hash(tuple(σ, Z))
    if Z in keys(cache)
        return cache[Z]
    else
        V  = mulαβαtc(Z, G, Diagonal(Z*σ))
        #V   = Z * G * Z' + Diagonal(Z*σ)
        cache[Z] = V
        return V
    end
end
function mvmat(G::AbstractMatrix, σ::Vector, Z::AbstractMatrix, mem, cache)::Matrix
    #h = hash(tuple(σ, Z))
    if Z in keys(cache)
        return cache[Z]
    else
        V  = mulαβαtc(Z, G, Diagonal(Z*σ), mem)
        #V   = Z * G * Z' + Diagonal(Z*σ)
        cache[Z] = V
        return V
    end
end

function mvmatall(G::AbstractMatrix, σ::AbstractVector{T}, Z::AbstractMatrix, mem, cache) where T
    if Z in keys(cache)
        return cache[Z]
    else
        local V::Symmetric{T,Array{T,2}}
        local V⁻¹::Symmetric{T,Array{T,2}}
        local log│V│::T
        #V = mulαβαtcupd!(mem.mvec, Z, G, σ, mem.svec)
        V   = mulαβαtc(Z, G, σ, mem)
        #V   = Z * G * Z' + Diagonal(Z*σ)

        if size(V, 1) <= 14
            sV     = SHermitianCompact(SMatrix{size(V, 1),size(V, 1)}(V))
            V⁻¹    = Symmetric(inv(sV))
            log│V│ = logdet(sV)
        else
            V⁻¹    = Symmetric(inv(V))
            log│V│ = logdet(V)
        end
        cache[Z] = (V⁻¹, log│V│)
        return cache[Z]
    end
end

function minv(G::AbstractMatrix, σ::Vector, Z::AbstractMatrix, cache::Dict)
    #h = hash(M)
    if Z in keys(cache)
        #return cache[h]
        return cache[Z]
    else
        V    = mulαβαtc(Z, G, Diagonal(Z*σ))
        #V   = Z * G * Z' + Diagonal(Z*σ)
        #V⁻¹  = inv(SHermitianCompact(SMatrix{size(V, 1),size(V, 1)}(V)))
        V⁻¹  = inv(V)
        #V⁻¹  = inv(V)
        #ldV = logdet(V)
        cache[Z] = V⁻¹
        return V⁻¹
    end
end

function mlogdet(M::Matrix, cache::Dict)
    #h = hash(M)
    if M in keys(cache)
        return cache[M]
    else
        iM = logdet(M)
        cache[M] = iM
        return iM
    end
end
#-------------------------------------------------------------------------------
#             REML FOR OPT ALGORITHM
#-------------------------------------------------------------------------------
"""
    -2 REML function for ForwardDiff
"""
function reml2(data::RBEDataStructure, θ::Vector{T}, β::Vector; memopt::Bool = true) where T
    #memory optimizations to reduse allocations (cache rebuild)
    #empty!(data.mem.dict)
    rebuildcache(data, promote_type(eltype(data.yv[1]), T))
    cache     = Dict{Matrix, Tuple{Matrix{T}, T}}()
    #cache     = data.mem.dict
    #---------------------------------------------------------------------------
    G         = gmat(view(θ, 3:5))
    θ₁        = zero(T)
    θ₂        = zeros(promote_type(eltype(data.yv[1]), eltype(θ)), data.p, data.p)
    θ₃        = zero(T)
    V⁻¹       = nothing
    #mVec      = pmap(x -> mvmatall(G, view(θ,1:2), x, first(data.mem.svec), cache), data.Zv)
    @simd for i = 1:data.n
        if MEMOPT && memopt

            V⁻¹, log│V│         = mvmatall(G, view(θ,1:2), data.Zv[i], first(data.mem.svec), cache)
            #V, V⁻¹, log│V│         =mVec[i]
            θ₁                    += log│V│
        else
            @inbounds V     = vmat(G, rmat(view(θ,1:2), data.Zv[i]), data.Zv[i])
            V⁻¹             = invchol(V)
            θ₁             += logdet(V)
        end
        #-----------------------------------------------------------------------
        #θ2 += Xv[i]'*iV*Xv[i]
        @inbounds mulαtβαinc!(θ₂, data.Xv[i], V⁻¹, first(data.mem.svec))
        #-----------------------------------------------------------------------
        #r    = yv[i] - Xv[i] * β
        #θ3  += r' * iV * r
        @inbounds θ₃  += mulθ₃(data.yv[i], data.Xv[i], β, V⁻¹, first(data.mem.svec))
    end
    return   θ₁ + logdet(θ₂) + θ₃ + data.remlc

end
"""
    -2 REML estimation with β recalculation for ForwardDiff
"""
function reml2bfd(data::RBEDataStructure, θ; memopt::Bool = true)
    return reml2b(data, θ; memopt = memopt)[1]
end
function reml2b(data::RBEDataStructure, θ::Vector{T}; memopt::Bool = true) where T

    rebuildcache(data, promote_type(eltype(data.yv[1]), T))
    cache     = Dict{Matrix, Tuple{Matrix{T}, T}}()
    #---------------------------------------------------------------------------
    G         = gmat(view(θ,3:5))
    V⁻¹       = Vector{Matrix{T}}(undef, data.n)                        # Vector of  V⁻¹ matrices
    #V         = nothing
    local log│V│::T                                                       # Vector log determinant of V matrix
    θ₁        = zero(T)
    θ₂        = zeros(promote_type(eltype(first(data.yv)), T), data.p, data.p)
    θ₃        = zero(T)
    βm        = zeros(promote_type(eltype(first(data.yv)), T), data.p)
    β         = zeros(promote_type(eltype(first(data.yv)), T), data.p)
    #mVec      = map(x -> mvmatall(G, θ[1:2], x, first(data.mem.svec), cache), data.Zv)
    @simd for i = 1:data.n
        if MEMOPT && memopt
            @inbounds V⁻¹[i], log│V│ = mvmatall(G, view(θ,1:2), data.Zv[i], first(data.mem.svec), cache)
            #V, V⁻¹[i], log│V│           = mVec[i]
            θ₁                         += log│V│
        else
            @inbounds R        = rmat(view(θ,1:2), data.Zv[i])
            @inbounds V        = vmat(G, R, data.Zv[i])
            @inbounds V⁻¹[i]   = invchol(V)
            θ₁                += logdet(V)
        end
        #-----------------------------------------------------------------------
        #θ2 += Xv[i]'*iVv[i]*Xv[i]
        #βm += Xv[i]'*iVv[i]*yv[i]
        mulθβinc!(θ₂, βm, data.Xv[i], V⁻¹[i], data.yv[i], first(data.mem.svec))
        #-----------------------------------------------------------------------
    end
    mul!(β, inv(θ₂), βm)
    for i = 1:data.n
        # r    = yv[i] - Xv[i] * β
        # θ3  += r' * iVv[i] * r
        @inbounds θ₃  += mulθ₃(data.yv[i], data.Xv[i], β, V⁻¹[i], first(data.mem.svec))
    end

    return   θ₁ + logdet(θ₂) + θ₃ + data.remlc,  β, θ₂
end
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# C matrix derivation
#-------------------------------------------------------------------------------
"""
non inverted C matrix gradient function
"""
function cmatgf(data, θ::Vector{T}; memopt::Bool = true) where T
    p      = size(data.Xv[1], 2)
    jC     = ForwardDiff.jacobian(x -> cmatvec(data, x; memopt = memopt),  θ)
    result = Vector{Matrix{T}}(undef, length(θ))
    for i in 1:length(θ)
        result[i] = reshape(view(jC, :, i), p, p) #<Opt
    end
    return result
end
"""
non inverted C matrix in vector form for gradient
"""
function cmatvec(data::RBEDataStructure, θ::Vector{T}; memopt::Bool = true) where T
    p     = size(data.Xv[1], 2)
    G     = gmat(view(θ, 3:5))
    C     = zeros(promote_type(eltype(data.Zv[1]), T), p, p)
    cache     = Dict()
    for i = 1:length(data.Xv)
        if memopt
            V⁻¹   = minv(G, θ[1:2], data.Zv[i], cache)
        else
            R   = rmat(θ[1:2], data.Zv[i])
            V⁻¹  = inv(mulαβαtc(data.Zv[i], G, R))
            #V⁻¹  = invchol(vmat(G, R, data.Zv[i]))
            #V⁻¹  = invchol!(mulαβαtc(data.Zv[i], G, R))

        end
        #C  += Xv[i]' * V⁻¹ * Xv[i]
        mulαtβαinc!(C, data.Xv[i], V⁻¹)
    end
    #return C[:]
    return C
end
"""
C matrix gradients
"""
function cmatg(data, θ::Vector, C::Matrix; memopt::Bool = true)
    g  = Vector{Matrix}(undef, length(θ))
    jC = cmatgf(data, θ; memopt = memopt)
    for i = 1:length(θ)
        g[i] = (- C * jC[i] * C)
    end
    return g
end
"""
L * C * L' for all C gradient marices
"""
function lclg(gradc, L)
    g  = Vector{eltype(gradc[1])}(undef, length(gradc))
    for i = 1:length(gradc)
        g[i] = (L * gradc[i] * L')[1]
    end
    return g
end
"""
"""
function contrastvec(data, res, L)
    lcl     = L*res.C*L'
    lclr    = rank(lcl)
    F       = res.β'*L'*inv(lcl)*L*res.β/lclr
    df      = sattdf(data, res.gradc, res.A, res.C, L, lcl)
    pval    = ccdf(FDist(lclr, df), F)
    return F, lclr, df, pval
end
function estimatevec(data, res, L)
    lcl     = L*res.C*L'
    β       = copy(res.β)
    est     = (L*β)[1]
    lclr    = rank(lcl)
    se      = sqrt((lcl)[1])
    t       = ((est)/se)
    return est, se, t
end
function sattdf(data, gradc, A, C, L, lcl)
    lclr    = rank(lcl)
    if lclr ≥ 2
        vm  = Vector{eltype(C)}(undef, lclr)
        # Spectral decomposition ?
        #ev  = eigen(lcl)
        #pl  = eigvecs(ev)
        #dm  = eigvals(ev)
        #ei  = pl * L
        for i = 1:lclr
            g         = lclg(gradc, L[i:i,:])
            #g         = lclg(res.gradc, ei[i:i, :])
            dfi       = 2*((L[i:i,:]*C*L[i:i,:]')[1])^2/(g'*A*g)
            #dfi       = 2*dm[i]^2/(g'*res.A*g)
            if dfi > 2
                vm[i] = dfi/(dfi-2)
            else
                vm[i] = 0
            end
        end
        E   = sum(vm)
        if E > lclr
            dfi = 2 * E / (E - lclr)
        else
            dfi = 0
        end
    else
        g   = lclg(gradc, L)
        dfi = 2*((lcl)[1])^2/(g'*A*g)
    end
    return max(1, dfi)
end
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
"""
    Initial variance computation
"""
function initvar(df::DataFrame, dv::Symbol, fac::Symbol)::Vector
    f  = unique(df[:, fac])
    fv = Array{eltype(df[!, dv]), 1}(undef, 0)
    for i in f
        push!(fv, var(df[df[!, fac] .== i, dv]))
    end
    return fv
end
#-------------------------------------------------------------------------------
function optimcallback(x)
    false
end
#-------------------------------------------------------------------------------
function vlink(σ::T) where T <: Real
    exp(σ)
end
function vlinkr(σ::T) where T <: Real
    log(σ)
end

function rholinkpsigmoid(ρ::T, m) where T <: Real
    return 1.0/(1.0 + exp(ρ * m))
end
function rholinkpsigmoidr(ρ::T, m) where T <: Real
    return log(1.0/ρ - 1.0)/m
end

function rholinksigmoid(ρ::T, m) where T <: Real
    return ρ/sqrt(1.0 + ρ^2)
end
function rholinksigmoidr(ρ::T, m) where T <: Real
    return sign(ρ)*sqrt(ρ^2/(1.0 - ρ^2))
end

function rholinksigmoid2(ρ::T, m) where T <: Real
    return atan(ρ)/pi*2.0
end
function rholinksigmoid2r(ρ::T, m) where T <: Real
    return tan(ρ*pi/2.0)
end

function varlinkmap(θ, r1::Union{Int, UnitRange}, r2::Union{Int, UnitRange}, f1::Function, f2::Function)
    #θl      = similar(θ)
    @inbounds @simd for i in r1
        θ[i]  = f1(θ[i])
    end
    #
    @inbounds @simd for i in r2
        θ[i]  = f2(θ[i])
    end
    return θ
end

@inline function lvecupd!(L::AbstractVector, fac)
    L .= 0.
    L[fac] .= 1.
end
