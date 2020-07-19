#-------------------------------------------------------------------------------
#                               GENERAL FUNCTIONS
#
#-------------------------------------------------------------------------------
"""
    Make X, Z matrices and vector y for each subject;
"""
function sortsubjects(df::DataFrame, sbj::Symbol, X::Matrix, Z::Matrix, y::Vector)
    u = unique(df[!, sbj])
    Xa = Vector{Matrix{eltype(X)}}(undef, length(u))
    Za = Vector{Matrix{eltype(Z)}}(undef, length(u))
    ya = Vector{Vector{eltype(y)}}(undef, length(u))
    for i = 1:length(u)
        v = findall(x->x==u[i], df[!, sbj])

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
    cov = σ[1] * σ[2] * σ[3]
    return Symmetric([σ[1]^2 cov; cov σ[2]^2])
end
"""
    G matrix  (memory pre-allocation)
"""
@inline function gmat!(G::Matrix{T}, σ::Vector) where T <: AbstractFloat
    G[1, 1] = σ[1]^2
    G[2, 2] = σ[2]^2
    G[1, 2] = G[2, 1] = σ[1] * σ[2] * σ[3]
    return
end

"""
    R matrix (ForwardDiff+)
"""
@inline function rmat(σ::Vector, Z::Matrix)::Matrix
    return Diagonal(Z*(σ.^2))
end
"""
    R matrix  (memory pre-allocation)
"""
@inline function rmat!(R::AbstractMatrix{T}, σ::Vector{T}, Z::Matrix{T}) where T <: AbstractFloat
    copyto!(R, Diagonal(Z*(σ.^2)))
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
    #h = hash(tuple(σ, Z))
    if Z in keys(cache)
        return cache[Z]
    else
        V  = mulαβαtc(Z, G, rmat(σ, Z))
        #V   = Z * G * Z' + Diagonal(Z*σ)
        cache[Z] = V
        return V
    end
end
function mvmat(G::AbstractMatrix, σ::Vector, Z::Matrix, mem, cache)::Matrix
    #h = hash(tuple(σ, Z))
    if Z in keys(cache)
        return cache[Z]
    else
        V  = mulαβαtc(Z, G, rmat(σ, Z), mem)
        #V   = Z * G * Z' + Diagonal(Z*σ)
        cache[Z] = V
        return V
    end
end

function mvmatall(G::AbstractMatrix, theta, Z::Matrix, vs, mem, cache)
    #h = hash(Z)
    #if h in keys(cache)
    if Z in keys(cache)
        #return cache[h]
        return cache[Z]
    else
        V   = mulαβαtc(Z, G, rmat(theta, Z, vs), mem)
        #V   = Z * G * Z' + Diagonal(Z*σ)
        V⁻¹ = nothing
        try
            V⁻¹  = invchol(V)
        catch e
            if typeof(e) <: PosDefException
                V⁻¹  = inv(V)
            else
                throw(e)
            end
        end
        #V⁻¹  = inv(V)
        log│V│   = logdet(V)
        cache[Z] = (V, V⁻¹, log│V│)
        return V, V⁻¹, log│V│
    end
end
#println("θ₁: ", θ1, " θ₂: ",  θ2,  " θ₃: ", θ3)

function minv(G::AbstractMatrix, theta, Z::Matrix, vs, cache::Dict)::Matrix
    #h = hash(M)
    if Z in keys(cache)
        #return cache[h]
        return cache[Z]
    else
        V    = mulαβαtc(Z, G, rmat(theta, Z, vs))
        #V   = Z * G * Z' + Diagonal(Z*σ)
        V⁻¹  = invchol(V)
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
function reml2(data::RBEDataStructure, θ, β::Vector; memopt::Bool = true)
    #memory optimizations to reduse allocations (cache rebuild)
    #empty!(data.mem.dict)
    rebuildcache(data, promote_type(eltype(data.yv[1]), eltype(θ)))
    cache     = Dict{Matrix, Tuple{Matrix, Matrix, Number}}()
    #cache     = data.mem.dict
    #---------------------------------------------------------------------------
    #G         = gmat(θ[3:5])
    G         = gmat(θ, data.vs)
    θ₁        = 0
    θ₂        = zeros(promote_type(eltype(data.yv[1]), eltype(θ)), data.p, data.p)
    θ₃        = 0
    V⁻¹       = nothing
    for i = 1:data.n
        if MEMOPT && memopt

            V, V⁻¹, log│V│         = mvmatall(G, θ, data.Zv[i], data.vs, first(data.mem.svec), cache)
            θ₁                    += log│V│
        else
            @inbounds V     = vmat(G, rmat(θ, data.Zv[i], data.vs), data.Zv[i])
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
function reml2b(data::RBEDataStructure, θ; memopt::Bool = true)

    rebuildcache(data, promote_type(eltype(data.yv[1]), eltype(θ)))
    cache     = Dict()
    #---------------------------------------------------------------------------
    #G         = gmat(θ[3:5])
    G         = gmat(θ, data.vs)
    V⁻¹       = Vector{Matrix{eltype(θ)}}(undef, data.n)                        # Vector of  V⁻¹ matrices
    V         = nothing
    log│V│    = nothing                                                         # Vector log determinant of V matrix
    θ₁        = 0
    θ₂        = zeros(promote_type(eltype(first(data.yv)), eltype(θ)), data.p, data.p)
    θ₃        = 0
    βm        = zeros(promote_type(eltype(first(data.yv)), eltype(θ)), data.p)
    β         = zeros(promote_type(eltype(first(data.yv)), eltype(θ)), data.p)
    for i = 1:data.n
        if MEMOPT && memopt
            @inbounds V, V⁻¹[i], log│V│ = mvmatall(G, θ, data.Zv[i], data.vs, first(data.mem.svec), cache)
            θ₁                         += log│V│
        else
            @inbounds R        = rmat(θ, data.Zv[i], data.vs)
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
function cmatgf(data::RBEDataStructure, θ::Vector; memopt::Bool = true)
    #p      = size(data.Xv[1], 2)
    jC     = ForwardDiff.jacobian(x -> cmatvec(data, x; memopt = memopt),  SVector{length(θ), eltype(θ)}(θ))
    result = Vector{Matrix}(undef, 0)
    for i in 1:length(θ)
        push!(result, reshape(jC[:,i], data.p, data.p))
    end
    return result
end
"""
non inverted C matrix in vector form for gradient
"""
function cmatvec(data::RBEDataStructure, θ; memopt::Bool = true)
    #p     = size(Xv[1], 2)
    #G     = gmat(θ[3:5])
    G     = gmat(θ, data.vs)
    C     = zeros(promote_type(eltype(data.Zv[1]), eltype(θ)), data.p, data.p)
    cache     = Dict()
    #cachem    = Dict()
    for i = 1:length(data.Xv)
        if memopt
            V⁻¹   = minv(G, θ, data.Zv[i], data.vs, cache)
        else
            #R    = rmat(θ[1:2], data.Zv[i])
            R    = rmat(θ, data.Zv[i], data.vs)
            V⁻¹  = invchol(vmat(G, R, data.Zv[i]))
        end
        #C  += Xv[i]' * V⁻¹ * Xv[i]
        mulαtβαinc!(C, data.Xv[i], V⁻¹)
    end
    return C[:]
end
"""
C matrix gradients
"""
function cmatg(data::RBEDataStructure, θ::Vector, C::Matrix; memopt::Bool = true)
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

function varlinkmap(θ, r1::Union{Int, UnitRange}, r2::Union{Int, UnitRange}, f1::Function, f2::Function)
    θl      = similar(θ)
    θl[r1]  = f1.(θ[r1])
    θl[r2]  = f2.(θ[r2])
    return θl
end
