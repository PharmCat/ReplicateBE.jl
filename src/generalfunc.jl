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
    #h = hash(Z)
    #if h in keys(cache)
    if Z in keys(cache)
        #return cache[h]
        return cache[Z]
    else
        V   = mulαβαtc(Z, G, Diagonal(Z*σ), mem)
        #V   = Z * G * Z' + Diagonal(Z*σ)
        iV  = inv(V)
        ldV = logdet(V)
        cache[Z] = (V, iV, ldV)
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
    -2 REML function for ForwardDiff
"""
function reml2(data::RBEDataStructure, θvec::Vector, β::Vector; memopt::Bool = true)

    #memory optimizations to reduse allocations (cache rebuild)
    #empty!(data.mem.dict)
    rebuildcache(data, promote_type(eltype(data.yv[1]), eltype(θvec)))
    cache     = Dict{Matrix, Tuple{Matrix, Matrix, Number}}()
    #cache     = data.mem.dict
    #---------------------------------------------------------------------------
    G         = gmat(θvec[3:5])
    θ1        = 0
    θ2        = zeros(promote_type(eltype(data.yv[1]), eltype(θvec)), data.p, data.p)
    θ3        = 0
    iV        = nothing
    for i = 1:data.n
        if MEMOPT && memopt

            V, iV, ldV         = mvmatall(G, θvec[1:2], data.Zv[i], first(data.mem.svec), cache)
            θ1                += ldV
        else
            @inbounds V    = vmat(G, rmat(θvec[1:2], data.Zv[i]), data.Zv[i])
            iV             = inv(V)
            θ1            += logdet(V)
        end
        #-----------------------------------------------------------------------
        #θ2 += Xv[i]'*iV*Xv[i]
        @inbounds mulαtβαinc!(θ2, data.Xv[i], iV, first(data.mem.svec))
        #-----------------------------------------------------------------------
        #r    = yv[i] - Xv[i] * β
        #θ3  += r' * iV * r
        @inbounds θ3  += mulθ₃(data.yv[i], data.Xv[i], β, iV, first(data.mem.svec))
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
    cache     = Dict()
    #---------------------------------------------------------------------------
    G         = gmat(θvec[3:5])
    iVv       = Array{Array{eltype(θvec), 2}, 1}(undef, data.n)
    V         = nothing
    ldV       = nothing
    θ1        = 0
    θ2        = zeros(promote_type(eltype(first(data.yv)), eltype(θvec)), data.p, data.p)
    θ3        = 0
    βm        = zeros(promote_type(eltype(first(data.yv)), eltype(θvec)), data.p)
    β         = zeros(promote_type(eltype(first(data.yv)), eltype(θvec)), data.p)
    for i = 1:data.n
        if MEMOPT && memopt

            @inbounds V, iVv[i], ldV = mvmatall(G, θvec[1:2], data.Zv[i], first(data.mem.svec), cache)
            θ1                      += ldV
        else
            @inbounds R        = rmat(θvec[1:2], data.Zv[i])
            @inbounds V        = vmat(G, R, data.Zv[i])
            @inbounds iVv[i]   = inv(V)
            θ1                += logdet(V)
        end
        #-----------------------------------------------------------------------
        #θ2 += Xv[i]'*iVv[i]*Xv[i]
        #βm += Xv[i]'*iVv[i]*yv[i]
        mulθβinc!(θ2, βm, data.Xv[i], iVv[i], data.yv[i], first(data.mem.svec))
        #-----------------------------------------------------------------------
    end
    mul!(β, inv(θ2), βm)
    for i = 1:data.n
        # r    = yv[i] - Xv[i] * β
        # θ3  += r' * iVv[i] * r
        @inbounds θ3  += mulθ₃(data.yv[i], data.Xv[i], β, iVv[i], first(data.mem.svec))
    end

    return   θ1 + logdet(θ2) + θ3 + data.remlc,  β, θ2
end
#-------------------------------------------------------------------------------

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
        mulαtβαinc!(C, Xv[i], iV)
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
    df      = sattdf(data, res, L, lcl)
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
function sattdf(data, res, L, lcl)
    lclr    = rank(lcl)
    if lclr ≥ 2
        vm  = Vector{eltype(res.C)}(undef, lclr)
        # Spectral decomposition ?
        #ev  = eigen(lcl)
        #pl  = eigvecs(ev)
        #dm  = eigvals(ev)
        #ei  = pl * L
        for i = 1:lclr
            g         = lclg(res.gradc, L[i:i,:])
            #g         = lclg(res.gradc, ei[i:i, :])
            dfi       = 2*((L[i:i,:]*res.C*L[i:i,:]')[1])^2/(g'*res.A*g)
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
        g   = lclg(res.gradc, L)
        dfi = 2*((lcl)[1])^2/(g'*res.A*g)
    end
    return max(1, dfi)
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

function varlinkmap(θ::Vector, r1::Union{Int, UnitRange}, r2::Union{Int, UnitRange}, f1::Function, f2::Function)
    θl      = similar(θ)
    θl[r1]  = f1.(θ[r1])
    θl[r2]  = f2.(θ[r2])
    return θl
end
