"""
A' * B * A
"""
function mulall(A::Vector, B::AbstractMatrix)
    q = size(B, 1)
    θ = 0
    #c .= 0
        for n = 1:q
            for m = 1:q
                #@inbounds c[n] += B[m, n] * A[m]
                @inbounds θ += B[m, n] * A[m] * A[n]
            end
            #@inbounds θ += A[n] * c[n]
        end
        #for n = 1:q
        #    @inbounds θ += A[n] * c[n]
        #end
    θ
end

"""
A * B * A'
"""
function mulαβαt(A::AbstractMatrix, B::AbstractMatrix)
    q  = size(B, 1)
    p  = size(A, 1)
    c  = zeros(eltype(B), q)
    mx = zeros(eltype(B), p, p)
    for i = 1:p
        c .= 0
        for n = 1:q
            for m = 1:q
                @inbounds c[n] +=  A[i, m] * B[n, m]
            end
        end
        for n = 1:p
            for m = 1:q
                 @inbounds mx[i, n] += A[n, m] * c[m]
            end
        end
    end
    mx
end

"""
A * B * A' + C -> O
"""
function mulαβαtc!(O::AbstractMatrix, A::AbstractMatrix, B::AbstractMatrix, C::AbstractMatrix)
    q  = size(B, 1)
    p  = size(A, 1)
    c  = zeros(eltype(B), q)
    O .= 0
    for i = 1:p
        c .= 0
        for n = 1:q
            for m = 1:q
                @inbounds c[n] +=  A[i, m] * B[n, m]
            end
        end
        for n = 1:p
            for m = 1:q
                 @inbounds O[i, n] += A[n, m] * c[m]
            end
        end
    end
    O .+= C
end
function mulαβαtc!(O::AbstractMatrix, A::AbstractMatrix, B::AbstractMatrix, C::AbstractMatrix, mem::MemCache)
    q  = size(B, 1)
    p  = size(A, 1)
    c  = mem.svec[p]
    O .= 0
    for i = 1:p
        c .= 0
        for n = 1:q
            for m = 1:q
                @inbounds c[n] +=  A[i, m] * B[n, m]
            end
        end
        for n = 1:p
            for m = 1:q
                 @inbounds O[i, n] += A[n, m] * c[m]
            end
        end
    end
    O .+= C
end
#-------------------------------------------------------------------------------
#Return term levels count by symbol
function termmodelleveln(MF::ModelFrame, symbol::Symbol)::Int
    id = findterm(MF, symbol)
    return length(MF.f.rhs.terms[id].contrasts.levels)
end
#=

"""
    lsm(rbe::RBE, L::Matrix)

Deprecated.
"""
#Deprecated
function lsm(rbe::RBE, L::Matrix)
    lcl  = L*rbe.C*L'
    return L*coef(rbe), sqrt.(lcl)
end
#
"""
    emm(obj::RBE, fm::Matrix, lm::Matrix)

Matrix mask.
"""
function emm(obj::RBE, fm::Matrix, lm::Matrix)
    La = lmean(obj::RBE)
    L  = La .* fm
    L  = L  .+ lm
    return lsm(obj, Matrix(L))
end
#General mean contrast L matrix 1xp
"""
    lmean(obj::RBE)

Return L-matrix for general mean.
"""
function lmean(obj::RBE)
    L    = zeros(1, length(fixed(obj).est))
    L[1] = 1.0
    it    = 2
    for f in obj.data.factors
        term = findterm(obj.model, f)
        len  = length(obj.model.f.rhs.terms[term].contrasts.termnames)
        dev  = 1/length(obj.model.f.rhs.terms[term].contrasts.levels)
        for i = 1:len
            L[it] = dev
            it  += 1
        end
    end
    return L
end

=#
#-------------------------------------------------------------------------------

#=
function randrbeds(n::Int, sequence::Vector,
    design::Matrix,
    inter::Real, intra::Vector,
    intercept::Real, seqcoef::Vector, periodcoef::Vector, formcoef::Vector,
    dropobs::Int, seed)
    if seed != 0
        rng = MersenneTwister(seed)
    else
        rng = MersenneTwister()
    end

    r = n/sum(sequence)
    sn = Array{Int, 1}(undef, length(sequence))
    for i = 1:(length(sequence)-1)
        sn[i] = round(r*sequence[i])
    end
    sn[length(sequence)] = n - sum(sn[1:(length(sequence)-1)])

    u      = unique(design)
    sqname = Array{String, 1}(undef,size(design)[1])
    sqnum  = size(design)[1]
    pnum   = size(design)[2]
    for i = 1:sqnum
        sqname[i] = join(design[i,:])
    end
    Zv = Array{Matrix, 1}(undef, sqnum)
    Vv = Array{Vector, 1}(undef, sqnum)
    for i = 1:size(design)[1]
        Z = Array{Int, 2}(undef, pnum, length(u))
        for c = 1:pnum
            for uc = 1:length(u)
                if design[i, c] == u[uc] Z[c, uc] = 1 else Z[c, uc] = 0 end
            end
        end
        Zv[i] = Z
        Vv[i] = Z * intra
    end
    Mv = Array{Array{Float64, 1}, 1}(undef, sqnum)
    for i = 1:sqnum
        Mv[i] = zeros(pnum) .+ intercept .+ seqcoef[i] + periodcoef + Zv[i]*formcoef
    end
    ndist  = Normal()
    subjds = DataFrame(subject = Int[], formulation = String[], period = Int[], sequence = String[], var = Float64[])
    subj   = 1
    subjmx = Array{Any, 2}(undef, pnum, 5)
    for i = 1:sqnum
        for sis = 1:sn[i]
            subjmx[:, 1] .= subj
            subjmx[:, 2]  = design[i,:]
            subjmx[:, 3]  = collect(1:pnum)
            subjmx[:, 4] .= sqname[i]
            subjmx[:, 5] .= 0
            subjmx[:, 5] .+= rand(rng, ndist)*sqrt(inter)
            subj += 1
            for c = 1:pnum
                subjmx[c, 5] += Mv[i][c] + rand(rng, ndist)*sqrt(Vv[i][c])
                push!(subjds, subjmx[c, :])
            end
        end
    end
    if dropobs > 0 && dropobs < size(subjds, 1)
        dellist = sample(rng, 1:size(subjds, 1), dropobs, replace = false)
        deleterows!(subjds, sort!(dellist))
    end
    categorical!(subjds, :subject);
    categorical!(subjds, :formulation);
    categorical!(subjds, :period);
    categorical!(subjds, :sequence);
    return subjds
end
=#

#-------------------------------------------------------------------------------
#=
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


struct MemAlloc{T <: AbstractFloat}
    mem1::Vector{Matrix{T}}
    mem2::Vector{Matrix{T}}
    mem3::Vector{Vector{T}}
    #mem4::Array{Float64, 2}
    function MemAlloc(p, zs, yv::Vector{Vector{T}}) where T <: AbstractFloat
        maxobs  = maximum(length.(yv))
        yvtype = eltype(yv[1])
        memc1 = Vector{Matrix{yvtype}}(undef, maxobs)
        memc2 = Vector{Matrix{yvtype}}(undef, maxobs)
        memc3 = Vector{Vector{yvtype}}(undef, maxobs)
        for i = 1:maxobs
            memc1[i] = zeros(i, zs)
            memc2[i] = zeros(p, i)
            memc3[i] = zeros(i)
        end
        #memc4 = zeros(p, p)
        new{T}(memc1, memc2, memc3)::MemAlloc
    end
end


"""
    Return C matrix
    var(β) p×p variance-covariance matrix
"""
@inline function cmat(Xv::Vector{Matrix{T}}, Zv::Vector, iVv::Vector, θ::Vector)::Matrix where T <: AbstractFloat
    p = size(Xv[1])[2]
    C = zeros(p, p)
    for i=1:length(Xv)
        #@inbounds C .+= Xv[i]' * iVv[i] * Xv[i]
        mulαtβαinc!(C, Xv[i], iVv[i])
    end
    return pinv(C)
end



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
    βm   = zeros(data.p)
    cache     = Dict()
    cachel    = Dict()
    @inbounds for i = 1:data.n
        #rmat!(Rv[i], θvec[1:2], data.Zv[i])
        #Memopt!
        Vv[i], iVv[i], ldV         = mvmatall(G, θvec[1:2], data.Zv[i], first(data.mem.svec), cache)
        θ1  += ldV
        #-
        #θ2 += Xv[i]'*iVv[i]*Xv[i]
        mulθβinc!(θ2, βm, data.Xv[i], iVv[i], data.yv[i], first(data.mem.svec))
    end
    mul!(β, inv(θ2), βm)
    for i = 1:data.n
        #Same:
        #r    = yv[i] - Xv[i]*β
        #θ3  += r'*iVv[i]*r
        @inbounds θ3  += mulθ₃(data.yv[i], data.Xv[i], β, iVv[i], first(data.mem.svec))
    end
    return   θ1 + logdet(θ2) + θ3 + data.remlc
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
        #C  += Xv[i]' * iV * Xv[i]
        mulall!(C, Xv[i], iV)
    end
    return (L * inv(C) * Lt)[1]
end
=#
#=
function mrmat(σ::Vector{S}, Z::Matrix{T}, cache)::Matrix where S <: Real where T <: Real
    h = hash(tuple(σ, Z))
    if h in keys(cache)
        return cache[h]
    else
        if σ[1] < 0.0 σ[1] = 1.0e-6 end
        if σ[2] < 0.0 σ[2] = 1.0e-6 end
        R = Matrix(Diagonal((Z*σ)[:,1]))
        cache[h] = R
        return R
    end
end
=#
#=
function Base.show(io::IO, obj::Tuple{Vararg{Tuple{Real, Real}}})
    for i in obj
        println(io, i)
    end
end
=#
# MEMOIZATION
#=
@memoize function memrmat(σ::Vector, Z::Matrix)::Matrix
    return rmat(σ, Z)
end
@memoize function memzgz(G::Matrix, Z::Matrix)::Matrix
    return Z*G*Z'
end
@memoize function memvmat(ZGZ::Matrix, R::Matrix)::Matrix
    return ZGZ + R
end
@memoize function meminv(m::Matrix)::Matrix
    return inv(m)
end
=#
#=

#M       = L'*inv(L*L')*L
#t1      = tr(M*C)
#v1      = t1^2/tr(M*C*M*C)
#F[i]    = β'*M*β/t1
#df[i]   = 2*(t1)^2/(g'*(A)*g)

"""
    Secondary param estimation:
    SE
    F
    DF
    C
"""
function ctrst(p, Xv, Zv, iVv, θ, β, A, C; memopt::Bool = true)
    se    = Array{Real, 1}(undef, p)
    F     = Array{Real, 1}(undef, p)
    df    = Array{Real, 1}(undef, p)
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
    Ra = Array{Array{Real,2}, 1}(undef, n)
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
    Va = Array{Array{Real,2}, 1}(undef, n)
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
function βcoef!(p::Int, n::Int, yv::Array{Array{Real, 1}, 1}, Xv::Array{Array{Real, 2}, 1}, iVv::Array{Array{Real, 2}, 1}, β::Array{Real, 1})
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
function reml2!(yv::S, Zv::T, p::Int, n::Int, N::Int, Xv::T, G::Array{Real, 2}, Rv::T, Vv::T, iVv::T, θvec::Array{Real, 1}, β::Array{Real, 1}, memc, memc2, memc3, memc4)::Real where T <: Array{Array{Real, 2}, 1} where S <: Array{Array{Real, 1}, 1}
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
    memc = Array{Array{Real, 2}, 1}(undef, maxobs)
    memc2 = Array{Array{Real, 2}, 1}(undef, maxobs)
    memc3 = Array{Array{Real, 1}, 1}(undef, maxobs)
    for i = 1:maxobs
        memc[i] = zeros(i, zs)
        memc2[i] = zeros(p, i)
        memc3[i] = zeros(i)
    end
    memc4 = zeros(p, p)
    return memc, memc2, memc3, memc4
end

#Model Frame utils
function lvec(mm::ModelMatrix, f::Int)
    l = zeros(length(mm.assign))
    for i = 1:length(l)
        if mm.assign == f l[i] = 1 end
    end
end

=#
