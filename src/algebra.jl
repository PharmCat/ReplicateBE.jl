"""
A' * B * A -> + θ (cache)
"""
function mulαtβαinc!(θ, A::AbstractMatrix, B::Matrix, c)
    q = size(B, 1)
    p = size(A, 2)
    for i = 1:p
        c .= 0
        for n = 1:q
            for m = 1:q
                @inbounds c[n] += B[m, n] * A[m, i]
            end
        end
        for n = 1:p
            for m = 1:q
                @inbounds θ[i, n] += A[m, n] * c[m]
            end
        end
    end
    θ
end
"""
A' * B * A -> + θ
"""
function mulαtβαinc!(θ, A::AbstractMatrix, B::AbstractMatrix)
    q = size(B, 1)
    p = size(A, 2)
    c = zeros(eltype(B), q)
    for i = 1:p
        c .= 0
        for n = 1:q
            for m = 1:q
                @inbounds c[n] += B[m, n] * A[m, i]
            end
        end
        for n = 1:p
            for m = 1:q
                @inbounds θ[i, n] += A[m, n] * c[m]
            end
        end
    end
    θ
end
#-------------------------------------------------------------------------------
"""
A' * B * A -> +θ
A' * B * C -> +β
"""
function mulθβinc!(θ, β, A::AbstractMatrix, B::AbstractMatrix, C::Vector, c)
    q = size(B, 1)
    p = size(A, 2)
    for i = 1:p
        c .= 0
        for n = 1:q
            for m = 1:q
                @inbounds c[n] += B[m, n] * A[m, i]
            end
            @inbounds β[i] += c[n] * C[n]
        end

        for n = 1:p
            for m = 1:q
                @inbounds θ[i, n] += A[m, n] * c[m]
            end
        end
    end
    θ, β
end
#-------------------------------------------------------------------------------
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
(y - X * β)' * V * (y - X * β)
"""
function mulθ₃(y::Vector, X::AbstractMatrix, β::Vector, V::AbstractMatrix, c)
    q = size(V, 1)
    p = size(X, 2)
    θ = 0
    c .= 0
    for n = 1:q
        for m = 1:p
            @inbounds c[n] += X[n, m] * β[m]
        end
    end
    for n = 1:q
        for m = 1:q
            @inbounds θ += V[m, n] * (y[m] - c[m]) * (y[n] - c[n])
        end
    end
    return θ
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
A * B * A' + C
"""
function mulαβαtc(A::AbstractMatrix, B::AbstractMatrix, C::AbstractMatrix)
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
    mx .+= C
end
function mulαβαtc(A::AbstractMatrix, B::AbstractMatrix, C::AbstractMatrix, c::AbstractVector)
    q  = size(B, 1)
    p  = size(A, 1)
    #c  = mem.svec[p]
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
    mx .+= C
    #SMatrix{p,p,eltype(mx)}(mx)
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

function invchol(M)
    q  = size(M, 1)
    v  = zeros(eltype(M), q, q)
    if q == 1
        v[1,1] = 1/M[1,1]
        return v
    end
    il = inv(cholesky(M).U)
    for n = 1:q
        for m = n:q
            @inbounds v[n, n] += il[n, m]^2
        end
    end
    for n = 1:(q-1)
        for m = (n+1):q
            for i = m:q
                @inbounds v[n, m] += il[n, i] * il[m, i]
            end
        end
    end
    for n = 1:q-1
        for m = 2:q
            @inbounds v[m, n] = v[n, m]
        end
    end
    return v
end
