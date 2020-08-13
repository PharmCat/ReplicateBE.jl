"""
A' * B * A -> + θ (cache)
"""
function mulαtβαinc!(θ, A::AbstractMatrix, B::AbstractMatrix, c)
    q = size(B, 1)
    p = size(A, 2)
    for i = 1:p
        c .= 0
        @inbounds for n = 1:q, m = 1:q
            c[n] += B[m, n] * A[m, i]
        end
        @inbounds for n = 1:p, m = 1:q
            θ[i, n] += A[m, n] * c[m]
        end
    end
    #θ
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
        @inbounds for n = 1:q, m = 1:q
            c[n] += B[m, n] * A[m, i]
        end
        @inbounds for n = 1:p, m = 1:q
            θ[i, n] += A[m, n] * c[m]
        end
    end
    #θ
end
#-------------------------------------------------------------------------------
"""
A' * B * A -> +θ (cache)
A' * B * C -> +β
"""
function mulθβinc!(θ, β, A::AbstractMatrix, B::AbstractMatrix, C::AbstractVector, c)
    q = size(B, 1)
    p = size(A, 2)
    for i = 1:p
        c .= 0
        @simd for n = 1:q
            @simd for m = 1:q
                @inbounds c[n] += B[m, n] * A[m, i]
            end
            @inbounds β[i] += c[n] * C[n]
        end
        @simd for n = 1:p
            @simd for m = 1:q
                @inbounds θ[i, n] += A[m, n] * c[m]
            end
        end
    end
    #θ, β
end
#-------------------------------------------------------------------------------
"""
(y - X * β)' * V * (y - X * β) (cache)
"""
function mulθ₃(y::AbstractVector, X::AbstractMatrix, β::AbstractVector, V::AbstractMatrix, c)
    q = size(V, 1)
    p = size(X, 2)
    θ = 0
    c .= 0
    @simd for n = 1:q
        @simd for m = 1:p
            @inbounds c[n] += X[n, m] * β[m]
        end
    end
    @simd for n = 1:q
        @simd for m = 1:q
            @inbounds θ += V[m, n] * (y[m] - c[m]) * (y[n] - c[n])
        end
    end
    return θ
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
        @simd for n = 1:q
            @simd for m = 1:q
                @inbounds c[n] +=  A[i, m] * B[n, m]
            end
        end
        @simd for n = 1:p
            @simd for m = 1:q
                 @inbounds mx[i, n] += A[n, m] * c[m]
            end
        end
    end
    mx .+= C
end
"""
A * B * A' + C (cache)
"""
function mulαβαtc(A::AbstractMatrix, B::AbstractMatrix, C::AbstractMatrix, c::AbstractVector)
    q  = size(B, 1)
    p  = size(A, 1)
    #c  = mem.svec[p]
    mx = zeros(eltype(B), p, p)
    for i = 1:p
        c .= 0
        @simd for n = 1:q
            @simd for m = 1:q
                @inbounds c[n] +=  A[i, m] * B[n, m]
            end
        end
        @simd for n = 1:p
            @simd for m = 1:q
                 @inbounds mx[i, n] += A[n, m] * c[m]
            end
        end
    end
    mx .+= C
end
"""
A * B * A' + Diagonal(A*C) (cache)
"""
function mulαβαtc(A::AbstractMatrix, B::AbstractMatrix, C::AbstractVector, c::AbstractVector)
    q  = size(B, 1)
    p  = size(A, 1)
    #c  = mem.svec[p]
    mx = zeros(eltype(B), p, p)
    for i = 1:p
        c .= 0
        @simd for n = 1:q
            @simd for m = 1:q
                @inbounds c[n] +=  A[i, m] * B[n, m]
            end
        end
        @simd for n = 1:p
            @simd for m = 1:q
                 @inbounds mx[i, n] += A[n, m] * c[m]
            end
        end
        @simd for m = 1:length(C)
             @inbounds mx[i, i] += A[i, m] * C[m]
        end
    end
    mx
end
"""
mx <- A * B * A' + Diagonal(A*C) (cache)
"""
function mulαβαtcupd!(mx::AbstractMatrix, A::AbstractMatrix, B::AbstractMatrix, C::AbstractVector, c::AbstractVector)
    q  = size(B, 1)
    p  = size(A, 1)
    mx .= 0
    for i = 1:p
        c .= 0
        @simd for n = 1:q
            @simd for m = 1:q
                @inbounds c[n] +=  A[i, m] * B[n, m]
            end
        end
        @simd for n = 1:p
            @simd for m = 1:q
                 @inbounds mx[i, n] += A[n, m] * c[m]
            end
        end
        @simd for m = 1:length(C)
             @inbounds mx[i, i] += A[i, m] * C[m]
        end
    end
    mx
end
"""
    Cholesky inverse
"""
function invchol(M::AbstractMatrix)
    q  = size(M, 1)
    v  = zeros(eltype(M), q, q)
    if q == 1
        v[1,1] = 1/M[1,1]
        return v
    end
    il = inv(cholesky(M).U)
    @simd for n = 1:q
        for m = n:q
            @inbounds v[n, n] += il[n, m]^2
        end
    end
    for n = 1:(q-1)
        for m = (n+1):q
            @simd for i = m:q
                @inbounds v[n, m] += il[n, i] * il[m, i]
            end
        end
    end
    #@simd for n = 1:q-1
    #    @simd for m = 2:q
    #        @inbounds v[m, n] = v[n, m] #Make Symmetric
    #    end
    #end
    return Symmetric(v, :U)
end
