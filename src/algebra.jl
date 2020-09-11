"""
A' * B * A -> + θ (cache)
"""
function mulαtβαinc!(θ, A::AbstractMatrix, B, c)
    q = size(B, 1)
    p = size(A, 2)
    for i = 1:p
        fill!(c, zero(eltype(θ)))
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
    c = zeros(eltype(θ), q)
    for i = 1:p
        fill!(c, zero(eltype(θ)))
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
A' * B * A -> +θ (cache)
A' * B * C -> +β
"""
function mulθβinc!(θ, β, A::AbstractMatrix, B::AbstractMatrix, C::AbstractVector, c)
    q = size(B, 1)
    p = size(A, 2)
    for i = 1:p
        fill!(c, zero(eltype(θ)))
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
(y - X * β)' * V * (y - X * β) (cache)
"""
function mulθ₃(y, X, β, V, c)
    q = size(V, 1)
    p = size(X, 2)
    θ = zero(eltype(β))
    #c .= 0
    fill!(c, zero(eltype(β)))
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
A * B * A' + C
"""
function mulαβαtc(A, B, C)::Symmetric
    q  = size(B, 1)
    p  = size(A, 1)
    c  = zeros(eltype(B), q)
    mx = zeros(eltype(B), p, p)
    for i = 1:p
        fill!(c, zero(eltype(c)))
        @simd for n = 1:q
            @simd for m = 1:q
                @inbounds c[n] +=  A[i, m] * B[n, m]
            end
        end
        @simd for n = i:p
            @simd for m = 1:q
                 @inbounds mx[i, n] += A[n, m] * c[m]
            end
            @inbounds mx[i, n] += C[i, n]
        end
    end
    #mx .+= C
    Symmetric(mx)
end
"""
A * B * A' + C (cache)
"""
function mulαβαtc(A::AbstractMatrix, B::AbstractMatrix, C::AbstractMatrix, c::AbstractVector)
    q  = size(B, 1)
    p  = size(A, 1)
    mx = zeros(eltype(B), p, p)
    for i = 1:p
        fill!(c, zero(eltype(c)))
        @simd for n = 1:q
            @simd for m = 1:q
                @inbounds c[n] +=  A[i, m] * B[n, m]
            end
        end
        @simd for n = i:p
            @simd for m = 1:q
                 @inbounds mx[i, n] += A[n, m] * c[m]
            end
            @inbounds mx[i, n] += C[i, n]
        end
    end
    #mx .+= C
    Symmetric(mx)
end
"""
A * B * A' + Diagonal(A*C) (cache)
"""
function mulαβαtc(A::AbstractMatrix, B::AbstractMatrix, C::AbstractVector, c::AbstractVector)::Symmetric
    q  = size(B, 1)
    p  = size(A, 1)
    mx = zeros(eltype(B), p, p)
    for i = 1:p
        fill!(c, zero(eltype(c)))
        @simd for n = 1:q
            @simd for m = 1:q
                @inbounds c[n] +=  A[i, m] * B[n, m]
            end
        end
        @simd for n = i:p
            @simd for m = 1:q
                 @inbounds mx[i, n] += A[n, m] * c[m]
            end
        end
        @simd for m = 1:length(C)
             @inbounds mx[i, i] += A[i, m] * C[m]
        end
    end
    #mx
    Symmetric(mx)
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
