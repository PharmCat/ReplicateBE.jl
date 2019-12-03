"""
A' * B * A -> + θ (cache)
"""
function mulall!(θ, A::Matrix, B::Matrix, c)
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
function mulall!(θ, A::Matrix, B::Matrix)
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
"""
A' * B * A -> θ
A' * B * C -> β
"""
function mulall!(θ, β, A::Matrix, B::Matrix, C::Vector, c)
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
"""
A' * B * A
"""
function mulall(A::Vector, B::Matrix)
    q = size(B, 1)
    θ = 0
    #c .= 0
        for n = 1:q
            for m = 1:q
                #@inbounds c[n] += B[m, n] * A[m]
                θ += B[m, n] * A[m] * A[n]
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
function mulall(y::Vector, X::Matrix, β::Vector, V::Matrix, c)
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
            θ += V[m, n] * (y[m] - c[m]) * (y[n] - c[n])
        end
    end
    return θ
end
