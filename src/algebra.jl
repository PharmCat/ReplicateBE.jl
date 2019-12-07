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
            θ += V[m, n] * (y[m] - c[m]) * (y[n] - c[n])
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
                c[n] +=  A[i, m] * B[n, m]
            end
        end
        for n = 1:p
            for m = 1:q
                 mx[i, n] += A[n, m] * c[m]
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
                c[n] +=  A[i, m] * B[n, m]
            end
        end
        for n = 1:p
            for m = 1:q
                 mx[i, n] += A[n, m] * c[m]
            end
        end
    end
    mx .+ C
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
                try
                c[n] +=  A[i, m] * B[n, m]
                catch
                    println("$q, $p, $i, $n, $m" )
                    error()
                end
            end
        end
        for n = 1:p
            for m = 1:q
                 mx[i, n] += A[n, m] * c[m]
            end
        end
    end
    mx .+ C
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
                c[n] +=  A[i, m] * B[n, m]
            end
        end
        for n = 1:p
            for m = 1:q
                 O[i, n] += A[n, m] * c[m]
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
                c[n] +=  A[i, m] * B[n, m]
            end
        end
        for n = 1:p
            for m = 1:q
                 O[i, n] += A[n, m] * c[m]
            end
        end
    end
    O .+= C
end
