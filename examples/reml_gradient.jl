function vmat(R, G, Z)
    V  = Z*G*Z' + R
    return V
end
function rmat(σ::Vector, Z::Matrix)
    return Matrix(Diagonal((Z*σ)))
end
function gmat(σ::Vector)
    cov = sqrt(σ[1] * σ[2]) * σ[3]
    return [σ[1] cov; cov σ[2]]
end
#-------------------------------------------------------------------------------
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
#-------------------------------------------------------------------------------
function gradvmat(V, vmatvecfunc, θvec, cache)
    h = hash(V)
    if h in keys(cache)
        return cache[h]
    else
        jV  = ForwardDiff.jacobian(vmatvecfunc, θvec)
        result = Vector{Matrix}(undef, 0)
        for i in 1:length(θvec)
            push!(result, reshape(jV[:,i], size(V,1), size(V,1)))
        end
        cache[h] = result
        return result
    end
end

vmatvec = x -> vmat(rmat(x[1:2], Zv[i]), gmat(x[3:end]), Zv[i])[:]
#-------------------------------------------------------------------------------


function greml(yv, Zv, p, Xv, θvec, β)
    n     = length(yv)
    G     = gmat(θvec[3:5])
    θ1    = zeros(length(θvec))
    θ2    = zeros(length(θvec))
    θ3    = zeros(length(θvec))
    iV    = nothing
    θ2m   = zeros(p, p)
    H     = zeros(p, p)
    cache = Dict()
    cache2 = Dict()
    for i = 1:n
        H += Xv[i]'*minv(vmat(rmat(θvec[1:2], Zv[i]), G, Zv[i]), cache)*Xv[i]
    end
    iH = inv(H)
    for i = 1:n
        #V matrix derivation vector
        vmatdvec = x -> vmat(rmat(x[1:2], Zv[i]), gmat(x[3:end]), Zv[i])[:]

        V   = vmat(rmat(θvec[1:2], Zv[i]), G, Zv[i])
        iV  = minv(V, cache)
        r   = yv[i] .- Xv[i]*β
        jV  = gradvmat(V, vmatdvec, θvec, cache2)
        #jV  = ForwardDiff.jacobian(vmatdvec, θvec)
        for i2 = 1:length(θvec)

            θ1[i2]  += tr(iV*jV[i2])

            A        = iV*jV[i2]*iV
            θ2[i2]  -= tr(iH*Xv[i]'*A*Xv[i])

            θ3[i2]  -= r'*A*r
        end
    end
    return - (θ1 .+ θ2 .+ θ3) ./ 2
end

gradremlf(x)  = -2*greml(yv, Zv, rank(X), Xv, x, β)
hessremlf(x)  = ForwardDiff.hessian(remlf, x)
function g!(G, x)
    copyto!(G, gradremlf(x))
end
function h!(H, x)
    copyto!(H, hessremlf(x))
end

O   = optimize(remlf, gradremlf, hessremlf, θvec0, method=Newton(); inplace = false)
