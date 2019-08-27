"""
    Return set of R matrices
"""
function rmatvec(σ₁, σ₂, Zvec)
    n  = length(Zvec)
    Ra = Array{Array{Float64,2}, 1}(undef, n)
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
    Va = Array{Array{Float64,2}, 1}(undef, n)
    for i = 1:length(Zvec)
        Va[i] = vmat(G, Rvec[i], Zvec[i])
    end
    return Va
end
