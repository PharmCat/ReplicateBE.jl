using DataFrames, CSV, GLM, StatsModels, Random, SparseArrays, ForwardDiff, LinearAlgebra, DiffResults, Optim, BenchmarkTools #, NLopt
filep = ".\\Mixed\\mbe.csv"
if isfile(filep) df = CSV.File(filep, delim=',') |> DataFrame end
categorical!(df, :subject);
categorical!(df, :sequence);
categorical!(df, :period);
categorical!(df, :formulation);
sort!(df,[:subject, :formulation, :period])
X = ModelMatrix(ModelFrame(@formula(var ~ sequence + period + formulation), df)).m
#Xf = ModelMatrix(ModelFrame(@formula(response ~ repeat + effect), df)).m
Z = ModelMatrix(ModelFrame(@formula(var ~ 0+formulation), df, contrasts = Dict(:effect => StatsModels.FullDummyCoding()))).m
#X = X[1:4,:]
y = df[:, :var]

#----

function getsubjects(df, sbj, X, Z, y)
    u = unique(df[:, sbj])
    Xa = Array{Array{Float64,2}, 1}(undef, 0)
    Za = Array{Array{Float64,2}, 1}(undef, 0)
    ya = Array{Array{Float64,1}, 1}(undef, 0)
    for i in u
        v = findall(x->x==i, df[:, sbj])
        Xs = Array{Float64, 1}(undef, 0)
        Zs = Array{Float64, 1}(undef, 0)
        ys = Array{Float64, 1}(undef, 0)
        for r in v
            append!(Xs, X[r, :])
            append!(Zs, Z[r, :])
            push!(ys, y[r])
        end
        push!(Xa, Matrix(reshape(Xs, size(X)[2], :)'))
        push!(Za, Matrix(reshape(Zs, size(Z)[2], :)'))
        push!(ya, ys)
    end
    return Xa, Za, ya
end
function rmatvec(σ₁, σ₂, Zvec)
    n  = length(Zvec)
    Ra = Array{Array{Float64,2}, 1}(undef, n)
    for i = 1:n
        Ra[i] = rmat([σ₁, σ₂], Zvec[i])
    end
    return Ra
end
function vmatvec(Zvec, G, Rvec)
    Va = Array{Array{Float64,2}, 1}(undef, 0)
    for i = 1:length(Zvec)
        push!(Va, cov(G, Rvec[i], Zvec[i]))
    end
    return Va
end
#----

function gmat(σ₁, σ₂, ρ)
    if ρ > 1.0 ρ = 1.0 end
    if σ₁ < 0.0 σ₁ = 1.0e-6 end
    if σ₂ < 0.0 σ₂ = 1.0e-6 end
    cov = sqrt(σ₁ * σ₂) * ρ
    return [σ₁ cov; cov σ₂]
end

function rmat(σ, Z)
    if σ[1] < 0.0 σ[1] = 1.0e-6 end
    if σ[2] < 0.0 σ[2] = 1.0e-6 end
    return Matrix(Diagonal((Z*σ)[:,1]))
end

function rmat3(σ₁, σ₂, Z)
    if σ₁ < 0.0 σ₁ = 1.0e-6 end
    if σ₂ < 0.0 σ₂ = 1.0e-6 end
    M = zeros(4, 4)
    for i = 1:4
        M[i, i] = Z[i,:]'*[σ₁, σ₂]
    end
    return M
end

function rmat2(σ₁, σ₂, Z)
    if σ₁ < 0.0 σ₁ = 1.0e-6 end
    if σ₂ < 0.0 σ₂ = 1.0e-6 end
    return [σ₁ 0 0 0; 0 σ₁ 0 0; 0 0 σ₂ 0; 0 0 0 σ₂]
end

function cov(G, R, Z)
    V  = Z*G*Z' + R
end

function cmat(Xv, Zv, θ)
    p = size(Xv[1])[2]
    C = zeros(p,p)
    G   = gmat(θ[3], θ[4], θ[5])
    Rv  = rmatvec(θ[1], θ[2], Zv)
    Vv  = vmatvec(Zv, G, Rv)
    iVv = inv.(Vv)
    for i=1:length(Xv)
        C = C + Xv[i]'*iVv[i]*Xv[i]
    end
    return pinv(C)
end


function βcoef(y, X, iV)
    p = rank(X)
    A = zeros(p,p)
    β = zeros(p)
    for i = 1:4
        s = ((i-1)*4 + 1)
        e = ((i-1)*4 + 4)
        A = A + (X[s:e,:]'*iV*X[s:e,:])
        β = β .+ X[s:e,:]'*iV*y[s:e]
    end
    return inv(A)*β
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

function βcoef!(yv, X, Xv, iVv, β)
    p = rank(X)
    n = length(yv)
    A = zeros(p,p)
    β0 = zeros(p)
    for i = 1:n
        A = A + (Xv[i]'*iVv[i]*Xv[i])
        β0 = β0 .+ Xv[i]'*iVv[i]*yv[i]
    end
    β0 = inv(A)*β0
    copyto!(β, β0)
    return nothing
    #β = inv(A)*β0
end

function reml(yv, Zv, X, Xv, θvec)
    n = length(yv)
    N = sum(length.(yv))
    G = gmat(θvec[3], θvec[4], θvec[5])
    p = rank(X)
    c  = (N-p)/2*log(2π)
    θ1 = 0
    θ2 = 0
    θ3 = 0
    θ2m  = zeros(p,p)
    Rv   = rmatvec(θvec[1], θvec[2], Zv)
    Vv   = vmatvec(Zv, G, Rv)
    iVv  = inv.(Vv)
    β    = βcoef(yv, X, Xv, iVv)
    for i = 1:n
        θ1  += log(det(Vv[i]))
        θ2m += Xv[i]'*iVv[i]*Xv[i]
        r    = yv[i]-Xv[i]*β
        θ3  += r'*iVv[i]*r
    end
    θ2       = log(det(θ2m))
    #println("θ₁: ", θ1, " θ₂: ",  θ2,  " θ₃: ", θ3)
    return   -(θ1/2 + θ2/2 + θ3/2 + c)
end


function reml!(yv, Zv, X, Xv, G, Rv, Vv, iVv, θvec, β)
    n = length(yv)
    N = sum(length.(yv))
    G = gmat(θvec[3], θvec[4], θvec[5])
    p = rank(X)
    c  = (N-p)/2*log(2π)
    θ1 = 0
    θ2 = 0
    θ3 = 0
    θ2m  = zeros(p,p)
    Rv   = rmatvec(θvec[1], θvec[2], Zv)
    Vv   = vmatvec(Zv, G, Rv)
    copyto!(iVv, inv.(Vv))
    for i = 1:n
        θ1  += log(det(Vv[i]))
        θ2m += Xv[i]'*iVv[i]*Xv[i]
        r    = yv[i]-Xv[i]*β
        θ3  += r'*iVv[i]*r
    end
    θ2       = log(det(θ2m))
    #println("θ₁: ", θ1, " θ₂: ",  θ2,  " θ₃: ", θ3)
    return   -(θ1/2 + θ2/2 + θ3/2 + c)
end

#=
function reml(y, Z, X, θvec)
    n = 4
    N = length(y)
    G = gmat(θvec[3], θvec[4], θvec[5])
    p = rank(X)
    c  = (N-p)/2*log(2π)
    θ1 = 0
    θ2 = 0
    θ3 = 0

    θ2m = zeros(p,p)

    R    = rmat(θvec[1], θvec[2], Z[1:4,:])
    V    = Z[1:4,:]*G*Z[1:4,:]' + R
    iV   = inv(V)
    β    = βcoef(y, X, iV)

    for i = 1:n
        s    = ((i-1)*4 + 1)
        e    = ((i-1)*4 + 4)

        R    = rmat(θvec[1], θvec[2], Z[s:e,:])
        V    = Z[s:e,:]*G*Z[s:e,:]' + R
        iV   = inv(V)

        θ1  += log(det(V))

        θ2m += X[s:e,:]'*iV*X[s:e,:]
        r    = y[s:e]-X[s:e,:]*β
        θ3  += r'*iV*r
    end
    θ2 = log(det(θ2m))
    println("θ₁: ", θ1, " θ₂: ",  θ2,  " θ₃: ", θ3)
    return -(θ1/2 + θ2/2 + θ3/2 + c)
end
=#
θvec = [0.01089, 0.03604, 0.3003, 0.2199, 1.0]
θvec0 = [0.02, 0.02, 0.2, 0.2, 0.001]
β = [1.4156, 0.2, -0.116, 0.118, 0.1304, 0.025]
β0 = [1.0, 0.2, -0.1, 0.1, 0.1, 0.02]



#=
R = rmat(0.01089, 0.03604, Z[1:4,:])
G = gmat(0.3003, 0.2199, 1.0)
V = cov(G, R, Z[1:4,:])
=#

Xv, Zv, yv = getsubjects(df, :subject, X, Z, y)
Rv   = rmatvec(θvec0[1], θvec0[2], Zv)
G    = gmat(θvec0[3], θvec0[4], θvec0[5])
Vv   = vmatvec(Zv, G, Rv)
iVv  = Array{Array{Float64,2}, 1}(undef, length(Xv))
#iVv  = inv.(Vv)
remlf!(x) = -2*reml!(yv, Zv, X, Xv, G, Rv, Vv, iVv, x, β)
#-2*reml(y, Z, X, θvec)

remlf(x) = -2*reml(yv, Zv, X, Xv, x)

#remlf(x)  = -2*reml(y, Z, X, x)
function ppp()
    println("!")
end

#O = optimize(remlf2, θvec0, method=Newton(),  g_tol= 1e-10, x_tol = 1e-10, store_trace = true, extended_trace = true, show_trace = false)

O = optimize(remlf!, θvec0, method=Newton(), g_tol= 1e-10, x_tol = 1e-10, store_trace = true, extended_trace = true, show_trace = false)

#callback = βcoef!(yv, X, Xv, iVv, β)

θ = Optim.minimizer(O)
H = Optim.trace(O)[end].metadata["h(x)"]
#h = ForwardDiff.hessian(remlf, θvec)

G    = gmat(θ[3], θ[4], θ[5])
Rv   = rmatvec(θ[1], θ[2], Zv)
Vv   = vmatvec(Zv, G, Rv)
iVv  = inv.(Vv)
β    = βcoef(yv, X, Xv, iVv)

L    = [0 0 0 0 0 0]
L    = [0 0 0 0 0 1]

C    = cmat(Xv, Zv, θ)
SEL6 = sqrt((L*C*L')[1])
F = (L*β)'*inv(L*C*L')*(L*β)



#lclg(x) = det(gmat(x[3], x[4], x[5]))
#lclg(x) = det(rmat(x[1], x[2], Zv[1]))
#det(sum(rmatvec(θ[1], θ[2], Zv)))
#lclg(x) = (L*cmat(Xv, Zv, x)*L')[1]

function rmat3(σ, Z)
    if σ[1] < 0.0 σ[1] = 1.0e-6 end
    if σ[2] < 0.0 σ[2] = 1.0e-6 end
    return Matrix(Diagonal((Z*σ)[:,1]))
end

function lcgf(L, Xv, Zv, θ)
    p = size(Xv[1])[2]
    C = zeros(p,p)
    G   = gmat(θ[3], θ[4], θ[5])
    for i=1:length(Xv)
        R   = rmat3([θ[1], θ[2]], Zv[i])
        iVv = inv(cov(G, R, Zv[i]))
        C = C + Xv[i]'*iVv*Xv[i]
    end
    return (L*inv(C)*L')[1]
end
lclg(x) = lcgf(L, Xv, Zv, x)
g       = ForwardDiff.gradient(lclg, θ)
v       = 2*((L*C*L')[1])^2/(g'*(2*pinv(H))*g)

#-------------------------------------------------------------------------------

function βcoef2(yv, p, Xv, Zv, θvec)
    n = length(yv)
    A = zeros(p,p)
    β = zeros(p)
    G   = gmat(θvec[3], θvec[4], θvec[5])
    for i = 1:n
        R   = rmat3([θvec[1], θvec[2]], Zv[i])
        iV = inv(cov(G, R, Zv[i]))
        A = A + (Xv[i]'*iV*Xv[i])
        β = β .+ Xv[i]'*iV*yv[i]
    end
    return inv(A)*β
end
#=
function reml2(yv, Zv, p, Xv, θvec)
    n = length(yv)
    N = sum(length.(yv))
    G = gmat(θvec[3], θvec[4], θvec[5])
    c  = (N-p)*log(2π)
    θ1 = 0
    θ2 = 0
    θ3 = 0
    iV   = nothing
    θ2m  = zeros(p,p)
    βm   = zeros(p)
    iVv = Array{Array{Float64,2}, 1}(undef, n)

    for i = 1:n
        R   = rmat([θvec[1], θvec[2]], Zv[i])
        V   = cov(G, R, Zv[i])
        iV  = inv(V)
        iVv[i] = iV
        θ1  += logdet(V)
        tm   = Xv[i]'*iV    #Temp matrix for Xv[i]'*iV*Xv[i] and Xv[i]'*iV*yv[i] calc
        θ2m += tm*Xv[i]
        βm  += tm*yv[i]
    end
    β        = inv(θ2m)*βm
    for i = 1:n
        r    = yv[i]-Xv[i]*β
        θ3  += r'*iVv[i]*r
    end
    θ2       = logdet(θ2m)
    return   -(θ1 + θ2 + θ3 + c)
end
=#
mem = zeros(6,4)
function reml2(yv, Zv, p, Xv, θvec, β)
    n = length(yv)
    N = sum(length.(yv))
    G = gmat(θvec[3], θvec[4], θvec[5])
    c  = (N-p)/2*log(2π)
    θ1 = 0
    θ2 = 0
    θ3 = 0
    iV   = nothing
    θ2m  = zeros(p,p)
    mem = zeros(6,4)
    for i = 1:n
        R   = rmat([θvec[1], θvec[2]], Zv[i])
        V   = cov(G, R, Zv[i])
        iV  = inv(V)
        θ1  += logdet(V)
        mul!(mem, Xv[i]', iV)
        θ2m = θ2m .+ mem*Xv[i]
        #θ2m = θ2m .+ Xv[i]'*iV *Xv[i]
        r    = yv[i] .- Xv[i]*β
        θ3  += r'*iV*r
    end
    θ2       = logdet(θ2m)
    return   -(θ1/2 + θ2/2 + θ3/2 + c)
end

remlf2(x) = -reml2(yv, Zv, 6, Xv, x, β)
θ[5] = 1
cfg1 = ForwardDiff.HessianConfig(remlf2, θ, ForwardDiff.Chunk{5}());
h2 = ForwardDiff.hessian(remlf2, θ, cfg1)
h1 = ForwardDiff.hessian(remlf2, θ)

hr = DiffResults.HessianResult(θ)
hr = ForwardDiff.hessian!(hr, remlf2, θ)
h3 = DiffResults.hessian(hr)

gr(x) =  ForwardDiff.gradient(remlf2, x)
function g!(storage, x)
    copyto!(storage, gr(x))
end
he(x) =  ForwardDiff.hessian(remlf2, x)
function h!(storage, x)
    copyto!(storage, he(x))
end
st = zeros(5)
st2 = zeros(5,5)
g!(st, θ)
h!(st2, θ)


td = TwiceDifferentiable(remlf!, θvec0; autodiff = :forward)
Ooo =  optimize(td, θvec0, Newton())

Ooo = optimize(remlf!,  g!, h!, θvec0, BFGS())

h[5,:] .= 0
h[:,5] .= 0
L    = [0 0 1 0 0 0; 0 0 0 1 0 0; 0 0 0 0 1 0]
L    = [0 0 0 0 0 1]
v       = 2*((L*C*L')[1])^2/(g'*(2*pinv(h))*g)
