
B = [1.0 0.0; 0.0 1.0]
A = [1.0 2.0; 2.0 4.0]
D = [1.0 2.0; 2.0 4.0]
C = [1.0, 2.0]
c = [0.0, 0.0]
θ = [0.0 0.0; 0.0 0.0]
β = [0.0, 0.0]
g = [1.0, 2.0, 3.0]
σ = [1.0, 2.0, 1.0, 1.5, 0.5]
t = [1.0, 2.0, 1.0, 1.5, 0.5]
@code_warntype ReplicateBE.mulαtβαinc!(θ, A, B, c)

@code_warntype ReplicateBE.mulαtβαinc!(θ, A, B)

@code_warntype ReplicateBE.mulθβinc!(θ, β, A, B, C, c)

@code_warntype ReplicateBE.gmat(g)

@code_warntype ReplicateBE.vmat(A, B, D)

mem   = ReplicateBE.MemCache(2, Float64)
@code_warntype ReplicateBE.mulαβαtcupd!(mem.mvec, A, B, C, mem.svec)


cache = Dict{Matrix, Tuple{Symmetric{Float64, Matrix{Float64}}, Float64}}()
@code_warntype ReplicateBE.mvmatall(A, σ, B, mem, cache)

@code_warntype ReplicateBE.rholinkpsigmoid(3.0, 1.5)

@code_warntype ReplicateBE.varlinkmap(t, 1:2, 3:5, ReplicateBE.rholinkpsigmoid, ReplicateBE.rholinkpsigmoid)
