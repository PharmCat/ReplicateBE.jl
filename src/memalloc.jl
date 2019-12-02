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
