struct MemAlloc{T <: AbstractFloat}
    mem1::Array{Array{T, 2}, 1}
    mem2::Array{Array{T, 2}, 1}
    mem3::Array{Array{T, 1}, 1}
    #mem4::Array{Float64, 2}
    function MemAlloc(p, zs, yv::Vector{Vector{T}}) where T <: AbstractFloat
        maxobs  = maximum(length.(yv))
        yvtype = eltype(yv[1])
        memc1 = Array{Array{yvtype, 2}, 1}(undef, maxobs)
        memc2 = Array{Array{yvtype, 2}, 1}(undef, maxobs)
        memc3 = Array{Array{yvtype, 1}, 1}(undef, maxobs)
        for i = 1:maxobs
            memc1[i] = zeros(i, zs)
            memc2[i] = zeros(p, i)
            memc3[i] = zeros(i)
        end
        #memc4 = zeros(p, p)
        new{T}(memc1, memc2, memc3)::MemAlloc
    end
end
