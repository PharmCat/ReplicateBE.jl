struct MemCache{T}
    svec::Vector{T}
    mvec::Matrix{T}
    #dict::Dict
    function MemCache(maxobs, t)
        #ForwardDiff.Dual{ForwardDiff.Tag{ReplicateBE.var"#23#39"{Bool,Float64,ReplicateBE.RBEDataStructure},Float64}}
        #svec    = Vector{Vector{t}}(undef, 1)
        svec = zeros(t, maxobs)
        #mvec    = Vector{Matrix{t}}(undef, 1)
        mvec = zeros(t, maxobs, maxobs)
        new{t}(svec, mvec)::MemCache{t}
    end
end

#function Base.Float64(n::T) where T <: ForwardDiff.Dual
#    n.value
#end

function rebuildcache(data, type)
        if !(eltype(data.mem.svec[1]) <: type)
            data.mem.svec[1] = zeros(type, data.maxobs)
            data.mem.mvec[1] = zeros(type, data.maxobs, data.maxobs)
            #for i = 1:data.maxobs
            #    data.mem.mvec[i] = Symmetric(zeros(type, i, i))
            #end
        end
end
