struct MemCache
    svec::Vector{Vector}
    mvec::Vector{Matrix} 
    #dict::Dict
    function MemCache(maxobs)
        #ForwardDiff.Dual{ForwardDiff.Tag{ReplicateBE.var"#23#39"{Bool,Float64,ReplicateBE.RBEDataStructure},Float64}}
        svec    = Vector{Vector}(undef, 1)
        svec[1] = zeros(maxobs)
        mvec    = Vector{Matrix}(undef, 1)
        mvec[1] = zeros(maxobs, maxobs)
        new(svec, mvec)::MemCache
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
