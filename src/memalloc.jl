struct MemCache
    svec::Vector{Vector}
    mvec::Vector{Symmetric} #vector of Symmetric matrices
    #dict::Dict
    function MemCache(maxobs)
        svec    = Vector{Vector}(undef, 1)
        svec[1] = zeros(1)
        mvec    = Vector{Symmetric}(undef, maxobs)
        new(svec, mvec)::MemCache
    end
end

function rebuildcache(data, type)
        if !(eltype(data.mem.svec[1]) <: type)
            data.mem.svec[1] = zeros(type, data.maxobs)
            for i = 1:data.maxobs
                data.mem.mvec[i] = Symmetric(zeros(type, i, i))
            end
        end
end
