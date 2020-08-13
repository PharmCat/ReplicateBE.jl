struct MemCache
    svec::Vector{Vector}
    mvec::Vector{Matrix}
    #dict::Dict
    function MemCache(maxobs)
        svec    = Vector{Vector}(undef, 1)
        svec[1] = zeros(1)
        mvec    = Vector{Matrix}(undef, maxobs)
        new(svec, mvec)::MemCache
    end
end
