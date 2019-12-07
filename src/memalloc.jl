struct MemCache
    svec::Vector{Vector}
    function MemCache(maxobs)
        svec  = Vector{Vector}(undef, 1)
        new(svec)::MemCache
    end
end
