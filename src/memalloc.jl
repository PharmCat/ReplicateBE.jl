struct MemCache
    svec::Vector{Vector}
    #dict::Dict
    function MemCache(maxobs)
        svec  = Vector{Vector}(undef, 1)
        new(svec, #=Dict{Matrix, Tuple{Matrix, Matrix, Number}}()=#)::MemCache
    end
end
