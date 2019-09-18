struct EffectTable
    name
    est
    se
    f
    df
    t
    p
    function EffectTable(name, est, se, f, df, t, p)
        new(name, est, se, f, df, t, p)
    end
end

function Base.getindex(t::EffectTable, r::Int, c::Int)
    return getfield(t, fieldnames(typeof(t))[c])[r]
end

function Base.show(io::IO, t::EffectTable)
    header      = ["Effect", "Value" , "SE",  "DF" , "F", "t", "P"]
    matrix      = Array{String, 2}(undef, length(t.name), length(header))
    matrix[:,1] = string.(t.name)
    for c = 2:length(header)
        for r = 1:length(t.name)
            matrix[r,c] = string(round(t[r,c], sigdigits=6))
        end
    end
    l = maximum(length.(matrix), dims = 1) .+ 3
    for c = 1:length(header)
        if l[c] < length(header[c]) l[c] = length(header[c]) + 3 end
        header[c] = addspace(header[c], l[c]-length(header[c]))
        for r = 1:length(t.name)
            matrix[r,c] = addspace(matrix[r,c], l[c]-length(matrix[r,c]))
        end
    end
    chl = Vector{Char}(undef, sum(l))
    chl .= '─'
    line = String(chl)
    println(io, line)
    println(io, header ...)
    println(io, line)
    for r = 1:length(t.name)
        println(io, matrix[r,:] ...)
    end
    print(io, line)
end
