function Base.show(io::IO, obj::RBE)
    println(io, "Bioequivalence Linear Mixed Effect Model")
    println(io, "")
    println(io, "REML: ", obj.reml)
    println(io, "")
    lines = 0
    for f in obj.factors
        lines += length(obj.model.contrasts[f].termnames)
    end
    pm = Array{Any,2}(undef, lines+2, 5)
    pm[1,1] = "Level"; pm[1,2] = "Value"; pm[1,3] = "SE"; pm[1,4] = "DF"; pm[1,5] = "F";
    pm[2,1] = "Intercept"; pm[2,2] = obj.β[1]; pm[2,3] = "-"; pm[2,4] = "-"; pm[2,5] = "-";
    it      = 2
    pmr     = 3
    for f in obj.factors
        for l in obj.model.contrasts[f].termnames
            pm[pmr,1] = string(f)*": "*string(l);
            pm[pmr,2] = obj.β[it];
            pm[pmr,3] = obj.se[it];
            pm[pmr,4] = obj.df[it];
            pm[pmr,5] = obj.f[it];
            it  += 1
            pmr += 1
        end
    end
    for r = 1:size(pm)[1]
        for c = 1:size(pm)[2]
            if isa(pm[r,c], Float64) pm[r,c] = round(pm[r,c], sigdigits=6) end
        end
    end
    #─┼┴┬│
    pm    = string.(pm)
    len   = Array{Int,1}(undef, 0)
    vch   = Array{Vector{Char},1}(undef, 0)
    for c = 1:size(pm)[2]
        ml = maximum(length.(pm[:,c]))
        push!(len, ml)
        ch = Vector{Char}(undef, ml)
        ch .= '─'
        push!(vch, ch)
    end
    for v in vch
        print(io, String(v))
        print(io, "┬")
    end
    println(io, "")
    for c = 1:size(pm)[2]
        print(io, pm[1,c])
        ch  = Vector{Char}(undef, len[c]-length(pm[1,c]))
        ch .= ' '
        print(io, String(ch))
        print(io, "│")
    end
    println(io, "")
    for r = 2:size(pm)[1]
        for v in vch
            print(io, String(v))
            print(io, "┼")
        end
        println(io, "")
        for c = 1:size(pm)[2]
            print(io, pm[r,c])
            ch  = Vector{Char}(undef, len[c]-length(pm[r,c]))
            ch .= ' '
            print(io, String(ch))
            print(io, "│")
        end
        println(io, "")
    end
end
