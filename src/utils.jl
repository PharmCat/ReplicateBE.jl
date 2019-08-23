

function StatsBase.confint(obj::RBE, alpha::Float64; expci::Bool = false, inv::Bool = false)
    ifv = 1
    if inv ifv = -1 end
    a = Array{Tuple{Float64, Float64},1}(undef, length(obj.β)-1)
    for i = 2:length(obj.β)
        a[i-1] = calcci(obj.β[i]*ifv, obj.se[i], obj.DF[i], alpha, expci)
    end
    return Tuple(a)
end
function calcci(x::Float64, se::Float64, df::Float64, alpha::Float64, expci::Bool)::Tuple{Float64, Float64}
    q = quantile(TDist(df), 1.0-alpha/2)
    if expci
        return x-q*se, x+q*se
    else
        return exp(x-q*se), exp(x+q*se)
    end
end
function Base.show(io::IO, obj::Tuple{Vararg{Tuple{Float64, Float64}}})
    for i in obj
        println(i)
    end
end
