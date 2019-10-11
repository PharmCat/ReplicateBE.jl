
mutable struct RandRBEDS
    n::Int
    sequence::Vector
    design::Matrix
    inter::Vector
    intra::Vector
    intercept::Real
    seqcoef::Vector
    periodcoef::Vector
    formcoef::Vector
    dropsubj::Float64
    dropobs::Int
    dataset::DataFrame
end

"""
    randrbeds(;n=24, sequence=[1,1],
        design = ["T" "R" "T" "R"; "R" "T" "R" "T"],
        inter=[0.5, 0.4, 0.9], intra=[0.1, 0.2],
        intercept = 0, seqcoef = [0.0, 0.0], periodcoef = [0.0, 0.0, 0.0, 0.0], formcoef = [0.0, 0.0],
        dropsubj = 0.0, dropobs::Int = 0, seed::Int = 0)

Random dataset for bioequivalence.
"""
function randrbeds(;n=24, sequence=[1,1],
    design = ["T" "R" "T" "R"; "R" "T" "R" "T"],
    inter=[0.5, 0.4, 0.9], intra=[0.1, 0.2],
    intercept = 0, seqcoef = [0.0, 0.0], periodcoef = [0.0, 0.0, 0.0, 0.0], formcoef = [0.0, 0.0],
    dropsubj = 0.0, dropobs::Int = 0, seed::Int = 0)
    return randrbeds(n, sequence, design, inter, intra, intercept, seqcoef, periodcoef, formcoef, dropsubj, dropobs, seed)
end

"""
    randrbeds(n::Int, sequence::Vector,
        design::Matrix,
        θinter::Vector, θintra::Vector,
        intercept::Real, seqcoef::Vector, periodcoef::Vector, formcoef::Vector,
        dropsubj::Float64, dropobs::Int, seed::Int)

Random dataset for bioequivalence.

"""
function randrbeds(n::Int, sequence::Vector,
    design::Matrix,
    θinter::Vector, θintra::Vector,
    intercept::Real, seqcoef::Vector, periodcoef::Vector, formcoef::Vector,
    dropsubj::Float64, dropobs::Int, seed::Int)

    rng = MersenneTwister()
    if seed == 0  Random.seed!(rng) else Random.seed!(seed) end

    r = n/sum(sequence)
    sn = Array{Int, 1}(undef, length(sequence))
    for i = 1:(length(sequence)-1)
        sn[i] = round(r*sequence[i])
    end
    sn[length(sequence)] = n - sum(sn[1:(length(sequence)-1)])

    u      = unique(design)
    sqname = Array{String, 1}(undef,size(design)[1])
    sqnum  = size(design)[1]
    pnum   = size(design)[2]
    for i = 1:sqnum
        sqname[i] = join(design[i,:])
    end
    Zv = Array{Matrix, 1}(undef, sqnum)
    Vv = Array{Matrix, 1}(undef, sqnum)
    G = gmat(θinter[1], θinter[2], θinter[3])
    for i = 1:size(design)[1]
        Z = Array{Int, 2}(undef, pnum, length(u))
        for c = 1:pnum
            for uc = 1:length(u)
                if design[i, c] == u[uc] Z[c, uc] = 1 else Z[c, uc] = 0 end
            end
        end
        Zv[i] = Z
        Vv[i] = vmat(G, rmat(θintra, Z), Z)
    end
    Mv = Array{Array{Float64, 1}, 1}(undef, sqnum)
    for i = 1:sqnum
        Mv[i] = zeros(pnum) .+ intercept .+ seqcoef[i] + periodcoef + Zv[i]*formcoef
    end

    subjds = DataFrame(subject = String[], formulation = String[], period = String[], sequence = String[], var = Float64[])
    subj = 1
    subjmx = Array{Any, 2}(undef, pnum, 5)
    for i = 1:sqnum
        mvnorm = MvNormal(PDMat(Vv[i]))
        for sis = 1:sn[i]
            subjmx[:, 1] .= string(subj)
            subjmx[:, 2]  = string.(design[i,:])
            subjmx[:, 3]  = string.(collect(1:pnum))
            subjmx[:, 4] .= sqname[i]
            subjmx[:, 5]  = rand(MvNormal(PDMat(Vv[i]))) + Mv[i]
            subj += 1
            for c = 1:pnum
                push!(subjds, subjmx[c, :])
            end
        end
    end
    if dropobs > 0 && dropobs < size(subjds, 1)
        dellist = sample(1:size(subjds, 1), dropobs, replace = false)
        deleterows!(subjds, sort!(dellist))
    end
    return subjds
end
