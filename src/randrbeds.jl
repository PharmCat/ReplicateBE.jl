"""
Random dataset task.

```julia
mutable struct RandRBEDS
    n::Int                       # Subjects number
    sequence::Vector             # Sequence distribution vector
    design::Matrix               # Design matrix
    inter::Vector                # Intersubject variance part
    intra::Vector                # Intrasubject variance part
                                 # Fixed effect:
    intercept::Real              # Intercept
    seqcoef::Vector              # Sequence
    periodcoef::Vector           # Period
    formcoef::Vector             # Formulation
    dropobs::Int                 # Drop observations
    seed                         # Seed
end
```
"""
mutable struct RandRBEDS
    n::Int
    sequence::Vector
    design::Matrix
    inter::Union{Vector, Real}
    intra::Vector
    intercept::Real
    seqcoef::Vector
    periodcoef::Vector
    formcoef::Vector
    #dropsubj::Real                #Deprecated
    dropobs::Int
    seed
    function RandRBEDS(n::Int, sequence::Vector,
        design::Matrix,
        θinter, θintra::Vector,
        intercept::Real, seqcoef::Vector, periodcoef::Vector, formcoef::Vector,
        dropobs::Int, seed)
        new(n, sequence,
            design,
            θinter, θintra,
            intercept, seqcoef, periodcoef, formcoef,
            dropobs, seed)::RandRBEDS
    end
end

"""
```julia
randrbetask(;n=24, sequence=[1,1],
        design = ["T" "R" "T" "R"; "R" "T" "R" "T"],
        inter=[0.5, 0.4, 0.9], intra=[0.1, 0.2],
        intercept = 0, seqcoef = [0.0, 0.0], periodcoef = [0.0, 0.0, 0.0, 0.0], formcoef = [0.0, 0.0],
        dropsubj = 0, dropobs::Int = 0, seed = 0)::RandRBEDS
```

Make task for random dataset generation.
"""
function randrbetask(;n=24, sequence=[1,1],
        design = ["T" "R" "T" "R"; "R" "T" "R" "T"],
        inter=[0.5, 0.4, 0.9], intra=[0.1, 0.2],
        intercept = 0, seqcoef = [0.0, 0.0], periodcoef = [0.0, 0.0, 0.0, 0.0], formcoef = [0.0, 0.0],
        dropobs::Int = 0, seed = 0)::RandRBEDS

    return RandRBEDS(n, sequence,
            design,
            inter, intra,
            intercept, seqcoef, periodcoef, formcoef,
            dropobs, seed)
end

struct RBEDSSimResult
    seed
    num
    seeds
    result
    errn
end
"""
```julia
    randrbeds(;n=24, sequence=[1,1],
        design = ["T" "R" "T" "R"; "R" "T" "R" "T"],
        inter=[0.5, 0.4, 0.9], intra=[0.1, 0.2],
        intercept = 0, seqcoef = [0.0, 0.0], periodcoef = [0.0, 0.0, 0.0, 0.0],
        formcoef = [0.0, 0.0],
        dropobs::Int = 0, seed::Int = 0)
```

Random dataset generation for bioequivalence.

#  Keywords

* ``n=24`` - number of subjects;
* ``sequence = [1,1]`` -  distribution in sequences [1,1] means 1:1, [1,3] - 1:4 etc.;
* ``design = ["T" "R" "T" "R"; "R" "T" "R" "T"]`` - desin matrix, each line is a sequence, each column - periods, cell - formulation id;
* ``inter=[0.5, 0.4, 0.9]`` - inter-subject variance vector for G matrix (length 3): [σ₁, σ₂, ρ], where σ₁, σ₂ - formulation inter-subject variance,  ρ - covariance coefficient;
* ``intra=[0.1, 0.2]`` - intra-subject variance vector for R matrix (length 2): [σ₁, σ₂], where σ₁, σ₂ - formulation intra-subject variance for each formulation;
* ``intercept = 0`` - intercept effect value;
* ``seqcoef = [0.0, 0.0]`` - coefficients of sequences, additive (length(sequence) == length(seqcoef) == size(design, 1));
* ``periodcoef = [0.0, 0.0, 0.0, 0.0]`` - coefficients of periods, additive  (length(periodcoef) == size(design, 2));
* ``formcoef = [0.0, 0.0]`` -  coefficients of formulations, additive ;
* ``dropobs = 0`` number of randomly dropped observations;
* ``seed = 0`` - seed for random number generator (0 - random seed).



Multivariate normal disribution:

```math
f(\\mathbf{x}; \\boldsymbol{\\mu}, \\boldsymbol{V}) = \\frac{1}{(2 \\pi)^{d/2} |\\boldsymbol{V}|^{1/2}}
\\exp \\left( - \\frac{1}{2} (\\mathbf{x} - \\boldsymbol{\\mu})^T V^{-1} (\\mathbf{x} - \\boldsymbol{\\mu}) \\right)
```

Where V:

```math
V_{i} = Z_{i}GZ_i'+R_{i}
```

"""
function randrbeds(;n=24, sequence=[1,1],
    design = ["T" "R" "T" "R"; "R" "T" "R" "T"],
    inter=[0.5, 0.4, 0.9], intra=[0.1, 0.2],
    intercept = 0, seqcoef = [0.0, 0.0], periodcoef = [0.0, 0.0, 0.0, 0.0], formcoef = [0.0, 0.0],
    dropobs::Int = 0, seed = 0)
    return randrbeds(n, sequence, design, inter, intra, intercept, seqcoef, periodcoef, formcoef, dropobs, seed)
end

"""
Another way to use:

```julia
    randrbeds(n::Int, sequence::Vector,
        design::Matrix,
        θinter::Vector, θintra::Vector,
        intercept::Real, seqcoef::Vector, periodcoef::Vector, formcoef::Vector,
        dropsubj::Float64, dropobs::Int, seed::Int)
```
"""
function randrbeds(n::Int, sequence::Vector,
    design::Matrix,
    θinter::Vector, θintra::Vector,
    intercept::Real, seqcoef::Vector, periodcoef::Vector, formcoef::Vector,
    dropobs::Int, seed)
    if seed != 0
        rng = MersenneTwister(seed)
    else
        rng = MersenneTwister()
    end

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
    G = gmat(θinter)
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

    subjds = DataFrame(subject = Int[], formulation = String[], period = Int[], sequence = String[], var = Float64[])
    subj = 1
    subjmx = Array{Any, 2}(undef, pnum, 5)
    for i = 1:sqnum
        for sis = 1:sn[i]
            subjmx[:, 1] .= subj
            subjmx[:, 2]  = design[i,:]
            subjmx[:, 3]  = collect(1:pnum)
            subjmx[:, 4] .= sqname[i]
            subjmx[:, 5]  = rand(rng, MvNormal(PDMat(Vv[i]))) + Mv[i]
            subj += 1
            for c = 1:pnum
                push!(subjds, subjmx[c, :])
            end
        end
    end
    if dropobs > 0 && dropobs < size(subjds, 1)
        dellist = sample(rng, 1:size(subjds, 1), dropobs, replace = false)
        deleterows!(subjds, sort!(dellist))
    end
    categorical!(subjds, :subject);
    categorical!(subjds, :formulation);
    categorical!(subjds, :period);
    categorical!(subjds, :sequence);
    return subjds
end
"""
Using with RandRBEDS object:

```julia
randrbeds(task::RandRBEDS)
```
"""
function randrbeds(task::RandRBEDS)
    return randrbeds(task.n, task.sequence, task.design,
                    task.inter, task.intra,
                    task.intercept, task.seqcoef, task.periodcoef, task.formcoef,
                    task.dropobs, task.seed)
end

"""
```julia
simulation(task::RandRBEDS; io = stdout, verbose = false, num = 100,
    l = log(0.8), u = log(1.25), seed = 0)
```

Count successful BE outcomes.

# Parameters

* task -  RandRBEDS object

# Keywords

* ``io = stdout`` - text output
* ``verbose = false`` - print messages to io
* ``num = 100`` - number of simulations
* ``l = log(0.8)`` - lower bound
* ``u = log(1.25)`` - upper bound
* ``seed = 0`` - seed for random number generator (0 - random seed)

"""
function simulation(task::RandRBEDS; io = stdout, verbose = false, num = 100, l = log(0.8), u = log(1.25), seed = 0, rsabe = false, rsabeconst = 0.760, reference = "R", alpha = 0.05)
    tl = l
    tu = u
    if rsabe
        if reference ∉ task.design
            rsabe = false
            @warn "Reference value not found in design. RSABE set false."
        end
    end

    #max range 69,84-143,19
    task.seed = 0
    #rng = MersenneTwister()
    if isa(seed, Array)
        seeds = seed
    else
        if seed != 0
            rng = MersenneTwister(seed)
        else
            rng = MersenneTwister()
        end
        seeds = Array{UInt32, 1}(undef, num)
        for i = 1:num
            seeds[i] = rand(rng, UInt32)
        end
    end

    n     = 0
    err   = 0
    cnt   = 0
    if verbose
        printstyled(io, "Start...\n"; color = :green)
        if isa(seed, Array)
            println(io, "Simulation seed: Array")
        else
        println(io, "Simulation seed: $(seed)")
        end
        println(io, "Task hash: $(hash(task))")
        println(io, "Alpha: $(alpha)")
        println(io, "RSABE: $(rsabe)")
        if rsabe
            println(io, "Regulatory const: $(rsabeconst)")
            println(io, "Reference formulation: $(reference)")
        end
    end

    for i = 1:num
        task.seed = seeds[i]
        rds       = randrbeds(task)
        try

            be        = rbe(rds, dvar = :var, subject = :subject, formulation = :formulation, period = :period, sequence = :sequence)
            #q         = quantile(TDist(be.fixed.df[end]), 1.0 - alpha)
            #ll        = be.fixed.est[end] - q*be.fixed.se[end]
            #ul        = be.fixed.est[end] + q*be.fixed.se[end]
            ll, ul    = confint(be, 2 * alpha)[end]
            #println("CI: $(exp(ll)) - $(exp(ul)) ")
            #!
            if verbose
                if !optstat(be) printstyled(io, "Iteration: ", i, ", seed ", seeds[i], ": unconverged! \n"; color = :yellow) end
                if !isposdef(Symmetric(be.H))
                    printstyled(io, "Iteration: ", i, ", seed ", seeds[i], ": Hessian not positive defined! \n"; color = :yellow)
                end
            end
            # If RSABE is true - calculating CI limits
            if rsabe
                ivar = intravar(be)[reference]
                if geocv(ivar) > 0.30
                    bconst = rsabeconst * sqrt(ivar)
                    tl  = -bconst
                    tu  = bconst
                    if tu > log(1.4319) || tl < log(0.6984)
                        tu = log(1.4319)
                        tl = log(0.6984)
                    end
                    if verbose
                        println("Reference scaled: $(exp(tl)*100) - $(exp(tu)*100)")
                    end
                else
                    tl  = l
                    tu  = u
                end

            end
            if ll > tl && ul < tu
                cnt += 1
            end
            #!
            if n > 1000
                println(io, "Iteration: $i")
                println(io, "Mem: $(Sys.free_memory()/2^20)")
                println(io, "Pow: $(cnt/i)")
                println(io, "-------------------------------")
                n = 0
            end
            n += 1

        catch
            err += 1
            printstyled(io, "Iteration: $i, seed $(seeds[i]): $(err): ERROR! \n"; color = :red)
        end

    end
    return RBEDSSimResult(seed, num, seeds, cnt/(num - err), err)
end



"""
```julia
simulation!(task::RandRBEDS, out, simfunc!::Function; io = stdout, verbose = false, num = 100, seed = 0)
```

Generalized simulation method.

"""
function simulation!(task::RandRBEDS, out, simfunc!::Function; io = stdout, verbose = false, num = 100, seed = 0)
    task.seed = 0
    if isa(seed, Array)
        seeds = seed
    else
        if seed != 0
            rng = MersenneTwister(seed)
        else
            rng = MersenneTwister()
        end
        seeds = Array{UInt32, 1}(undef, num)
        for i = 1:num
            seeds[i] = rand(rng, UInt32)
        end
    end
    n     = 0
    err   = 0
    if verbose
        printstyled(io, "Custom simulation start...\n"; color = :green)
        if isa(seed, Array)
            println(io, "Simulation seed: Array")
        else
            println(io, "Simulation seed: $(seed)")
        end
        println(io, "Task hash: $(hash(task))")
    end
    for i = 1:num
        task.seed = seeds[i]
        rds       = randrbeds(task)
        try
            be        = rbe(rds, dvar = :var, subject = :subject, formulation = :formulation, period = :period, sequence = :sequence)
            simfunc!(out, be)
            if n > 1000
                println(io, "Iteration: $i")
                println(io, "Mem: $(Sys.free_memory()/2^20)")
                println(io, "Pow: $(cnt/i)")
                println(io, "-------------------------------")
                n = 0
            end
            n += 1
        catch
            err += 1
            printstyled(io, "Iteration: $i, seed $(seeds[i]): $(err): ERROR! \n"; color = :red)
        end
    end
    return out
end

function Base.show(io::IO, obj::RBEDSSimResult)
    if isa(obj.seed, Array)
        println(io, "Simulation seed: Array")
    else
        println(io, "Seed: $(obj.seed)")
    end
    println(io, "Seed: $(obj.seed)")
    println(io, "Number: $(obj.num)")
    println(io, "Errors: $(obj.errn)")
    println(io, "Result: $(obj.result)")
end


function randrbeds(n::Int, sequence::Vector,
    design::Matrix,
    inter::Real, intra::Vector,
    intercept::Real, seqcoef::Vector, periodcoef::Vector, formcoef::Vector,
    dropobs::Int, seed)
    if seed != 0
        rng = MersenneTwister(seed)
    else
        rng = MersenneTwister()
    end

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
    Vv = Array{Vector, 1}(undef, sqnum)
    for i = 1:size(design)[1]
        Z = Array{Int, 2}(undef, pnum, length(u))
        for c = 1:pnum
            for uc = 1:length(u)
                if design[i, c] == u[uc] Z[c, uc] = 1 else Z[c, uc] = 0 end
            end
        end
        Zv[i] = Z
        Vv[i] = Z * intra
    end
    Mv = Array{Array{Float64, 1}, 1}(undef, sqnum)
    for i = 1:sqnum
        Mv[i] = zeros(pnum) .+ intercept .+ seqcoef[i] + periodcoef + Zv[i]*formcoef
    end
    ndist  = Normal()
    subjds = DataFrame(subject = Int[], formulation = String[], period = Int[], sequence = String[], var = Float64[])
    subj   = 1
    subjmx = Array{Any, 2}(undef, pnum, 5)
    for i = 1:sqnum
        for sis = 1:sn[i]
            subjmx[:, 1] .= subj
            subjmx[:, 2]  = design[i,:]
            subjmx[:, 3]  = collect(1:pnum)
            subjmx[:, 4] .= sqname[i]
            subjmx[:, 5] .= 0
            subjmx[:, 5] .+= rand(rng, ndist)*sqrt(inter)
            subj += 1
            for c = 1:pnum
                subjmx[c, 5] += Mv[i][c] + rand(rng, ndist)*sqrt(Vv[i][c])
                push!(subjds, subjmx[c, :])
            end
        end
    end
    if dropobs > 0 && dropobs < size(subjds, 1)
        dellist = sample(rng, 1:size(subjds, 1), dropobs, replace = false)
        deleterows!(subjds, sort!(dellist))
    end
    categorical!(subjds, :subject);
    categorical!(subjds, :formulation);
    categorical!(subjds, :period);
    categorical!(subjds, :sequence);
    return subjds
end
