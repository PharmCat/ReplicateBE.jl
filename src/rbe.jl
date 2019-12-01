#
# GENARAL REPLICATE BIOEQUIVALENCE STRUCTURE
#
"""
```julia
    struct RBE
        model::ModelFrame               #Model frame
        rmodel::ModelFrame              #Random effect model
        design::Design                  #Design description
        factors::Array{Symbol, 1}       #Factor list
        θ0::Array{Float64, 1}           #Initial variance paramethers
        vlm::Real
        θ::Tuple{Vararg{Float64}}       #Final variance paramethers
        reml::Float64                   #-2REML
        fixed::EffectTable              #Fixed Effect table
        typeiii::ContrastTable          #Type III table
        R::Array{Matrix{Float64},1}     #R matrices for each subject
        V::Array{Matrix{Float64},1}     #V matrices for each subject
        G::Matrix{Float64}              #G matrix
        C::Matrix{Float64}              #C var(β) p×p variance-covariance matrix
        A::Matrix{Float64}              #asymptotic variance-covariance matrix of b θ
        H::Matrix{Float64}              #Hessian matrix
        X::Matrix                       #Matrix for fixed effects
        Z::Matrix                       #Matrix for random effects
        Xv::Array{Matrix{Float64},1}    #X matrices for each subject
        Zv::Array{Matrix{Float64},1}    #Z matrices for each subject
        yv::Array{Array{Float64, 1},1}  #responce vectors for each subject
        detH::Float64                   #Hessian determinant
        preoptim::Union{Optim.MultivariateOptimizationResults, Nothing}         #Pre-optimization result object
        optim::Optim.MultivariateOptimizationResults                            #Optimization result object
    end
```

Replicate bioequivalence structure.

"""
struct RBE
    model::ModelFrame               #Model frame
    rmodel::ModelFrame              #Random effect model
    design::Design
    factors::Array{Symbol, 1}       #Factor list
    θ0::Array{Float64, 1}           #Initial variance paramethers
    vlm::Real
    θ::Tuple{Vararg{Float64}}       #Final variance paramethers
    reml::Float64                   #-2REML
    fixed::EffectTable
    typeiii::ContrastTable
    R::Array{Matrix{Float64},1}     #R matrices for each subject
    V::Array{Matrix{Float64},1}     #V matrices for each subject
    G::Matrix{Float64}              #G matrix
    C::Matrix{Float64}              #C var(β) p×p variance-covariance matrix
    A::Matrix{Float64}              #asymptotic variance-covariance matrix ofb θ
    H::Matrix{Float64}              #Hessian matrix
    X::Matrix                       #Matrix for fixed effects
    Z::Matrix                       #Matrix for random effects
    Xv::Array{Matrix{Float64},1}    #X matrices for each subject
    Zv::Array{Matrix{Float64},1}    #Z matrices for each subject
    yv::Array{Array{Float64, 1},1}  #responce vectors for each subject
    detH::Float64                   #Hessian determinant
    preoptim::Union{Optim.MultivariateOptimizationResults, Nothing}        #Pre-optimization result object
    optim::Optim.MultivariateOptimizationResults           #Optimization result object
end

"""
```julia
rbe(df; dvar::Symbol,
    subject::Symbol,
    formulation::Symbol,
    period::Symbol,
    sequence::Symbol,
    g_tol::Float64 = 1e-8, x_tol::Float64 = 0.0, f_tol::Float64 = 0.0,
    iterations::Int = 100,
    store_trace = false, extended_trace = false, show_trace = false,
    memopt = true,
    init = [],
    postopt = false, maxopttry = 100)

```
Mixed model fitting function for replicate bioequivalence without data preparation (apply categorical! for each factor and sort! to dataframe).

Mixed model in matrix form:

```math
y = X\\beta + Zu + \\epsilon
```

with covariance matrix for each subject:

```math
V_{i} = Z_{i}GZ_i'+R_{i}
```

# Arguments

* df - DataFrame with data;

# Keywords:

* dvar - dependent variable;
* subject - subject factor;
* formulation - formulation factor;
* period - period factor;
* sequence - sequence factor;
* g_tol = 1e-8
* x_tol = 0.0
* f_tol = 0.0
* iterations = 100 - maximum iteration for optimization
* store_trace = false
* extended_trace = false
* show_trace = false
* memopt = true - memory optimization (function cache)
* init = [] - initial variance paremeters
* postopt = false - post optimization
* maxopttry = 100 - maximum attempts to optimize

"""
function rbe(df; dvar::Symbol,
    subject::Symbol,
    formulation::Symbol,
    period::Symbol,
    sequence::Symbol,
    g_tol::Real = 1e-12, x_tol::Real = 0.0, f_tol::Real = 0.0, iterations::Int = 100,
    store_trace = false, extended_trace = false, show_trace = false,
    memopt = true,
    init = [],
    postopt = false, vlm = 1.0, maxopttry = 100, rhoadjstep = 0.15,
    rholink = :psigmoid,
    singlim = 1e-10)

    #Check
    if any(x -> x ∉ names(df), [subject, formulation, period, sequence]) throw(ArgumentError("Names not found in DataFrame!")) end
    if !(eltype(df[!,dvar]) <: AbstractFloat)
        @warn "Responce variable ∉ Float!"
    end
    vartype = eltype(df[!,dvar])
    if !(typeof(df[!,subject]) <: CategoricalArray)
        @warn "Subject variable not Categorical, use rbe!()!"
    end
    if !(typeof(df[!,formulation]) <: CategoricalArray)
        @warn "Formulation variable not Categorical, use rbe!()!"
    end
    if !(typeof(df[!,period]) <: CategoricalArray)
        @warn "Period variable not Categorical, use rbe!()!"
    end
    if !(typeof(df[!,sequence]) <: CategoricalArray)
        @warn "Sequence variable not Categorical, use rbe!()!"
    end

    #Model
    Xf  = @eval(@formula($dvar ~ $sequence + $period + $formulation))
    Zf  = @eval(@formula($dvar ~ 0 + $formulation))
    MF  = ModelFrame(Xf, df)
    RMF = ModelFrame(Zf, df, contrasts = Dict(formulation => StatsModels.FullDummyCoding()))
    MM  = ModelMatrix(MF)
    X   = MM.m
    Z   = ModelMatrix(RMF).m
    p   = rank(X)
    zxr = rank(ModelMatrix(ModelFrame(@eval(@formula($dvar ~ $sequence + $period + $subject*$formulation)), df)).m)
    y   = df[:, dvar]                                                           #Dependent variable
    #Make pre located arrays with matrices for each subject
    Xv, Zv, yv = sortsubjects(df, subject, X, Z, y)
    n  = length(Xv)
    N  = sum(length.(yv))
    pn = termmodellen(MF, period)
    sn = termmodellen(MF, sequence)
    fn = termmodellen(MF, formulation)
    numdf = Dict(period => pn, sequence => sn, formulation => fn)
    sbf = Array{Int, 1}(undef, 0)
    fl  = MF.f.rhs.terms[findterm(MF, formulation)].contrasts.levels
    for i = 1:length(fl)
        push!(sbf, sbjnbyf(df, subject, formulation, fl[i]))
    end
    #Memory pre-allocation arrays for matrix computations
    memalloc = MemAlloc(p, 2, yv)
    #Check data
    checkdata(X, Z, Xv, Zv, y)
    #Calculate initial fixed parameters
    qro   = qr(X)
    β     = inv(qro.R)*qro.Q'*y
    #Calculate initial variance
    if length(init) == 5
        θvec0 = init
    else
        iv = initvar(df, dvar, formulation, subject)
        iv = iv .+ eps()
        θvec0 = [iv[3], iv[3], iv[1], iv[2], 0.5]
    end

    #Variance link function
    if rholink == :psigmoid
        varlink  = (x, y) ->  varlinkmap(x, 1:4, 5,  vlink,  z -> rholinkpsigmoid(z, y))
        rvarlink = (x, y) ->  varlinkmap(x, 1:4, 5,  vlinkr, z -> rholinkpsigmoidr(z, y))
    elseif rholink == :sigmoid
        varlink  = (x, y) ->  varlinkmap(x, 1:4, 5,  vlink,  z -> rholinksigmoid(z, y))
        rvarlink = (x, y) ->  varlinkmap(x, 1:4, 5,  vlinkr, z -> rholinksigmoidr(z, y))
    elseif rholink == :arctgsigmoid
        varlink  = (x, y) ->  varlinkmap(x, 1:4, 5,  vlink,  z -> rholinksigmoid2(z, y))
        rvarlink = (x, y) ->  varlinkmap(x, 1:4, 5,  vlinkr, z -> rholinksigmoidr2(z, y))
    else
        varlink  = (x, y) ->  varlinkmap(x, 1:4, 5,  vlink,  z -> rholinkpsigmoid(z, y))
        rvarlink = (x, y) ->  varlinkmap(x, 1:4, 5,  vlinkr, z -> rholinkpsigmoidr(z, y))
    end


    θvec0 = rvarlink(θvec0, vlm)
    #Prelocatiom for G, R, V, V⁻¹ matrices
    G     = zeros(2, 2)
    Rv    = Array{Array{vartype,2}, 1}(undef, n)
    Vv    = Array{Array{vartype,2}, 1}(undef, n)
    iVv   = Array{Array{vartype,2}, 1}(undef, n)
    matvecz!(Rv, Zv)
    matvecz!(Vv, Zv)
    matvecz!(iVv, Zv)
    #Optimization
    pO      = nothing
    td      = TwiceDifferentiable(x -> -2*remlb(yv, Zv, p, Xv, varlink(x, vlm), β; memopt = memopt), θvec0; autodiff = :forward)
    opttry  = true
    optnum  = 0
    rng     = MersenneTwister(hash(θvec0))
    while opttry
        try
            O       = optimize(td, θvec0, method=Newton(),  g_tol=g_tol, x_tol=x_tol, f_tol=f_tol, allow_f_increases = true, store_trace = store_trace, extended_trace = extended_trace, show_trace = show_trace)
            opttry  = false
        catch

            θvec0 = rvarlink(abs.(varlink(θvec0, vlm) .+ (rand(rng)-0.5)/20 .* varlink(θvec0, vlm) .+ eps()), vlm)[1:4]
            push!(θvec0, rand(rng))
            #θvec0[5] = θvec0[5] - rhoadjstep

        end
        optnum += 1
        if optnum > maxopttry
            opttry = false
            throw(ErrorException("Optimization faild! Iteration $(optnum), θvec = $(varlink(θvec0, vlm))"))
        end

    end
    θ          = Optim.minimizer(O)

    #Post optimization
    if postopt
        pO     = O
        od     = OnceDifferentiable(x -> -2*reml(yv, Zv, p, Xv, varlink(x, vlm), β; memopt = memopt), θ; autodiff = :forward)
        method = BFGS(linesearch = LineSearches.HagerZhang(), alphaguess = LineSearches.InitialStatic())
        O      = optimize(od, [-Inf, -Inf, -Inf, -Inf, -Inf], [Inf, Inf, Inf, Inf, Inf], θ,  Fminbox(method), Optim.Options(g_tol=g_tol, x_tol=x_tol, f_tol=f_tol))
        θ      = Optim.minimizer(O)
    end
    θ = varlink(θ, vlm)
    #Get reml
    remlv       = -reml2b!(yv, Zv, p, n, N, Xv, G, Rv, Vv, iVv, θ, β, memalloc)
    #remlv       = -reml2b!(yv, Zv, p, n, N, Xv, G, Rv, Vv, iVv, varlink(θ, vlm), β, memalloc)
    #Get Hessian matrix (H) with ForwardDiff
    H           = ForwardDiff.hessian(x -> -2*reml(yv, Zv, p, Xv, x, β), θ)
    #H           = ForwardDiff.hessian(x -> -2*reml(yv, Zv, p, Xv, varlink(x, vlm), β), θ)
    #H           = O.trace[end].metadata["h(x)"]
    #θ           = varlink(θ, vlm)
    A = nothing
    #If rho is near to 1.0 it leads to singular hessian matrix, and rho should be removed from variance-covariance matrix
    #It can be done another way: using varlink everywhere, but it leads to problems of calling varlink after RBE object creation with other methods
    if abs(θ[5]) > 1 - singlim
        θ[5]    = 1.0
        H[:,5] .= 0
        H[5,:] .= 0
        A       = 2 * pinv(H)
    else
        A       = 2 * inv(H)
    end

    dH          = det(H)
    #Secondary parameters calculation
    #=
    if abs(dH) > singlim
        A       = 2 * inv(H)
    else
        A       = 2 * pinv(H)
    end
    =#
    C           = cmat(Xv, Zv, iVv, θ)
    se          = Array{vartype, 1}(undef, p)
    F           = Array{vartype, 1}(undef, p)
    df          = Array{vartype, 1}(undef, p)
    t           = Array{vartype, 1}(undef, p)
    pval        = Array{vartype, 1}(undef, p)
    for i = 1:p
        L       = zeros(1, p)
        L[i]    = 1
        Lt      = L'
        lcl     = L*C*Lt                                                        #lcl     = L*C*L'
        lclr    = rank(lcl)
        se[i]   = sqrt((lcl)[1])
        F[i]    = β'*L'*inv(lcl)*L*β/lclr                                       #F[i]    = (L*β)'*inv(L*C*L')*(L*β)/lclr
        g       = ForwardDiff.gradient(x -> lclgf(L, Lt, Xv, Zv, x; memopt = memopt), θ)
        df[i]   = max(1, 2*((lcl)[1])^2/(g'*(A)*g))
        t[i]    = ((L*β)/se[i])[1]
        pval[i] = ccdf(TDist(df[i]), abs(t[i]))*2
        #LinearAlgebra.eigen(L*C*L')
    end
    fixed       = EffectTable(coefnames(MF), β, se, F, df, t, pval)
    fac         = [sequence, period, formulation]
    F           = Array{vartype, 1}(undef, length(fac))
    df          = Array{vartype, 1}(undef, length(fac))
    ndf         = Array{vartype, 1}(undef, length(fac))
    pval        = Array{vartype, 1}(undef, length(fac))
    for i = 1:length(fac)
        L       = lmatrix(MF, fac[i])
        lcl     = L*C*L'
        lclr    = rank(lcl)
        F[i]    = β'*L'*inv(lcl)*L*β/lclr
        if lclr ≥ 2
            vm  = Array{vartype, 1}(undef, lclr)
            for i = 1:lclr
                g        = ForwardDiff.gradient(x -> lclgf(L[i:i,:], L[i:i,:]', Xv, Zv, x; memopt = memopt), θ)
                dfi      = 2*((L[i:i,:]*C*L[i:i,:]')[1])^2/(g'*A*g)
                vm[i]    = dfi/(dfi-2)
            end
            dfi = 2*sum(vm)/(sum(vm)-lclr)
        else
            g   = ForwardDiff.gradient(x -> lclgf(L, L', Xv, Zv, x; memopt = memopt), θ)
            dfi = 2*((lcl)[1])^2/(g'*(A)*g)
        end
        df[i]   = max(1, dfi)
        ndf[i]  = numdf[fac[i]]
        pval[i] = ccdf(FDist(ndf[i], df[i]), F[i])
    end
    typeiii     = ContrastTable(fac, F, ndf, df, pval)
    design      = Design(N, n,
    termmodelleveln(MF, sequence),
    termmodelleveln(MF, period),
    termmodelleveln(MF, formulation),
    sbf,
    p, zxr)
    return RBE(MF, RMF, design, fac, varlink(θvec0, vlm), vlm, Tuple(θ), remlv, fixed, typeiii, Rv, Vv, G, C, A, H, X, Z, Xv, Zv, yv, dH, pO, O)

end #END OF rbe()
"""
This function can convert non-categorical to categorical and sort dataframe.


It can takes more time, but can help to avoid some errors like: "ERROR: type ContinuousTerm has no field contrasts".
"""
function rbe!(df; dvar::Symbol,
    subject::Symbol,
    formulation::Symbol,
    period::Symbol,
    sequence::Symbol,
    g_tol::Real = 1e-12, x_tol::Real = 0.0, f_tol::Real = 0.0, iterations::Int = 100,
    store_trace = false, extended_trace = false, show_trace = false,
    memopt = true,
    init = [],
    postopt = false, vlm = 1.0, maxopttry = 50, rhoadjstep = 0.15,
    rholink = :psigmoid,
    singlim = 1e-6)


    if any(x -> x ∉ names(df), [subject, formulation, period, sequence]) throw(ArgumentError("Names not found in DataFrame!")) end
    if !(eltype(df[!,dvar]) <: AbstractFloat)
        @warn "Responce variable ∉ AbstractFloat!"
        df[!,dvar] = float.(df[!,dvar])
    end
    if !(typeof(df[!,subject]) <: CategoricalArray)
        categorical!(df, subject);
    end
    if !(typeof(df[!,formulation]) <: CategoricalArray)
        categorical!(df, formulation);
    end
    if !(typeof(df[!,period]) <: CategoricalArray)
        categorical!(df, period);
    end
    if !(typeof(df[!,sequence]) <: CategoricalArray)
        categorical!(df, sequence);
    end
    sort!(df, [subject, formulation, period])
    return rbe(df, dvar=dvar, subject=subject, formulation=formulation, period=period, sequence=sequence,
    g_tol=g_tol, x_tol=x_tol, f_tol=f_tol, iterations=iterations,
    store_trace=store_trace, extended_trace=extended_trace, show_trace=show_trace,
    memopt=memopt, init=init, postopt=postopt, vlm = vlm, maxopttry = maxopttry, rhoadjstep = rhoadjstep,
    rholink = rholink, singlim = singlim)

end
#-------------------------------------------------------------------------------
#returm -2REML
"""
    reml2(rbe::RBE, θ::Array{T, 1}) where T <: AbstractFloat

Returm -2logREML for rbe model with θ variance vector.
"""
function reml2(rbe::RBE, θ::Array{T, 1}) where T <: AbstractFloat
    return -2*reml(rbe.yv, rbe.Zv, rank(ModelMatrix(rbe.model).m), rbe.Xv, θ, coef(rbe))
end
"""
    reml2(rbe::RBE)

Returm -2logREML for rbe model

```math
logREML(\\theta,\\beta) = -\\frac{N-p}{2} - \\frac{1}{2}\\sum_{i=1}^nlog|V_{i}|-\n

-\\frac{1}{2}log|\\sum_{i=1}^nX_i'V_i^{-1}X_i|-\\frac{1}{2}\\sum_{i=1}^n(y_i - X_{i}\\beta)'V_i^{-1}(y_i - X_{i}\\beta)
```

"""
function reml2(rbe::RBE)
    return rbe.reml
end
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
"""
    StatsBase.nobs(rbe::RBE)

Return number of statistically independent observations (subjects).
"""
function StatsBase.nobs(rbe::RBE)
    return rbe.design.subj
end
"""
    StatsBase.stderror(rbe::RBE)

Return the standard errors for the coefficients of the model.

```math
se = \\sqrt{LCL'}
```

where

```math
C = (\\sum_{i=1}^{n} X_i'V_i^{-1}X_i)^{-1}
```
"""
function StatsBase.stderror(rbe::RBE)
    return collect(rbe.fixed.se)
end
"""
    coef(rbe::RBE)

Return model coefficients.

```math
\\beta = {(\\sum_{i=1}^n X_{i}'V_i^{-1}X_{i})}^{-1}(\\sum_{i=1}^n X_{i}'V_i^{-1}y_{i})
```
"""
function StatsBase.coef(rbe::RBE)
    return collect(rbe.fixed.est)
end
"""
    dof(rbe::RBE)

Return the number of degrees of freedom for the coefficients of the model.

"""
function StatsBase.dof(rbe::RBE)
    return collect(rbe.fixed.df)
end
#Confidence interval
"""
```julia
confint(obj::RBE; level::Real=0.95, expci::Bool = false, inv::Bool = false, df = :sat)
```

Compute confidence intervals for coefficients, with confidence level ```level``` (by default 95%).

# Arguments

* ```expci = true```: return exponented CI.

* ```inv = true```: return ```-estimate ± t(alpha, df)*SE```

* ```df = :sat```: use Satterthwaite DF approximation.

* ```df = :df3``` or ```df = :cont```: DF (contain) = N - rank(ZX).

```math
CI = estimate ± t(alpha, df)*SE
```

"""
function StatsBase.confint(obj::RBE; level::Real=0.95, expci::Bool = false, inv::Bool = false, df = :sat)
    confint(obj, 1 - level; expci = expci, inv = inv, df = df)
end
function StatsBase.confint(obj::RBE, alpha::Real; expci::Bool = false, inv::Bool = false, df = :sat)
    ifv = 1
    if inv ifv = -1 end
    if isa(df, Array)
        if length(obj.fixed.df) != length(df)
            df = obj.fixed.df
        else
            @warn "length(df) not equal parameters count, default df used!"
            df = obj.fixed.df
        end
    elseif isa(df, Symbol)
        if df == :df2
            df  = zeros(length(obj.fixed.df))
            df .= obj.design.df2
        elseif df == :df3 || df == :cont
            df  = zeros(length(obj.fixed.df))
            df .= obj.design.df3
        elseif df == :contw
            df  = zeros(length(obj.fixed.df))
            df .= sum(obj.design.sbf) - length(obj.design.sbf)*obj.design.sqn
        elseif df == :sat
            df = obj.fixed.df
        else
            @warn "df unknown, default df used!"
            df = obj.fixed.df
        end
    end
    a = Array{Tuple{Float64, Float64},1}(undef, length(obj.fixed.est))
    for i = 1:length(obj.fixed.est)
        a[i] = calcci(obj.fixed.est[i]*ifv, obj.fixed.se[i], df[i], alpha, expci)
    end
    return Tuple(a)
end
function calcci(x::T, se::T, df::T, alpha::T, expci::Bool) where T <: AbstractFloat
    q = quantile(TDist(df), 1.0-alpha/2)
    if !expci
        return x-q*se, x+q*se
    else
        return exp(x-q*se), exp(x+q*se)
    end
end
#-------------------------------------------------------------------------------
"""
    theta(rbe::RBE)

Return theta (θ) vector (vector of variation parameters from optimization procedure).
"""
function theta(rbe::RBE)
    return collect(rbe.θ)
end
"""
    coefnum(rbe::RBE)

Return number of coefficients (length β).
"""
function coefnum(rbe::RBE)
    return length(rbe.fixed.se)
end
"""
    design(rbe::RBE)::Design

Return design information object.
"""
function design(rbe::RBE)::Design
    return rbe.design
end
"""
    fixed(rbe::RBE)

Return fixed effect table (β).
"""
function fixed(rbe::RBE)
    return rbe.fixed
end
"""
    typeiii(rbe::RBE)

Return TYPE III table.

(see contrast)
"""
function typeiii(rbe::RBE)
    return rbe.typeiii
end
"""
    optstat(rbe::RBE)

Return optimization status.
* true - converged
* false - not converged
"""
function optstat(rbe::RBE)
    return Optim.converged(rbe.optim)
end
#-------------------------------------------------------------------------------
function Base.show(io::IO, rbe::RBE)
    rcoef = coefnames(rbe.rmodel);
    θ     = theta(rbe)
    println(io, "Bioequivalence Linear Mixed Effect Model (status: $(Optim.converged(rbe.optim) ? "converged" : printstyled(io, "not converged"; color = :red)))")
    if !isposdef(Symmetric(rbe.H))
        printstyled(io, "Hessian not positive defined!"; color = :yellow)
        println(io, "")
    end
    if θ[end] >= 1.0 - eps()
        printstyled(io, "Rho ~ 1.0!"; color = :yellow)
        println(io, "")
    end
    println(io, "")
    println(io, "-2REML: $(round(rbe.reml, sigdigits=6))    REML: $(round(-rbe.reml/2, sigdigits=6))")
    println(io, "")
    println(io, "Fixed effect:")
    println(io, rbe.fixed)
    println(io, "Intra-individual variance:")

    printmatrix(io,[rcoef[1] round(θ[1], sigdigits=6) "CVᵂ:" round(geocv(θ[1])*100, digits=2) "%";
                    rcoef[2] round(θ[2], sigdigits=6) "CVᵂ:" round(geocv(θ[2])*100, digits=2) "%"])
    println(io, "")

    println(io, "Inter-individual variance:")

    printmatrix(io,[rcoef[1] round(θ[3], sigdigits=6) "";
                    rcoef[2] round(θ[4], sigdigits=6) "";
                    "ρ:"     round(θ[5], sigdigits=6) "Cov: $(round(sqrt(θ[4]*θ[3])*θ[5], sigdigits=6))"])
    println(io, "")

    println(io, "Confidence intervals(90%):")
    ci = confint(rbe, 0.1, expci = true, inv = true)
    println(io, rcoef[1], " / ", rcoef[2])
    println(io, "Ratio: $(round(exp(-coef(rbe)[end])*100, digits=2)), CI: ", round(ci[end][1]*100, digits=2), " - ", round(ci[end][2]*100, digits=2), " (%)")
    ci = confint(rbe, 0.1, expci = true, inv = false)
    println(io, rcoef[2], " / ", rcoef[1])
    print(io,"Ratio: $(round(exp(coef(rbe)[end])*100, digits=2)), CI: ", round(ci[end][1]*100, digits=2), " - ", round(ci[end][2]*100, digits=2), " (%)")
end #─┼┴┬│
