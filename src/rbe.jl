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

    rbe(df; dvar::Symbol,
        subject::Symbol,
        formulation::Symbol,
        period::Symbol,
        sequence::Symbol,
        g_tol::Float64 = 1e-8, x_tol::Float64 = 0.0, f_tol::Float64 = 0.0, iterations::Int = 100,
        store_trace = false, extended_trace = false, show_trace = false,
        memopt = true)

Mixed model fitting function for replicate bioequivalence without data preparation (apply categorical! for each factor and sort! to dataframe).

Mixed model in matrix form:

```math
y = X\\beta + Zu + \\epsilon
```


with covariance matrix for each subject:

```math
V_{i} = Z_{i}GZ_i'+R_{i}
```

"""
function rbe(df; dvar::Symbol,
    subject::Symbol,
    formulation::Symbol,
    period::Symbol,
    sequence::Symbol,
    g_tol::Float64 = 1e-8, x_tol::Float64 = 0.0, f_tol::Float64 = 0.0, iterations::Int = 100,
    store_trace = false, extended_trace = false, show_trace = false,
    memopt = true,
    init = [],
    twostep = true,
    postopt = false)

    if any(x -> x ∉ names(df), [subject, formulation, period, sequence]) throw(ArgumentError("Names not found in DataFrame!")) end
    if !(eltype(df[!,dvar]) <: Real)
        @warn "Responce variable ∉ Real!"
    end

    #should not change initial DS
    #categorical!(df, subject);
    #categorical!(df, formulation);
    #categorical!(df, period);
    #categorical!(df, sequence);
    #sort!(df, [subject, formulation, period])
    Xf  = @eval(@formula($dvar ~ $sequence + $period + $formulation))
    Zf  = @eval(@formula($dvar ~ 0 + $formulation))
    MF  = ModelFrame(Xf, df)
    RMF = ModelFrame(Zf, df, contrasts = Dict(formulation => StatsModels.FullDummyCoding()))
    MM  = ModelMatrix(MF)
    X   = MM.m
    Z   = ModelMatrix(RMF).m
    p   = rank(X)
    zxr = rank(ModelMatrix(ModelFrame(@eval(@formula($dvar ~ $sequence + $period + $subject*$formulation)), df)).m)
    y   = df[:, dvar]                                                            #Dependent variable
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
    iv = initvar(df, dvar, formulation, subject)
    if iv[1] < iv[3] || iv[2] < iv[3] iv[1] = iv[2] = 2*iv[3] end
    θvec0 = [iv[3], iv[3], iv[1]-iv[3], iv[2]-iv[3], 0.501]
    #Prelocatiom for G, R, V, V⁻¹ matrices
    G     = zeros(2, 2)
    Rv    = Array{Array{Float64,2}, 1}(undef, n)
    Vv    = Array{Array{Float64,2}, 1}(undef, n)
    iVv   = Array{Array{Float64,2}, 1}(undef, n)
    matvecz!(Rv, Zv)
    matvecz!(Vv, Zv)
    matvecz!(iVv, Zv)
    #First step optimization (pre-optimization)
    od = OnceDifferentiable(x -> -2*reml(yv, Zv, p, Xv, x, β; memopt = memopt), θvec0; autodiff = :forward)
    method = BFGS(linesearch = LineSearches.HagerZhang(), alphaguess = LineSearches.InitialStatic())
    #method = optm
    #method = ConjugateGradient()
    #method = NelderMead()
    limeps  = eps()
    pO      = nothing
    if twostep
        pO = optimize(od, [limeps, limeps, limeps, limeps, limeps], [Inf, Inf, Inf, Inf, 1.0], θvec0,  Fminbox(method), Optim.Options(g_tol = 1e-1))
    #pO = optimize(remlf,  [limeps, limeps, limeps, limeps, limeps], [Inf, Inf, Inf, Inf, 1.0], θvec0,  Fminbox(method), Optim.Options(g_tol = 1e-3))
        θ  = Optim.minimizer(pO)
    else
        #@warn "First optimization step failed. Start step two with initial θ..."
        θ  = θvec0
    end
    #Final optimization
    #Provide gradient function for Optim
    #Not used yet
    #g!(storage, θx) = copyto!(storage, ForwardDiff.gradient(x -> -2*reml(yv, Zv, p, Xv, x, β), θx))
    #REML function for optimization
    td = TwiceDifferentiable(x -> -2*remlb(yv, Zv, p, Xv, x, β; memopt = memopt), θvec0; autodiff = :forward)
    #remlfb(x) = -reml2b!(yv, Zv, p, n, N, Xv, G, Rv, Vv, iVv, x, β, memc, memc2, memc3, memc4)
    #@timeit to "o2" O  = optimize(remlfb, θ, method=Newton(),  g_tol=g_tol, x_tol=x_tol, f_tol=f_tol, allow_f_increases = true, store_trace = store_trace, extended_trace = extended_trace, show_trace = show_trace)
    O  = optimize(td, θ, method=Newton(),  g_tol=g_tol, x_tol=x_tol, f_tol=f_tol, allow_f_increases = true, store_trace = store_trace, extended_trace = extended_trace, show_trace = show_trace)
    θ  = copy(Optim.minimizer(O))
    #Get reml
    #remlv = -2*reml(yv, Zv, p, Xv, θ, β)
    remlv = -reml2b!(yv, Zv, p, n, N, Xv, G, Rv, Vv, iVv, θ, β, memalloc)

    #θ[5] can not be more than 1.0
    if θ[5] >= 1.0
        θ[5] = 1.0 - eps()
        if !twostep && !postopt
            @warn "ρ is more than 1.0, and no twostep or postopt used. Results may be incorrect, use twostep = true or postopt = true"
        end
        if postopt
            O  = optimize(od, [limeps, limeps, limeps, limeps, limeps], [Inf, Inf, Inf, Inf, 1.0], θ,  Fminbox(method), Optim.Options(g_tol=g_tol, x_tol=x_tol, f_tol=f_tol))
            θ  = copy(Optim.minimizer(O))
            remlv = -reml2b!(yv, Zv, p, n, N, Xv, G, Rv, Vv, iVv, θ, β, memalloc)
        end
    end
    #Get Hessian matrix (H) with ForwardDiff
    #H           = Optim.trace(O)[end].metadata["h(x)"]
    H           = ForwardDiff.hessian(x -> -2*reml(yv, Zv, p, Xv, x, β), θ)
    dH          = det(H)


    H[5,:]     .= eps()
    H[:,5]     .= eps()

    #Secondary parameters calculation
    A           = 2*pinv(H)
    #A[5,:]     .= 0
    #A[:,5]     .= 0
    C           = cmat(Xv, Zv, iVv, θ)
    se          = Array{Float64, 1}(undef, p)
    F           = Array{Float64, 1}(undef, p)
    df          = Array{Float64, 1}(undef, p)
    t           = Array{Float64, 1}(undef, p)
    pval        = Array{Float64, 1}(undef, p)
    for i = 1:p
        L       = zeros(1, p)
        L[i]    = 1
        Lt      = L'
        lcl     = L*C*Lt                                                        #lcl     = L*C*L'
        lclr    = rank(lcl)
        se[i]   = sqrt((lcl)[1])
        #Lβ      = L*β
        F[i]    = β'*L'*inv(lcl)*L*β/lclr                                       #F[i]    = (L*β)'*inv(L*C*L')*(L*β)/lclr
        #lclg    = lclgf(L, Lt, Xv, Zv, x; memopt = memopt)
        g       = ForwardDiff.gradient(x -> lclgf(L, Lt, Xv, Zv, x; memopt = memopt), θ)
        df[i]   = max(1, 2*((lcl)[1])^2/(g'*(A)*g))
        t[i]    = ((L*β)/se[i])[1]
        pval[i] = ccdf(TDist(df[i]), abs(t[i]))*2
        #LinearAlgebra.eigen(L*C*L')
    end
    fixed       = EffectTable(coefnames(MF), β, se, F, df, t, pval)
    fac         = [sequence, period, formulation]
    F           = Array{Float64, 1}(undef, length(fac))
    df          = Array{Float64, 1}(undef, length(fac))
    ndf         = Array{Float64, 1}(undef, length(fac))
    pval        = Array{Float64, 1}(undef, length(fac))
    for i = 1:length(fac)
        L       = lmatrix(MF, fac[i])
        lcl     = L*C*L'
        lclr    = rank(lcl)
        F[i]    = β'*L'*inv(lcl)*L*β/lclr
        if lclr ≥ 2
            vm  = Array{Float64, 1}(undef, lclr)
            for i = 1:lclr
                g        = ForwardDiff.gradient(x -> lclgf(L[i:i,:], L[i:i,:]', Xv, Zv, x; memopt = memopt), θ)
                dfi      = 2*((L[i:i,:]*C*L[i:i,:]')[1])^2/(g'*(A)*g)
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
    return RBE(MF, RMF, design, fac, θvec0, Tuple(θ), remlv, fixed, typeiii, Rv, Vv, G, C, A, H, X, Z, Xv, Zv, yv, dH, pO, O)
end #END OF rbe()
"""
This function apply following code for each factor before executing:

```julia
    categorical!(df, subject);
    categorical!(df, formulation);
    categorical!(df, period);
    categorical!(df, sequence);
    sort!(df, [subject, formulation, period])
```

It can takes more time, but can help to avoid some errors like: "ERROR: type ContinuousTerm has no field contrasts".
"""
function rbe!(df; dvar::Symbol,
    subject::Symbol,
    formulation::Symbol,
    period::Symbol,
    sequence::Symbol,
    g_tol::Float64 = 1e-8, x_tol::Float64 = 0.0, f_tol::Float64 = 0.0, iterations::Int = 100,
    store_trace = false, extended_trace = false, show_trace = false,
    memopt = true,
    init = [],
    twostep = true,
    postopt = false)

    if any(x -> x ∉ names(df), [subject, formulation, period, sequence]) throw(ArgumentError("Names not found in DataFrame!")) end
    if !(eltype(df[!,dvar]) <: Real)
        @warn "Responce variable ∉ Real!"
        df[!,dvar] = float.(df[!,dvar])
    end

    categorical!(df, subject);
    categorical!(df, formulation);
    categorical!(df, period);
    categorical!(df, sequence);
    sort!(df, [subject, formulation, period])

    return rbe(df, dvar=dvar, subject=subject, formulation=formulation, period=period, sequence=sequence,
    g_tol=g_tol, x_tol=x_tol, f_tol=f_tol, iterations=iterations,
    store_trace=store_trace, extended_trace=extended_trace, show_trace=show_trace,
    memopt=memopt, init=init, twostep=twostep, postopt=postopt)
end
#-------------------------------------------------------------------------------
#returm -2REML
"""
    reml2(rbe::RBE, θ::Array{Float64, 1})

Returm -2REML for rbe model
"""
function reml2(rbe::RBE, θ::Array{Float64, 1})
    return -2*reml(rbe.yv, rbe.Zv, rank(ModelMatrix(rbe.model).m), rbe.Xv, θ, coef(rbe))
end
"""
    reml2(rbe::RBE)

Returm -2logREML for rbe model

```
logREML(\\theta,\\beta) = -\\frac{N-p}{2} - \\frac{1}{2}\\sum_{i=1}^nlog|V_{i}|-

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
    confint(obj::RBE; level::Real=0.95, expci::Bool = false, inv::Bool = false, df = :sat)

Compute confidence intervals for coefficients, with confidence level ```level``` (by default 95%).

# Arguments

```expci = true```: return exponented CI.

```inv = true```: return ```-estimate ± t(alpha, df)*SE```

```df = :sat```: use Satterthwaite DF approximation.

```df = :df3``` or ```df = :cont```: DF (contain) = N - rank(ZX).

```math
CI = estimate ± t(alpha, df)*SE
```

"""
function StatsBase.confint(obj::RBE; level::Real=0.95, expci::Bool = false, inv::Bool = false, df = :sat)
    confint(obj, 1 - level; expci = expci, inv = inv, df = df)
end
function StatsBase.confint(obj::RBE, alpha::Float64; expci::Bool = false, inv::Bool = false, df = :sat)
    ifv = 1
    if inv ifv = -1 end
    if isa(df, Array{Float64, 1})
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
function calcci(x::Float64, se::Float64, df::Float64, alpha::Float64, expci::Bool)::Tuple{Float64, Float64}
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

Return design information object, where:

```julia
    struct Design
        obs::Int          # Number of observations
        subj::Int         # Number of statistica independent subjects
        sqn::Int          # Number of sequences
        pn::Int           # Number of periods
        fn::Int           # Number of formulations
        sbf::Vector{Int}  # Subjects in each formulation level
        rankx::Int        # Rank of fixed effect matrix
        rankxz::Int       # Rank of XZ (fixed+random) effect matrix
        df2::Int          # subj - sqn         (Robust DF)
        df3::Int          # obs  - rankxz      (Contain DF for sequence and period)
        df4::Int          # obs  - rankxz + p
    end
```
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

see contrast
"""
function typeiii(rbe::RBE)
    return rbe.typeiii
end
"""
    optstat(rbe::RBE)

Return optimization status.
"""
function optstat(rbe::RBE)
    return Optim.converged(rbe.optim)
end
#-------------------------------------------------------------------------------
function Base.show(io::IO, rbe::RBE)
    rcoef = coefnames(rbe.rmodel);
    println(io, "Bioequivalence Linear Mixed Effect Model (status: $(Optim.converged(rbe.optim) ? "converged" : "not converged"))")
    println(io, "")
    println(io, "-2REML: $(round(rbe.reml, sigdigits=6))    REML: $(round(-rbe.reml/2, sigdigits=6))")
    println(io, "")
    println(io, "Fixed effect:")
    println(io, rbe.fixed)
    println(io, "Intra-individual variation:")

    printmatrix(io,[rcoef[1] round(rbe.θ[1], sigdigits=6) "CVᵂ:" round(geocv(rbe.θ[1]), sigdigits=6);
                    rcoef[2] round(rbe.θ[2], sigdigits=6) "CVᵂ:" round(geocv(rbe.θ[2]), sigdigits=6)])
    println(io, "")

    println(io, "Inter-individual variation:")

    printmatrix(io,[rcoef[1] round(rbe.θ[3], sigdigits=6) "";
                    rcoef[2] round(rbe.θ[4], sigdigits=6) "";
                    "ρ:"     round(rbe.θ[5], sigdigits=6) "Cov: $(round(sqrt(rbe.θ[4]*rbe.θ[3])*rbe.θ[5], sigdigits=6))"])
    println(io, "")

    println(io, "Confidence intervals(90%):")
    ci = confint(rbe, 0.1, expci = true, inv = false)
    println(io, rcoef[1], " / ", rcoef[2])
    println(io, round(ci[end][1]*100, digits=4), " - ", round(ci[end][2]*100, digits=4), " (%)")
    ci = confint(rbe, 0.1, expci = true, inv = true)
    println(io, rcoef[2], " / ", rcoef[1])
    print(io, round(ci[end][1]*100, digits=4), " - ", round(ci[end][2]*100, digits=4), " (%)")
end #─┼┴┬│
