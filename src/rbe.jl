#
# GENARAL REPLICATE BIOEQUIVALENCE STRUCTURE
#
struct RBE
    model::ModelFrame               #Model frame
    rmodel::ModelFrame              #Random effect model
    design::Design
    factors::Array{Symbol, 1}       #Factor list
    #β::Array{Float64, 1}            #β coefficients (fixed effect)
    θ0::Array{Float64, 1}           #Initial variance paramethers
    θ::Array{Float64, 1}            #Final variance paramethers
    reml::Float64                   #-2REML
    fixed::EffectTable
    typeiii::ContrastTable
    #se::Array{Float64, 1}           #SE for each β level
    #f::Array{Float64, 1}            #F for each β level
    #df::Array{Float64, 1}           #DF (degree of freedom) for each β level (Satterthwaite)
    #df2::Float64                    #DF N / pn - sn
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
    preoptim::Optim.MultivariateOptimizationResults        #Pre-optimization result object
    optim::Optim.MultivariateOptimizationResults           #Optimization result object
end

"""
    Mixed model fitting function for replicate bioequivalence
"""
function rbe(df; dvar::Symbol,
    subject::Symbol,
    formulation::Symbol,
    period::Symbol,
    sequence::Symbol,
    g_tol::Float64 = 1e-8, x_tol::Float64 = 0.0, f_tol::Float64 = 0.0, iterations::Int = 100,
    store_trace = false, extended_trace = false, show_trace = false,
    memopt = true)

    if any(x -> x ∉ names(df), [subject, formulation, period, sequence]) throw(ArgumentError("Names not found in DataFrame!")) end
    if !(eltype(df[!,dvar]) <: Real) println("Responce variable not a float!") end

    to = TimerOutput()
    categorical!(df, subject);
    categorical!(df, formulation);
    categorical!(df, period);
    categorical!(df, sequence);
    @timeit to "sort" sort!(df, [subject, formulation, period])
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
    @timeit to "sortsubj" Xv, Zv, yv = sortsubjects(df, subject, X, Z, y)
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
    @timeit to "check" checkdata(X, Z, Xv, Zv, y)
    #Calculate initial fixed parameters
    qro   = qr(X)
    @timeit to "lm"  β     = inv(qro.R)*qro.Q'*y
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
    method = LBFGS()
    #method = ConjugateGradient()
    #method = NelderMead()
    limeps=eps()
    @timeit to "o1" pO = optimize(od, [limeps, limeps, limeps, limeps, limeps], [Inf, Inf, Inf, Inf, 1.0], θvec0,  Fminbox(method), Optim.Options(g_tol = 1e-1))
    #pO = optimize(remlf,  [limeps, limeps, limeps, limeps, limeps], [Inf, Inf, Inf, Inf, 1.0], θvec0,  Fminbox(method), Optim.Options(g_tol = 1e-3))
    θ  = Optim.minimizer(pO)
    #Final optimization
    #Provide gradient function for Optim
    #Not used yet
    #g!(storage, θx) = copyto!(storage, ForwardDiff.gradient(x -> -2*reml(yv, Zv, p, Xv, x, β), θx))
    #REML function for optimization
    td = TwiceDifferentiable(x -> -2*remlb(yv, Zv, p, Xv, x, β; memopt = memopt), θvec0; autodiff = :forward)
    #remlfb(x) = -reml2b!(yv, Zv, p, n, N, Xv, G, Rv, Vv, iVv, x, β, memc, memc2, memc3, memc4)
    #@timeit to "o2" O  = optimize(remlfb, θ, method=Newton(),  g_tol=g_tol, x_tol=x_tol, f_tol=f_tol, allow_f_increases = true, store_trace = store_trace, extended_trace = extended_trace, show_trace = show_trace)
    @timeit to "o2" O  = optimize(td, θ, method=Newton(),  g_tol=g_tol, x_tol=x_tol, f_tol=f_tol, allow_f_increases = true, store_trace = store_trace, extended_trace = extended_trace, show_trace = show_trace)
    θ  = Optim.minimizer(O)
    #Get reml
    #remlv = -2*reml(yv, Zv, p, Xv, θ, β)
    @timeit to "remlv"  remlv = -reml2b!(yv, Zv, p, n, N, Xv, G, Rv, Vv, iVv, θ, β, memalloc)
    #θ[5] can not be more than 1.0
    if θ[5] > 1 θ[5] = 1 end
    #Get Hessian matrix (H) with ForwardDiff
    @timeit to "H" H         = ForwardDiff.hessian(x -> -2*reml(yv, Zv, p, Xv, x, β), θ)
    dH          = det(H)
    H[5,:]     .= 0
    H[:,5]     .= 0
    #Secondary parameters calculation
    A           = 2*pinv(H)
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
        df[i]   = 2*((lcl)[1])^2/(g'*(A)*g)
        t[i]    = ((L*β)/se[i])[1]
        pval[i] = ccdf(TDist(df[i]), abs(t[i]))*2
        #LinearAlgebra.eigen(L*C*L')
    end
    fixed       = EffectTable(coefnames(MF), β, se, F, df, t, pval)
    fac         = [sequence, period, formulation]
    F           = Array{Float64, 1}(undef, length(fac))
    df          = Array{Float64, 1}(undef, length(fac))
    pval        = Array{Float64, 1}(undef, length(fac))
    for i = 1:length(fac)
        L       = lmatrix(MF, fac[i])
        Lt      = L'
        lcl     = L*C*Lt
        lclr    = rank(lcl)
        #M       = L'*inv(L*L')*L
        #t1      = tr(M*C)
        #v1      = t1^2/tr(M*C*M*C)
        #F[i]    = β'*M*β/t1
        #df[i]   = 2*(t1)^2/(g'*(A)*g)
        F[i]    = β'*L'*inv(lcl)*L*β/lclr
        g       = ForwardDiff.gradient(x -> lclgf(L, Lt, Xv, Zv, x; memopt = memopt), θ)
        df[i]   = 2*((lcl)[1])^2/(g'*(A)*g)
        pval[i] = ccdf(FDist(numdf[fac[i]], df[i]), F[i])
    end
    typeiii     = ContrastTable(fac, F, df, pval)
    design      = Design(N, n,
    termmodelleveln(MF, sequence),
    termmodelleveln(MF, period),
    termmodelleveln(MF, formulation),
    sbf,
    p, zxr, n - termmodelleveln(MF, sequence), N - zxr, N - zxr + p)
    #println(to)
    return RBE(MF, RMF, design, fac, θvec0, θ, remlv, fixed, typeiii, Rv, Vv, G, C, A, H, X, Z, Xv, Zv, yv, dH, pO, O)
end #END OF rbe()
#-------------------------------------------------------------------------------
#returm -2REML
"""
    Returm -2REML for rbe model object and θ variance parameters vector
"""
function reml2(rbe::RBE, θ::Array{Float64, 1})
    return -2*reml(rbe.yv, rbe.Zv, rank(ModelMatrix(rbe.model).m), rbe.Xv, θ, rbe.fixed.est)
end
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
"""
    Return model coefficients
"""
function StatsBase.coef(rbe::RBE)
    return copy(rme.fixed.est)
end
#Confidence interval
function StatsBase.confint(obj::RBE, alpha::Float64; expci::Bool = false, inv::Bool = false, df = :sat)
    ifv = 1
    if inv ifv = -1 end
    if isa(df, Array{Float64, 1})
        if length(obj.fixed.df) != length(df)
            df = obj.fixed.df
        end
    elseif isa(df, Symbol)
        if df == :df2
            df  = zeros(length(obj.fixed.df))
            df .= obj.design.df2
        else
            df = obj.fixed.df
        end
    end
    a = Array{Tuple{Float64, Float64},1}(undef, length(obj.fixed.est)-1)
    for i = 2:length(obj.fixed.est)
        a[i-1] = calcci(obj.fixed.est[i]*ifv, obj.fixed.se[i], df[i], alpha, expci)
    end
    return Tuple(a)
end

#-------------------------------------------------------------------------------
function Base.show(io::IO, rbe::RBE)
    rcoef = coefnames(rbe.rmodel);
    println(io, "Bioequivalence Linear Mixed Effect Model")
    println(io, "")
    println(io, "-2REML: $(round(rbe.reml, sigdigits=6))    REML: $(round(-rbe.reml/2, sigdigits=6))")
    println(io, "")
    println(io, "Fixed effect:")
    println(io, rbe.fixed)
    println(io, "Intra-individual variation:")
    println(io, rcoef[1], "  ", round(rbe.θ[1], sigdigits=6), "   CVᵂ: ", round(geocv(rbe.θ[1]), sigdigits=6))
    println(io, rcoef[2], "  ", round(rbe.θ[2], sigdigits=6), "   CVᵂ: ", round(geocv(rbe.θ[2]), sigdigits=6))
    println(io, "")
    println(io, "Inter-individual variation:")
    println(io, rcoef[1], "  ", round(rbe.θ[3], sigdigits=6))
    println(io, rcoef[2], "  ", round(rbe.θ[4], sigdigits=6))
    println(io,   "Cov:", "  ", round(sqrt(rbe.θ[4]*rbe.θ[3])*rbe.θ[5], sigdigits=6))
    println(io, "")
    println(io, "Confidence intervals(90%):")
    ci = confint(rbe, 0.1, expci = true, inv = false)
    println(io, rcoef[1], " / ", rcoef[2])
    println(io, round(ci[end][1]*100, digits=4), " - ", round(ci[end][2]*100, digits=4), " (%)")
    ci = confint(rbe, 0.1, expci = true, inv = true)
    println(io, rcoef[2], " / ", rcoef[1])
    print(io, round(ci[end][1]*100, digits=4), " - ", round(ci[end][2]*100, digits=4), " (%)")
end #─┼┴┬│
