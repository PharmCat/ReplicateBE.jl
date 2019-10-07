*Alpha version!* This program comes with absolutely no warranty. No liability is accepted for any loss and risk to public health resulting from use of this software.

<p align="center">
  <img src="https://github.com/PharmCat/ReplicateBE.jl/blob/master/docs/ReplicateBE-LogoNoSpace.png">
</p>

Mixed model solution for replicate designed bioequivalence study. This can be used to obtained results with methods C (random effects with interaction), given by the EMA in [Annex I](https://www.ema.europa.eu/en/documents/other/31-annex-i-statistical-analysis-methods-compatible-ema-bioequivalence-guideline_en.pdf "EMA/582648/2016, 21 September 2016"). Statistical model formed with accordance [FDA Guidance for Industry: Statistical Approaches to Establishing Bioequivalence](https://www.fda.gov/media/70958/download), APPENDIX F.

[![Build Status](https://api.travis-ci.com/PharmCat/ReplicateBE.jl.svg?branch=master)](https://travis-ci.com/PharmCat/ReplicateBE.jl)
[![codecov](https://codecov.io/gh/PharmCat/ReplicateBE.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/PharmCat/ReplicateBE.jl)


Install:
```
using Pkg; Pkg.add("ReplicateBE")
```
or:
```
using Pkg; Pkg.clone("https://github.com/PharmCat/ReplicateBE.jl.git")
```

Using:
```
using ReplicateBE
be = ReplicateBE.rbe(df, dvar = :var, subject = :subject, formulation = :formulation, period = :period, sequence = :sequence);
ci = confint(be, 0.1)
```
Where:

- dvar::Symbol - dependent variable;
- subject::Symbol - subject;
- formulation::Symbol - formulation/drug;
- period::Symbol - study period;
- sequence::Symbol - sequence;

Get results:
```
be.fixed
#or
fixed(be)

be.typeiii
#or
typeiii(be)
```

Output example:
```
Bioequivalence Linear Mixed Effect Model

-2REML: 164.613    REML: -82.3067

Fixed effect:
─────────────────────────────────────────────────────────────────────────────────────────
Effect           Value        SE         F           DF        t           P|t|
─────────────────────────────────────────────────────────────────────────────────────────
(Intercept)      1.27158      0.191994   43.8645     26.6243   6.62303     4.47193e-7*
sequence: TRTR   -0.114458    0.242518   0.222743    30.4503   -0.471956   0.640323
period: 2        0.035864     0.139526   0.0660707   83.1868   0.257042    0.797781
period: 3        0.0168505    0.100284   0.0282329   41.619    0.168027    0.867376      
period: 4        0.0552216    0.139526   0.156642    83.1868   0.39578     0.693281
formulation: T   -0.0165797   0.120162   0.0190378   62.661    -0.137978   0.890701      
─────────────────────────────────────────────────────────────────────────────────────────
Intra-individual variation:
formulation: R  0.187712   CVᵂ: 0.454407
formulation: T  0.0889289   CVᵂ: 0.304964

Inter-individual variation:
formulation: R  0.318313
formulation: T  0.422413
Cov:  0.266255

Confidence intervals(90%):
formulation: R / formulation: T
80.4773 - 120.2059 (%)
formulation: T / formulation: R
83.1906 - 124.2586 (%)
```

# Validation

Validation information: [here](https://github.com/PharmCat/ReplicateBE.jl/blob/master/validation/validation.md)

# Methods

### StatsBase.confint

```
    StatsBase.confint(obj::RBE, alpha::Float64; expci::Bool = false, inv::Bool = false)
```
Return (1-alpha)×100% confidence intervals for β.

* obj::RBE - bioequivalence struct;
* alpha::Float64 - alpha;
* expci::Bool - exp(ci)
* inv::Bool - β = -β

### StatsBase.coef

```
    StatsBase.coef(rbe::RBE)
```
Return model coefficients.

### ReplicateBE.coefse

```
    coefse(rbe::RBE)
```
Return standard error for coefficients.

### ReplicateBE.design

```
    design(rbe::RBE)::Design
```
Return design information.

### ReplicateBE.fixed

```
    fixed(rbe::RBE)
```
Return fixed effect table.

### ReplicateBE.typeiii

```
    typeiii(rbe::RBE)
```
Return type III effect table.

### ReplicateBE.reml2

```
    reml2(obj::RBE, θ::Array{Float64, 1})
```
Return -2REML for vector θ.

```
    reml2(obj::RBE)
```
Return -2REML for model.

### ReplicateBE.contrast

```
    contrast(obj::RBE, L::Matrix{T}) where T <: Real
```
Return F for L matrix. L matrix should be 1×p or l×p.

### ReplicateBE.lsm

```
    lsm(obj::RBE, L::Matrix{T}) where T <: Real
```
Return mean & se for L matrix. L matrix should be 1×p.

### ReplicateBE.emm

```
    emm(obj::RBE, fm, lm)
```
Return marginal means. fm and lm arrays p length.

fm - factor map;
lm - level map.

### ReplicateBE.emm
```
    lmean(obj::RBE)
```
Return L matrix for general mean.

# Random Dataset

```
randrbeds(;n=24, sequence=[1,1],
    design = ["T" "R" "T" "R"; "R" "T" "R" "T"],
    inter=[0.5, 0.4, 0.9], intra=[0.1, 0.2],
    intercept = 0, seqcoef = [0.0, 0.0], periodcoef = [0.0, 0.0, 0.0, 0.0], formcoef = [0.0, 0.0], seed = 0)
```
Generate DataFrame with random multivariate data. Where:

 - n - Subject number;
 - sequence - sequence subject distribution, [1,1] is equal 1:1, [1,2] - 1:2, [2,3,5] - 2:3:5 ets.;
 - design - design matrix (sXp, where s - number of sequences, p - number of period), cells contains formulation label;
 - inter - Inter-subject variation vector for G matrix: [σ₁, σ₂, ρ], where σ₁, σ₂ - formulation inter-subject variance,  ρ - covariance coefficient;
 - intra - Intra-subject variation vector for R matrix:[σ₁, σ₂], where σ₁, σ₂ - formulation intra-subject variance;
 - intercept - model intercept value;
 - seqcoef - model sequence coefficient values (additive): length = s (number of sequences);
 - periodcoef - model period coefficient values (additive): length = p (number of periods);
 - formcoef - model formulation coefficient values (additive): length = 2;

## Structures

### RBE

```
struct RBE
    model::ModelFrame               #Model frame
    rmodel::ModelFrame              #Random effect model
    design::Design
    factors::Array{Symbol, 1}       #Factor list
    θ0::Array{Float64, 1}           #Initial variance paramethers
    θ::Array{Float64, 1}            #Final variance paramethers
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
    preoptim::Optim.MultivariateOptimizationResults        #Pre-optimization result object
    optim::Optim.MultivariateOptimizationResults           #Optimization result object
end
```

### Design

```
struct Design
    obs::Int
    subj::Int
    sqn::Int
    pn::Int
    fn::Int
    sbf::Vector{Int}
    rankx::Int
    rankxz::Int
    df2::Int
    df3::Int
    df4::Int
end
```

### EffectTable

```
struct EffectTable <: RBETable
    name::Vector
    est::Vector
    se::Vector
    f::Vector
    df::Vector
    t::Vector
    p::Vector
end
```

### ContrastTable

```
struct ContrastTable <: RBETable
    name::Vector
    f::Vector
    df::Vector
    p::Vector
end
```

### EstimateTable

```
struct EstimateTable <: RBETable
    name::Vector
    f::Vector
    df::Vector
    p::Vector
end
```

## Acknowledgments

Best acknowledgments to D.Sc. in Physical and Mathematical Sciences Anastasia Shitova <a.shitova@qayar.ru> for support, datasets and testing procedures.

## References

- LINDSTROM & J.; BATES, M. (1988). Newton—Raphson and EM Algorithms for Linear Mixed-Effects Models for Repeated-Measures Data. Journal of the American Statistical Association. 83. 1014. 10.1080/01621459.1988.10478693.
- Gurka, Matthew. (2006). Selecting the Best Linear Mixed Model under REML. The American Statistician. 60. 19-26. 10.1198/000313006X90396.
- Wolfinger, Russ. (1993). Covariance structure selection in general mixed models. Communications in Statistics-simulation and Computation - COMMUN STATIST-SIMULAT COMPUT. 22. 1079-1106. 10.1080/03610919308813143.
- R Henderson, C. (1984). Application of Linear Models in Animal Breeding.
- Giesbrecht, F.G. & Burns, Joseph. (1985). Two-Stage Analysis Based on a Mixed Model: Large-Sample Asymptotic Theory and Small-Sample Simulation Results. Biometrics. 41. 10.2307/2530872.
- G. Kenward, Michael & Roger, James. (1997). Small Sample Inference for Fixed Effects From Restricted Maximum Likelihood. Biometrics. 53. 983-97. 10.2307/2533558.
- Wright, Stephen, and Jorge Nocedal (2006) "Numerical optimization." Springer

Author: Vladimir Arnautov aka PharmCat
Copyright © 2019 Vladimir Arnautov aka PharmCat <mail@pharmcat.net>
