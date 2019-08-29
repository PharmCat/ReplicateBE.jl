# ReplicateBE.jl

Alpha version! This program comes with absolutely no warranty. No liability is accepted for any loss and risk to public health resulting from use of this software.

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
be.β            #For β vector
be.θ            #For variance components vector
be.se           #For SE vector
be.reml         #REML value
be.df           #Degree of freedom
be.f            #F Statistics vector
be.G            #G matrix
```

Other:

```
struct RBE
    model::ModelFrame               #Model frame
    rmodel::ModelFrame              #Random effect model
    factors::Array{Symbol, 1}       #Factor list
    β::Array{Float64, 1}            #β coefficients (fixed effect)
    θ0::Array{Float64, 1}           #Initial variance paramethers
    θ::Array{Float64, 1}            #Final variance paramethers
    reml::Float64                   #-2REML
    se::Array{Float64, 1}           #SE for each β level
    f::Array{Float64, 1}            #F for each β level
    df::Array{Float64, 1}           #DF (degree of freedom) for each β level (Satterthwaite)
    df2::Float64                    #DF N / pn - sn
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

### ReplicateBE.reml2

```
reml2(obj::RBE, θ::Array{Float64, 1})
```
Return -2REML for vector θ.

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
