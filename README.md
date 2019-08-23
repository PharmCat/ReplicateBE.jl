# ReplicateBE.jl
Mixed model solution for replicate designed bioequivalence study.

[![Build Status](https://api.travis-ci.com/PharmCat/ReplicateBE.jl.svg?branch=master)](https://travis-ci.com/PharmCat/ReplicateBE.jl)


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

Get results (as vectors):

```
be.β            #For β
be.se           #For SE
be.reml         #REML value
be.DF           #Degree of freedom
be.F            #F Statistics
```

Other:

```
struct RBE
    model::ModelFrame
    factors::Array{Symbol, 1}
    β::Array{Float64, 1}
    θ::Array{Float64, 1}
    reml::Float64
    se::Array{Float64, 1}
    F::Array{Float64, 1}
    DF::Array{Float64, 1}
    R::Array{Matrix{Float64},1}
    V::Array{Matrix{Float64},1}
    G::Matrix{Float64}
    A::Matrix{Float64}
    H::Matrix{Float64}
    detH::Float64
    optim::Optim.MultivariateOptimizationResults
end
```
# Methods

StatsBase.confint(obj::RBE, alpha::Float64; expci::Bool = false, inv::Bool = false)

* obj::RBE - bioequivalence struct;
* alpha::Float64 - alpha;
* expci::Bool - exp(ci)
* inv::Bool - β = -β



Author: Vladimir Arnautov aka PharmCat

Copyright © 2019 Vladimir Arnautov aka PharmCat <mail@pharmcat.net>
