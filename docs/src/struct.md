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
