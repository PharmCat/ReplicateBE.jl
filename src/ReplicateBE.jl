# ReplicateBE
# Author: Vladimir Arnautov aka PharmCat
# Copyright © 2019 Vladimir Arnautov aka PharmCat <mail@pharmcat.net>
# Licence: GNU General Public License v3.0

module ReplicateBE

using DataFrames, Distributions, StatsModels, StatsBase, ForwardDiff, LinearAlgebra, Random, Optim, LineSearches, CategoricalArrays, Printf, StaticArrays

    export rbe, rbe!, reml2, nobs, coef, stderror, dof, coefnum, fixed, theta, typeiii, design, show, confint, contrast, estimate, optstat, randrbeds, randrbetask
    import Base.show
    import StatsBase.confint, StatsBase.coef, StatsBase.nobs, StatsBase.dof, StatsBase.stderror, StatsBase.coeftable, StatsBase.vcov
    import Statistics.var

const LOG2PI = log(2π)
const MEMOPT = true

include("rbesettings.jl")
include("rbetable.jl")
include("memalloc.jl")
include("design.jl")
include("randrbeds.jl")
include("deprecated.jl")
include("rbe.jl")
include("utils.jl")
include("generalfunc.jl")
include("algebra.jl")

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
end # module
