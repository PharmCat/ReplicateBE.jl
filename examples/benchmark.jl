using BenchmarkTools, CSV, DataFrames
using ReplicateBE

datapath = joinpath(dirname(pathof(ReplicateBE)))*"/../test/csv/df6.csv"

df6     = CSV.File(datapath) |> DataFrame

#@benchmark be = ReplicateBE.rbe!(df6, dvar = :var, subject = :subject, formulation = :formulation, period = :period, sequence = :sequence, g_tol = 1e-10)
@benchmark be = ReplicateBE.rbe!($df6, dvar = :var, subject = :subject, formulation = :formulation, period = :period, sequence = :sequence, g_tol = 1e-10) seconds = 25 samples = 200 evals = 5

#=
julia> @benchmark be = ReplicateBE.rbe!(datadf, dvar = :var, subject = :subject, formulation = :formulation, period = :period, sequence = :sequence, g_tol = 1e-10)
BenchmarkTools.Trial:
  memory estimate:  6.15 MiB
  allocs estimate:  66530
  --------------
  minimum time:     24.307 ms (0.00% GC)
  median time:      27.732 ms (0.00% GC)
  mean time:        29.258 ms (3.78% GC)
  maximum time:     47.647 ms (0.00% GC)
  --------------
  samples:          172
  evals/sample:     1

julia>
=#

#v1.5
#=
julia> @benchmark be = ReplicateBE.rbe!(df6, dvar = :var, subject = :subject, formulation = :formulation, period = :period, sequence = :sequence, g_tol = 1e-10)
BenchmarkTools.Trial:
  memory estimate:  5.40 MiB
  allocs estimate:  48564
  --------------
  minimum time:     16.528 ms (0.00% GC)
  median time:      18.803 ms (0.00% GC)
  mean time:        20.431 ms (4.65% GC)
  maximum time:     53.102 ms (50.43% GC)
  --------------
  samples:          245
  evals/sample:     1
=#

#v1.0.11 Julia 1.4
#=
julia> @benchmark be = ReplicateBE.rbe!(df6, dvar = :var, subject = :subject, formulation = :formulation, period = :period, sequence = :sequence, g_tol
= 1e-10)
BenchmarkTools.Trial:
  memory estimate:  4.34 MiB
  allocs estimate:  46670
  --------------
  minimum time:     10.112 ms (0.00% GC)
  median time:      18.371 ms (0.00% GC)
  mean time:        19.176 ms (2.40% GC)
  maximum time:     83.157 ms (74.30% GC)
  --------------
  samples:          261
  evals/sample:     1
=#

#v1.0.11d Julia 1.5.1
#=
julia> @benchmark be = ReplicateBE.rbe!(df6, dvar = :var, subject = :subject, formulation = :formulation, period = :period, sequence = :sequence, g_tol
= 1e-10)
BenchmarkTools.Trial:
  memory estimate:  4.48 MiB
  allocs estimate:  43228
  --------------
  minimum time:     17.604 ms (0.00% GC)
  median time:      19.636 ms (0.00% GC)
  mean time:        20.833 ms (2.70% GC)
  maximum time:     47.499 ms (0.00% GC)
  --------------
  samples:          240
  evals/sample:     1
=#

#=
julia> @benchmark be = ReplicateBE.rbe!($df6, dvar = :var, subject = :subject, formulation = :formulation, period = :period, sequence = :sequence, g_tol = 1e-10) seconds = 25 samples = 200 evals = 5
BenchmarkTools.Trial:
  memory estimate:  4.47 MiB
  allocs estimate:  43136
  --------------
  minimum time:     17.257 ms (0.00% GC)
  median time:      19.049 ms (0.00% GC)
  mean time:        19.897 ms (3.34% GC)
  maximum time:     29.502 ms (32.83% GC)
  --------------
  samples:          200
  evals/sample:     5
=#
