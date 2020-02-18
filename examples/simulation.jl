using ReplicateBE

task = ReplicateBE.randrbetask(;n=24,
sequence=[1,2], design = ["T" "R" "T" "R"; "R" "T" "R" "T"],
inter=[0.05, 0.04, 0.6], intra=[0.02, 0.04],
intercept = 1.0, seqcoef = [0.0, 0.0], periodcoef = [0.0, 0.0, 0.0, 0.0],
formcoef = [0.0, log(0.8)], seed = 0, dropobs = 0)


@time result =  ReplicateBE.simulation(task; num = 100000, seed = 730795390628834530, verbose = true)
#=
1837.914923 seconds (3.99 G allocations: 1.092 TiB, 9.91% gc time)
Seed: 730795390628834530
Number: 100000
Result: 0.05078660225829358
=#

task = ReplicateBE.randrbetask(;n=24,
sequence=[1,2], design = ["T" "R" "T"; "R" "T" "R"],
inter=[0.05, 0.04, 0.4], intra=[0.02, 0.04],
intercept = 1.0, seqcoef = [0.0, 0.0], periodcoef = [0.0, 0.0, 0.0],
formcoef = [0.0, log(0.8)], seed = 0, dropobs = 0)
@time result =  ReplicateBE.simulation(task; num = 100000, seed = 730795390628834530, verbose = true)

#=
1638.566348 seconds (3.77 G allocations: 907.478 GiB, 9.77% gc time)
Seed: 730795390628834530
Number: 100000
Result: 0.04997148203368122
#E63
=#

task = ReplicateBE.randrbetask(;n=24,
sequence=[1,2], design = ["T" "R" "T" "R"; "R" "T" "R" "T"],
inter=[0.05, 0.04, 0.3], intra=[0.02, 0.04],
intercept = 1.0, seqcoef = [0.0, 0.0], periodcoef = [0.0, 0.0, 0.0, 0.0],
formcoef = [0.0, log(0.8)], seed = 0, dropobs = 10)

result =  ReplicateBE.simulation(task; num = 2000, seed = 6564683774737, verbose = true)


task = ReplicateBE.randrbetask(;n=24,
sequence=[1,2], design = ["T" "R" "T" "R"; "R" "T" "R" "T"],
inter=[0.01, 0.01, 0.9], intra=[0.09, 0.09],
intercept = 1.0, seqcoef = [0.0, 0.0], periodcoef = [0.0, 0.0, 0.0, 0.0],
formcoef = [log(0.95), 0.0], seed = 0)

rds    =  ReplicateBE.randrbeds(task)
be = ReplicateBE.rbe!(rds, dvar = :var, subject = :subject, formulation = :formulation, period = :period, sequence = :sequence)

result =  ReplicateBE.simulation(task; num = 1002, verbose = true, rsabe = false)



task = ReplicateBE.randrbetask(;n=24,
sequence=[1,2], design = ["T" "R" "T" "R"; "R" "T" "R" "T"],
inter=2000.0, intra=[0.09, 0.09],
intercept = 1.0, seqcoef = [0.0, 0.0], periodcoef = [0.0, 0.0, 0.0, 0.0],
formcoef = [log(0.95), 0.0], seed = 0)

rds    =  ReplicateBE.randrbeds(task)
result =  ReplicateBE.simulation(task; num = 1002, verbose = true, rsabe = true)




task = ReplicateBE.randrbetask(;n=24,
sequence=[1,1], design = ["T" "R" "T" "R"; "R" "T" "R" "T"],
inter=[0.04, 0.04, 1.0], intra=[0.02, 0.02],
intercept = 1.0, seqcoef = [0.0, 0.0], periodcoef = [0.0, 0.0, 0.0, 0.0],
formcoef = [0.0, log(0.8)], seed = 0, dropobs = 0)

out = zeros(Float64, 0)
function simfunc!(out, be)
    push!(out, ReplicateBE.theta(be)[5])
end
result =  ReplicateBE.simulation!(task, out, simfunc!; num = 10, seed = 10, verbose = false)
m      = mean(result)

using Plots

histogram(result)
