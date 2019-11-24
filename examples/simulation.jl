using ReplicateBE

task = ReplicateBE.RandRBEDS(;n=24,
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

task = ReplicateBE.RandRBEDS(;n=24,
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

task = ReplicateBE.RandRBEDS(;n=24,
sequence=[1,2], design = ["T" "R" "T" "R"; "R" "T" "R" "T"],
inter=[0.05, 0.04, 0.3], intra=[0.02, 0.04],
intercept = 1.0, seqcoef = [0.0, 0.0], periodcoef = [0.0, 0.0, 0.0, 0.0],
formcoef = [0.0, log(0.8)], seed = 0, dropobs = 10)

result =  ReplicateBE.simulation(task; num = 10, seed = 6564683774737, verbose = true)
