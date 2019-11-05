using ReplicateBE

task = ReplicateBE.RandRBEDS(;n=24,
sequence=[1,1], design = ["T" "R" "T" "R"; "R" "T" "R" "T"],
inter=[0.05, 0.04, 0.6], intra=[0.02, 0.02],
intercept = 1.0, seqcoef = [0.0, 0.0], periodcoef = [0.0, 0.0, 0.0, 0.0],
formcoef = [0.0, log(0.8)], seed = 0, dropobs = 10)

result =  ReplicateBE.simulation(task; num = 200, seed = 123)
