using ReplicateBE


task = ReplicateBE.RandRBEDS(;n=24,
sequence=[1,1], design = ["T" "R" "T" "R"; "R" "T" "R" "T"],
inter=[0.5, 0.4, 0.5], intra=[0.1, 0.2],
intercept = 1.0, seqcoef = [0.0, 0.0], periodcoef = [0.0, 0.0, 0.0, 0.0],
formcoef = [0.0, log(0.8)], seed = 10001)


pow =  ReplicateBE.simulation(task; num = 200, seed = 123)
#df = ReplicateBE.randrbeds(task)
