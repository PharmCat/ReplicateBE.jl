task = ReplicateBE.RandRBEDS(;n=24,
sequence=[1,1], design = ["T" "R" "T" "R"; "R" "T" "R" "T"],
inter=[0.05, 0.04, 0.6], intra=[0.02, 0.02],
intercept = 1.0, seqcoef = [0.0, 0.0], periodcoef = [0.0, 0.0, 0.0, 0.0],
formcoef = [0.0, log(0.8)], seed = 2283587414, dropobs = 10)

df = ReplicateBE.randrbeds(task)

be  = ReplicateBE.rbe!(df, dvar = :var, subject = :subject, formulation = :formulation, period = :period, sequence = :sequence)
