## Examples

### Bioequivalence example

```julia
using CSV, DataFrames, ReplicateBE, StatsBase

#Load Dataframe...

df=CSV.read(IOBuffer("""subject,sequence,period,formulation,var
1,1,1,1,1.0
1,1,2,2,1.1
1,1,3,1,1.2
1,1,4,2,1.3
2,1,1,1,2.0
2,1,2,2,2.1
2,1,3,1,2.4
2,1,4,2,2.2
3,2,1,2,1.3
3,2,2,1,1.5
3,2,3,2,1.6
3,2,4,1,1.4
4,2,1,2,1.5
4,2,2,1,1.7
4,2,3,2,1.3
4,2,4,1,1.4
5,2,1,2,1.5
5,2,2,1,1.7
5,2,3,2,1.2
5,2,4,1,1.8""")) |> DataFrame

#Execute BE

be = rbe!(df, dvar = :var, subject = :subject, formulation = :formulation, period = :period, sequence = :sequence)

#Get fixed effect object

fixed(be)

#Get Type III object

ypeiii(be)

#Get model coefficients

coef(be)

#Get Standard Error for coefficients

stderror(be)

#Get confidence intervals for all coefficients

confint(be, 0.1, expci = false, inv = false)

#Get -2 REML for model

reml2(be)

#Design information

design(be)
```

### Random dataset generation

```julia
using ReplicateBE

task = ReplicateBE.RandRBEDS(;n=24,
sequence=[1,1], design = ["T" "R" "T" "R"; "R" "T" "R" "T"],
inter=[0.05, 0.04, 0.6], intra=[0.02, 0.02],
intercept = 1.0, seqcoef = [0.0, 0.0], periodcoef = [0.0, 0.0, 0.0, 0.0],
formcoef = [0.0, log(0.8)], seed = 1234, dropobs = 10)

df = ReplicateBE.randrbeds(task)

#or

df = ReplicateBE.randrbeds(;n=24,
sequence=[1,1], design = ["T" "R" "T" "R"; "R" "T" "R" "T"],
inter=[0.05, 0.04, 0.6], intra=[0.02, 0.02],
intercept = 1.0, seqcoef = [0.0, 0.0], periodcoef = [0.0, 0.0, 0.0, 0.0],
formcoef = [0.0, log(0.8)], seed = 1234, dropobs = 10)

```

### Simulation

```julia

using ReplicateBE

task = ReplicateBE.RandRBEDS(;n=24,
sequence=[1,1], design = ["T" "R" "T" "R"; "R" "T" "R" "T"],
inter=[0.05, 0.04, 0.6], intra=[0.02, 0.02],
intercept = 1.0, seqcoef = [0.0, 0.0], periodcoef = [0.0, 0.0, 0.0, 0.0],
formcoef = [0.0, log(0.8)], seed = 0, dropobs = 10)

result =  ReplicateBE.simulation(task; num = 200, seed = 123)

```
