
## Basic using

```
using CSV, DataFrames, ReplicateBE

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

be = ReplicateBE.rbe(df, dvar = :var, subject = :subject, formulation = :formulation, period = :period, sequence = :sequence)

#Get fixed effect object

ReplicateBE.fixed(be)

#Get Type III object

ReplicateBE.typeiii(be)

#Get model coefficients

coef(be)

#Get Standard Error for coefficients

ReplicateBE.coefse(be)

#Get confidence intervals

ReplicateBE.confint(be, 0.1, expci = false, inv = false)

#Get -2 REML for model

ReplicateBE.reml2(be)

#Design information

ReplicateBE.design(be)
```
