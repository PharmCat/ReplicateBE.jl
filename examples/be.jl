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
