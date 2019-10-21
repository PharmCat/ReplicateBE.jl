# Package validation

Validation program include public datasets and generated datasets. All public datasets include in test/testdata.jl.

Generated datasets made with *randrbeds* function with fixed seed (could be reproduced at any time) and include following disignes:

### 24 subjects, balanced

  * 1  TRTR/RTRT
  * 2  TRRT/RTTR
  * 3  TTRR/RRTT
  * 4  TRTR/RTRT/TRRT/RTTR
  * 5  TRRT/RTTR/TTRR/RRTT
  * 6  TRTR/RTRT/TTRR/RRTT
  * 7  TRT/RTR
  * 8  TRR/RTT
  * 9  TR/RT/TT/RR
  * 10 TRR/RTR/RRT
  * 11 TRR/RTR*

### 48 subjects, unbalanced, 20 dropped observations

    * 12  TRTR/RTRT
    * 13  TRRT/RTTR
    * 14  TTRR/RRTT
    * 15  TRTR/RTRT/TRRT/RTTR
    * 16  TRRT/RTTR/TTRR/RRTT
    * 17  TRTR/RTRT/TTRR/RRTT
    * 18  TRT/RTR
    * 19  TRR/RTT
    * 20  TR/RT/TT/RR
    * 21 TRR/RTR/RRT
    * 22 TRR/RTR*

### Special cases

    * 101 SP1: TRTR/RTRT 1024 subjects, 2000 dropped observations
    * 102 SP2: TRT/RTR 4096 subjects, 2000 dropped observations (total 10288 observations)

SP1 output:
```
rds = ReplicateBE.randrbeds(;n=1024, sequence=[1,2], design = ["T" "R" "T" "R"; "R" "T" "R" "T"], inter=[0.5, 0.4, 0.9], intra=[0.1, 0.2], intercept = 1.0, seqcoef = [0.0, 0.0], periodcoef = [0.0, 0.0, 0.0, 0.0], formcoef = [0.0, 0.1], dropobs = 2000, seed = 10101)
@time be = ReplicateBE.rbe(rds, dvar = :var, subject = :subject, formulation = :formulation, period = :period, sequence = :sequence)
  6.107196 seconds (14.51 M allocations: 981.933 MiB, 8.28% gc time)
Bioequivalence Linear Mixed Effect Model (status: converged)

-2REML: 4032.92    REML: -2016.46

Fixed effect:
─────────────────────────────────────────────────────────────────────────────────────────────
Effect           Value         SE          F           DF        t           P|t|
─────────────────────────────────────────────────────────────────────────────────────────────
(Intercept)      1.17091       0.0342916   1165.94     1215.96   34.1458     9.58974e-180*
sequence: TRTR   -0.0230415    0.0497151   0.214805    976.825   -0.463471   0.64313
period: 2        -0.0118358    0.0300503   0.15513     1397.78   -0.393865   0.69374
period: 3        -0.00522645   0.0288399   0.0328418   1139.11   -0.181223   0.856225
period: 4        -0.0271755    0.0300382   0.81848     1398.08   -0.904699   0.365781
formulation: T   -0.102651     0.0231198   19.7131     1445.78   -4.43995    9.68179e-6*
─────────────────────────────────────────────────────────────────────────────────────────────
Intra-individual variation:
formulation: R   0.202459
formulation: T   0.105328

Inter-individual variation:
formulation: R   0.396337
formulation: T   0.52727
ρ:               0.933991   Cov: 0.426964

Confidence intervals(90%):
formulation: R / formulation: T
86.8747 - 93.7445 (%)
formulation: T / formulation: R
106.673 - 115.1083 (%)
```

SP1 output:    
```
rds = ReplicateBE.randrbeds(;n=4096, sequence=[1,4], design = ["T" "R" "T"; "R" "T" "R"], inter=[0.5, 0.4, 0.9], intra=[0.1, 0.9], intercept = 1.0, seqcoef = [10.0, 0.0], periodcoef = [0.0, 0.0, 0.0], formcoef = [0.0, 1.0], dropobs = 2000, seed = 10102)
be = ReplicateBE.rbe(rds, dvar = :var, subject = :subject, formulation = :formulation, period = :period, sequence = :sequence)
312.655087 seconds (269.38 M allocations: 14.273 GiB, 1.47% gc time)
Bioequivalence Linear Mixed Effect Model (status: converged)

-2REML: 26633.2    REML: -13316.6

Fixed effect:
───────────────────────────────────────────────────────────────────────────────────────
Effect           Value         SE          F          DF        t           P|t|
───────────────────────────────────────────────────────────────────────────────────────
(Intercept)      1.99206       0.0184652   11638.4    4779.03   107.882     0.0*
sequence: TRT    10.0134       0.0322104   96643.1    5028.62   310.875     0.0*
period: 2        0.0193164     0.0223827   0.744777   6811.0    0.863005    0.388165
period: 3        -0.00501682   0.0152331   0.108463   1577.27   -0.329338   0.741944
formulation: T   -0.995365     0.0210042   2245.71    6298.6    -47.3889    0.0*
Intra-individual variation:
formulation: R   0.885914
formulation: T   0.106905

Inter-individual variation:
formulation: R   0.380248
formulation: T   0.491359
ρ:               0.902941   Cov: 0.390294

Confidence intervals(90%):
formulation: R / formulation: T
35.7036 - 38.2583 (%)
formulation: T / formulation: R
261.3815 - 280.0838 (%)
```

  SAS procedures for generated datasets can be found in [*validation/sas/*](https://github.com/PharmCat/ReplicateBE.jl/tree/master/validation/sas) folder.

  All validation datasets (except special - SPX) include in package test procedure.

  *SAS WARNING*: Stopped because of infinite likelihood.

  ## Simulation study

  Following simulation was performed for version v0.1.4:
  ```
  using Distributions, ReplicateBE
  function simulation(num)
    n   = 0
    err = 0
    cnt = 0
    b = l = log(0.8)
    u = log(1.25)
    println("Start...")
    for i = 1:num
        try
            rds   = ReplicateBE.randrbeds(;n=24, sequence=[1,1], design = ["T" "R" "T" "R"; "R" "T" "R" "T"], inter=[0.2, 0.2, 0.5], intra=[0.05, 0.05], intercept = 1.0, seqcoef = [0.0, 0.0], periodcoef = [0.0, 0.0, 0.0, 0.0], formcoef = [b, 0.0])
            be    = ReplicateBE.rbe(rds, dvar = :var, subject = :subject, formulation = :formulation, period = :period, sequence = :sequence, memopt = false)
            q     = quantile(TDist(be.fixed.df[6]), 0.95)
            ll    = be.fixed.est[6] - q*be.fixed.se[6]
            ul    = be.fixed.est[6] + q*be.fixed.se[6]
            if ll > l && ul < u cnt += 1 end
            if n > 1000
                println("Iteration: $i")
                println("Mem: $(Sys.free_memory()/2^20)")
                println("Pow: $(cnt/i)")
                println("-------------------------------")
                n = 0
            end
            n += 1
        catch
            err += 1
            println("Error $(err)!")
        end
    end
    return cnt/num
end

simulation(100000)
```

 Rusults:
```
 7104.944474 seconds (9.51 G allocations: 2.230 TiB, 44.90% gc time)
0.04939
```
 Cofinence interval (95%) for power: 0.048047 - 0.050733. No statistically significant difference found.
