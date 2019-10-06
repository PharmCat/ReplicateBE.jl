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

  SAS procedures for generated datasets can be found in [*validation/sas/*](https://github.com/PharmCat/ReplicateBE.jl/tree/master/validation/sas) folder.

  All validation datasets include in package test procedure.

  *SAS WARNING: Stopped because of infinite likelihood.
  
  ## Simulation study
  
  Following simulation was performed:
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
