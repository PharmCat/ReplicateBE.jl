# ReplicateBE.jl
Mixed model solution for replicate designed bioequivalence study.

Install
```
using Pkg; Pkg.clone("https://github.com/PharmCat/ReplicateBE.jl.git")
```

Using:
```
be = ReplicateBE.rbe(df, var = :var, subject = :subject, formulation = :formulation, period = :period, sequence = :sequence);
```
Where:

- var::Symbol - dependent variable;
- subject::Symbol - subject;
- formulation::Symbol - formulation/drug;
- period::Symbol - study period;
- sequence::Symbol - sequence;

Get resuls:

```
be.β            #For β
be.se           #For SE
be.reml         #REML value
```

Other:

```
struct RBE
    model
    factors
    β
    θ
    reml
    se
    F
    DF
    R
    V
    G
    A
    H
    detH
end
```

Author: Vladimir Arnautov aka PharmCat

Copyright © 2019 Vladimir Arnautov aka PharmCat <mail@pharmcat.net>
