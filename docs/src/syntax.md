## Syntax

### Simple syntax

Simple syntax can be used in general purpose.

```julia
rbe!(df; dvar::Symbol, subject::Symbol, formulation::Symbol, period::Symbol, sequence::Symbol)
```

### Full syntax

With keywords you can define additional options.

```julia
rbe!(df; dvar::Symbol,
    subject::Symbol,
    formulation::Symbol,
    period::Symbol,
    sequence::Symbol,
    g_tol::Float64 = 1e-8, x_tol::Float64 = 0.0, f_tol::Float64 = 0.0, iterations::Int = 100,
    store_trace = false, extended_trace = false, show_trace = false,
    memopt = true)
```

### Arguments and keywords

- df::DataFrame - DataFrame with data
- dvar::Symbol - variable
- subject::Symbol - Subject column
- formulation::Symbol - Formulation clumn
- period::Symbol - Period column
- sequence::Symbol - Sequence column
- g_tol
- x_tol
- f_tol
- iterations:: - Maximum iterations
- store_trace
- extended_trace
- show_trace
- memopt::Bool - memory optimization, can increase performance  

### Not modifying syntax

If you use ```rbe()``` function no data transformations done with ```df```such as ```categorical!()``` and ```sort!()```.

```julia
rbe(df; dvar::Symbol, subject::Symbol, formulation::Symbol, period::Symbol, sequence::Symbol)
```

### Random dataset

Random dataset function is made for generation validation datasets and simulation data. 

```
randrbeds(;n=24, sequence=[1,1],
    design = ["T" "R" "T" "R"; "R" "T" "R" "T"],
    inter=[0.5, 0.4, 0.9], intra=[0.1, 0.2],
    intercept = 0, seqcoef = [0.0, 0.0], periodcoef = [0.0, 0.0, 0.0, 0.0], formcoef = [0.0, 0.0], seed = 0)
```
Where:

 - n - Subject number;
 - sequence - sequence subject distribution, [1,1] is equal 1:1, [1,2] - 1:2, [2,3,5] - 2:3:5 ets.;
 - design - design matrix (sXp, where s - number of sequences, p - number of period), cells contains formulation label;
 - inter - Inter-subject variation vector for G matrix: [σ₁, σ₂, ρ], where σ₁, σ₂ - formulation inter-subject variance,  ρ - covariance coefficient;
 - intra - Intra-subject variation vector for R matrix:[σ₁, σ₂], where σ₁, σ₂ - formulation intra-subject variance;
 - intercept - model intercept value;
 - seqcoef - model sequence coefficient values (additive): length = s (number of sequences);
 - periodcoef - model period coefficient values (additive): length = p (number of periods);
 - formcoef - model formulation coefficient values (additive): length = 2;
