
## Syntax

```
rbe(df; dvar::Symbol,
    subject::Symbol,
    formulation::Symbol,
    period::Symbol,
    sequence::Symbol,
    g_tol::Float64 = 1e-8, x_tol::Float64 = 0.0, f_tol::Float64 = 0.0, iterations::Int = 100,
    store_trace = false, extended_trace = false, show_trace = false,
    memopt = true)
```

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
