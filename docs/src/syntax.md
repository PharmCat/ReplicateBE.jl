## Syntax

### Simple syntax

Simple syntax can be used in general purpose. Can modify df DataFrame.

```julia
rbe!(df; dvar::Symbol, subject::Symbol, formulation::Symbol, period::Symbol, sequence::Symbol)
```

### Full syntax (Not modifying syntax)

```@docs
ReplicateBE.rbe
```

### Modifying syntax

```@docs
ReplicateBE.rbe!
```

### Random dataset

Random dataset function is made for generation validation datasets and simulation data.

```@docs
ReplicateBE.randrbeds
```

### Simulation

Simulation based on work with task - RandRBEDS object.
Purpose of simulation to count successful BE outcomes for power or alpha calculation.
Also simulation used to test package stability.

```@docs
ReplicateBE.simulation
```
