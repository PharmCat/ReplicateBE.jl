# ReplicateBE.jl Documentation

```@meta
CurrentModule = ReplicateBE
```

*ReplicateBE* is a Julia package providing mixed model solution for replicate designed bioequivalence study. This can be used to obtained results with methods C (random effects with interaction), given by the EMA in [Annex I](https://www.ema.europa.eu/en/documents/other/31-annex-i-statistical-analysis-methods-compatible-ema-bioequivalence-guideline_en.pdf "EMA/582648/2016, 21 September 2016"). Statistical model formed with accordance [FDA Guidance for Industry: Statistical Approaches to Establishing Bioequivalence](https://www.fda.gov/media/70958/download), APPENDIX F.


## Installation

Install:
```
using Pkg; Pkg.add("ReplicateBE")
```

or:
```
using Pkg; Pkg.clone("https://github.com/PharmCat/ReplicateBE.jl.git")
```


```@contents
Pages = [
        "examples.md",
        "syntax.md",
        "details.md",
        "testval.md",
        "struct.md",
        "api.md"]
Depth = 3
```
