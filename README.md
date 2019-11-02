**Public validated Beta version!** This program comes with absolutely no warranty. No liability is accepted for any loss and risk to public health resulting from use of this software.

<p align="center">
  <img src="https://github.com/PharmCat/ReplicateBE.jl/blob/master/docs/ReplicateBE-LogoNoSpace.png">
</p>

Mixed model solution for replicate designed bioequivalence study. This can be used to obtained results with methods C (random effects with interaction), given by the EMA in [Annex I](https://www.ema.europa.eu/en/documents/other/31-annex-i-statistical-analysis-methods-compatible-ema-bioequivalence-guideline_en.pdf "EMA/582648/2016, 21 September 2016"). Statistical model formed with accordance [FDA Guidance for Industry: Statistical Approaches to Establishing Bioequivalence](https://www.fda.gov/media/70958/download), APPENDIX F.

[![GitHub version](https://badge.fury.io/gh/PharmCat%2FReplicateBE.jl.svg)](https://badge.fury.io/gh/PharmCat%2FReplicateBE.jl)
[![Build Status](https://api.travis-ci.com/PharmCat/ReplicateBE.jl.svg?branch=master)](https://travis-ci.com/PharmCat/ReplicateBE.jl)
[![codecov](https://codecov.io/gh/PharmCat/ReplicateBE.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/PharmCat/ReplicateBE.jl)
[![Latest docs](https://img.shields.io/badge/docs-latest-blue.svg)](https://pharmcat.github.io/ReplicateBE.jl/build/)
[![doi](https://img.shields.io/badge/doi-10.13140%2FRG.2.2.27418.39363-blue)](https://doi.org/10.13140/RG.2.2.27418.39363)

http://doi.org/10.13140/RG.2.2.27418.39363
Install:
```
using Pkg; Pkg.add("ReplicateBE")
```
or:
```
using Pkg; Pkg.clone("https://github.com/PharmCat/ReplicateBE.jl.git")
```

Using:
```
using ReplicateBE
be = ReplicateBE.rbe!(df, dvar = :var, subject = :subject, formulation = :formulation, period = :period, sequence = :sequence);
ci = confint(be, 0.1)
```
Where:

- dvar::Symbol - dependent variable;
- subject::Symbol - subject;
- formulation::Symbol - formulation/drug;
- period::Symbol - study period;
- sequence::Symbol - sequence.

How to get results?
```
#Fixed effect table:
fixed(be)

#Type III table
typeiii(be)
```

Output example:
```
Bioequivalence Linear Mixed Effect Model

-2REML: 164.613    REML: -82.3067

Fixed effect:
─────────────────────────────────────────────────────────────────────────────────────────
Effect           Value        SE         F           DF        t           P|t|
─────────────────────────────────────────────────────────────────────────────────────────
(Intercept)      1.27158      0.191994   43.8645     26.6243   6.62303     4.47193e-7*
sequence: TRTR   -0.114458    0.242518   0.222743    30.4503   -0.471956   0.640323
period: 2        0.035864     0.139526   0.0660707   83.1868   0.257042    0.797781
period: 3        0.0168505    0.100284   0.0282329   41.619    0.168027    0.867376      
period: 4        0.0552216    0.139526   0.156642    83.1868   0.39578     0.693281
formulation: T   -0.0165797   0.120162   0.0190378   62.661    -0.137978   0.890701      
─────────────────────────────────────────────────────────────────────────────────────────
Intra-individual variation:
formulation: R  0.187712   CVᵂ: 0.454407
formulation: T  0.0889289   CVᵂ: 0.304964

Inter-individual variation:
formulation: R  0.318313
formulation: T  0.422413
Cov:  0.266255

Confidence intervals(90%):
formulation: R / formulation: T
80.4773 - 120.2059 (%)
formulation: T / formulation: R
83.1906 - 124.2586 (%)
```

# Validation

Validation information: [here](https://pharmcat.github.io/ReplicateBE.jl/build/testval/)

# Basic methods

All API docs see [here](https://pharmcat.github.io/ReplicateBE.jl/build/api/).

# Random Dataset

Random dataset function is made for generation validation datasets and simulation data.  

```
randrbeds(;n=24, sequence=[1,1],
    design = ["T" "R" "T" "R"; "R" "T" "R" "T"],
    inter=[0.5, 0.4, 0.9], intra=[0.1, 0.2],
    intercept = 0, seqcoef = [0.0, 0.0], periodcoef = [0.0, 0.0, 0.0, 0.0], formcoef = [0.0, 0.0], seed = 0)
```
Generate DataFrame with random multivariate data. Where:

 - n - Subject number;
 - sequence - sequence subject distribution, [1,1] is equal 1:1, [1,2] - 1:2, [2,3,5] - 2:3:5 ets.;
 - design - design matrix (sXp, where s - number of sequences, p - number of period), cells contains formulation label;
 - inter - Inter-subject variation vector for G matrix: [σ₁, σ₂, ρ], where σ₁, σ₂ - formulation inter-subject variance,  ρ - covariance coefficient;
 - intra - Intra-subject variation vector for R matrix:[σ₁, σ₂], where σ₁, σ₂ - formulation intra-subject variance;
 - intercept - model intercept value;
 - seqcoef - model sequence coefficient values (additive): length = s (number of sequences);
 - periodcoef - model period coefficient values (additive): length = p (number of periods);
 - formcoef - model formulation coefficient values (additive): length = 2;

## Structures

Struct information see [here](https://pharmcat.github.io/ReplicateBE.jl/build/struct/).

## Acknowledgments

Best acknowledgments to D.Sc. in Physical and Mathematical Sciences Anastasia Shitova <a.shitova@qayar.ru> for support, datasets and testing procedures.

## References

* FDA Guidance for Industry: Statistical Approaches to Establishing Bioequivalence, 2001
* Fletcher, Roger (1987), Practical methods of optimization (2nd ed.), New York: John Wiley & Sons, ISBN 978-0-471-91547-8
* Giesbrecht, F. G., and Burns, J. C. (1985), "Two-Stage Analysis Based on a Mixed Model: Large-sample Asymptotic Theory and Small-Sample Simulation Results," Biometrics, 41, 853-862.
* Gurka, Matthew. (2006). Selecting the Best Linear Mixed Model under REML. The American Statistician. 60. 19-26. 10.1198/000313006X90396.
* Henderson, C. R., et al. “The Estimation of Environmental and Genetic Trends from Records Subject to Culling.” Biometrics, vol. 15, no. 2, 1959, pp. 192–218. JSTOR, www.jstor.org/stable/2527669.
* Hrong-Tai Fai & Cornelius (1996) Approximate F-tests of multiple degree of freedom hypotheses in generalized least squares analyses of unbalanced split-plot experiments, Journal of Statistical Computation and Simulation, 54:4, 363-378, DOI: 10.1080/00949659608811740
* Jennrich, R., & Schluchter, M. (1986). Unbalanced Repeated-Measures Models with Structured Covariance Matrices. Biometrics, 42(4), 805-820. doi:10.2307/2530695
* Laird, Nan M., and James H. Ware. “Random-Effects Models for Longitudinal Data.” Biometrics, vol. 38, no. 4, 1982, pp. 963–974. JSTOR, www.jstor.org/stable/2529876.
* Lindstrom & J.; Bates, M. (1988). Newton—Raphson and EM Algorithms for Linear Mixed-Effects Models for Repeated-Measures Data. Journal of the American Statistical Association. 83. 1014. 10.1080/01621459.1988.10478693.
* Mogensen et al., (2018). Optim: A mathematical optimization package for Julia. Journal of Open Source Software, 3(24), 615,doi: 10.21105/joss.00615
* Patterson, S. D. and Jones, B. (2002), Bioequivalence and the pharmaceutical industry. Pharmaceut. Statist., 1: 83-95. doi:10.1002/pst.15
* Revels, Jarrett & Lubin, Miles & Papamarkou, Theodore. (2016). Forward-Mode Automatic Differentiation in Julia.
* Schaalje GB, McBride JB, Fellingham GW. Adequacy of approximations to distributions of test statistics in complex mixed linear models. J Agric Biol Environ Stat. 2002;7:512–24.
* Van Peer, A. (2010), Variability and Impact on Design of Bioequivalence Studies. Basic & Clinical Pharmacology & Toxicology, 106: 146-153. doi:10.1111/j.1742-7843.2009.00485.x
* Wolfinger et al., (1994) Computing gaussian likelihoods and their derivatives for general linear mixed models doi: 10.1137/0915079
* Wright, Stephen, and Jorge Nocedal (2006) "Numerical optimization." Springer

Author: Vladimir Arnautov aka PharmCat
Copyright © 2019 Vladimir Arnautov <mail@pharmcat.net>
