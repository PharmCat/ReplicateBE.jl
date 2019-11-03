## Details

[![doi](https://img.shields.io/badge/doi-10.13140%2FRG.2.2.27418.39363-blue)](https://doi.org/10.13140/RG.2.2.27418.39363)

### Introduction

The replicate designed bioequivalence is a powerful approach to get more information about variation. In some cases the number of subjects required to demonstrate bioequivalence can be reduced by up to about 50% (Van Peer, A., 2010). For a high variability product, replication can really improve the precision and provide more complete intra-individual variation estimate. Also replicate design could be used for reference-scaled average bioequivalence (RSABE) to demonstrate bioequivalence for highly variable drugs (HVDs).
With accordance to US FDA guideline linear mixed-effects model procedures, available in PROC MIXED in SAS or equivalent software, should be used for the analysis of replicated crossover studies for average BE (US FDA).
At this moment linear mixed model effect analysis can be done with proprietary (SPSS, SAS, Stata) and open source (R:nlme, R:lme4, Julia:MixedModels) software. But not all statistical mixed models packages support flexible covariance structure fitting with structures like “heterogeneous compound symmetry” (CSH), FA0(2). This doesn’t means that lme4 or MixedModels can’t be used for bioequivalence estimation,  but CSH structure not available in this packages and comparison of results performed in SAS/SPSS with lme4 can be problematically.
Objective of this work is: to provide instrument to make bioequivalence analysis with type C model and a development of a demonstrative code for step-by-step clarification of mixed model computation procedure for any interested developers.

### Materials and Methods

FDA recommended model can be described with following equation (Patterson, 2002; US FDA):

```math
 Y_{ijkl} = \mu_k + \gamma_{ikl} + \delta_{ijk} + \varepsilon_{ijkl}
```

Where  ``i=1,…s`` indicates sequence, ``j=1,…n_i`` - subjects,  ``k=R,T`` – treatment, ``l=1,2`` indicates replicate on treatment ``k`` for subjects within sequence ``i``. ``Y_ijkl`` is the response of replicate ``l`` on treatment ``k`` for subject ``j`` in sequence ``i``, ``γ_ikl`` represents the fixed effect of replicate ``l`` on treatment ``k`` in sequence ``i``, ``δ_ijk`` is the random subject effect for subject ``j``  in sequence ``i`` on treatment ``k``, and ``ε_{ijkl}`` is the random error for subject ``j`` within sequence ``i`` on replicate ``l`` of treatment ``k``.  The ``ε_{ijkl}`` are assumed to be mutually independent and identically distributed as

```math
\varepsilon_{ijkl} \sim N(0, \sigma_{Wk}^2)
```

And the random subject effects are assumed to be mutually independent and distributed as

```math
\delta_{ij} \sim N_2
\begin{bmatrix}
    \begin{pmatrix}  \mu_R  \\ \mu_T \end
    {pmatrix}\begin{pmatrix}\sigma_{BR}^{2} &  \rho\sigma_{BT}\sigma_{BR} \\
    \rho\sigma_{BT}\sigma_{BR} & \sigma_{BR}^{2}
    \end{pmatrix}
\end{bmatrix}
```

Following code illustrates an example of program statements to run the average bioequivalence analysis using PROC MIXED in SAS:

```sas
PROC MIXED;
CLASSES SEQ SUBJ PER TRT;
MODEL  Y = SEQ PER TRT/ DDFM=SATTERTH;
RANDOM  TRT/TYPE=FA0(2) SUB=SUBJ G;
REPEATED/GRP=TRT SUB=SUBJ;
ESTIMATE 'T vs. R' TRT 1 -1/CL ALPHA=0.1;
```
Statement TYPE=CSH also can be used to match the model described above.

In matrix notation a mixed effect model can be represented as:

```math
y = X\beta + Zu + \epsilon
```

And gives Henderson's «mixed model equations»:

```math
\begin{pmatrix}X'R^{-1}X&X'R^{-1}Z\\Z'R^{-1}X&Z'R^{-1}Z+G_{-1}\end{pmatrix}  \begin{pmatrix}\widehat{\beta} \\ \widehat{u} \end{pmatrix}= \begin{pmatrix}X'R^{-1}y\\Z'R^{-1}y\end{pmatrix}
```

The solution to the mixed model equations is a maximum likelihood estimate when the distribution of the errors is normal. PROC MIXED in SAS used restricted maximum likelihood (REML) approach by default. REML equation can be described with following (Henderson,  1959;Laird et.al. 1982; Jennrich 1986; Lindstrom & Bates, 1988; Gurka et.al 2006):

```math
logREML(\theta,\beta) = -\frac{N-p}{2} - \frac{1}{2}\sum_{i=1}^nlog|V_{i}|-

-\frac{1}{2}log|\sum_{i=1}^nX_i'V_i^{-1}X_i|-\frac{1}{2}\sum_{i=1}^n(y_i - X_{i}\beta)'V_i^{-1}(y_i - X_{i}\beta)
```

Where

```math
\beta = {(\sum_{i=1}^n X_{i}'V_i^{-1}X_{i})}^{-1}(\sum_{i=1}^n X_{i}'V_i^{-1}y_{i})
```

Where

```math
V_{i} = Z_{i}GZ_i'+R_{i}
```

Where ``N`` – total number of observations, ``n`` – number of independent sampling units (subjects), ``y_i`` individual response vector,  ``X_i`` individual design matrix of fixed effects, ``β`` vector of fixed effects parameters, ``V_i`` individual covariance matrix for the response vector, ``Z_i`` individual design matrix of random effects,  ``G`` covariance matrix of ``u`` (random effect), ``R_i`` individual covariance matrix  of  ``ϵ`` (residual error).
Finding solution for minimization  -2logL(θ) respectively to θ can be done with Newton’s family methods. In ReplicateBE used optimization with Optim.jl package (Newton's Method). In some cases post-optimization step can be performed with  Broyden–Fletcher–Goldfarb–Shanno  method ((L)-BFGS)(Fletcher & Roger, 1987; Wright, 2006). Because variance have only positive values and ρ is limited as -1 ≤ ρ ≤1 in CSH (SAS implementation) and 0 ≤ ρ ≤1 in ReplicateBE "link" function is used. Exponential values is optimizing in variance part and ρ is linked with sigmoid function.
All steps perform with differentiable functions with forward automatic differentiation using ForwardDiff package. ForwardDiff is a Julia package for forward-mode automatic differentiation (AD) featuring performance competitive with low-level languages like C++. Unlike recently developed AD tools in other popular high-level languages such as Python and MATLAB, ForwardDiff takes advantage of just-in-time (JIT) compilation to transparently recompile AD-unaware user code, enabling efficient support for higher-order differentiation and differentiation using custom number types (including complex numbers). The field of automatic differentiation provides methods for automatically computing exact derivatives (up to floating-point error) given only the function itself (Revels et al., 2016; Mogensen et al., 2018).

After solving optimization problem other statistical parameters can be found (Giesbrecht & Burns, 1985; Hrong-Tai Fai & Corneliu 1996; Schaalje et al 2002):

```math
F = \frac{\beta'L'(LCL')^{-1}L\beta}{rank(LCL')}
```

Where

```math
C = \sum_{i=1}^{n} X_i'V_i^{-1}X_i
```

And

```math
se = \sqrt{LCL'}
```

Degree of freedom (DF) computed with Satterthwaite approximation or with “contain” method (N – rank(XZ)).

```math
df = \frac{2(LCL')^{2}}{g'Ag}
```

Where  ``A = 2H^{-1}``; ``g = \triangledown_{\theta}(LC^{-1}_{\theta}L')``


Where L is a vector of known constant, C – variance-covariance matrix of fixed effects (var(β)), H – hessian matrix of REML function, N – total number of observations.

### Validation

ReplicateBE was validated with 6 reference public datasets, 24 generated datasets and simulation study. ReplicateBE version 0.1.4 and 0.2.0 is compliant to SAS/SPSS, values checked: REML estimate, variance components estimate, fixed effect estimate, standard error of fixed effect estimate. Validation procedures included in package test procedure and perform each time when new version released or can be done at any time on user machine. Confidence interval (95%) for type I error (alpha) is 0.048047 - 0.050733 (10000 iterations). No statistically significant difference found with acceptable rate (0.05) found (version 0.1.4). Testing procedures cover approximately [![codecov](https://codecov.io/gh/PharmCat/ReplicateBE.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/PharmCat/ReplicateBE.jl) of code, and perform for each release on Travis CI platform: [![Build Status](https://api.travis-ci.com/PharmCat/ReplicateBE.jl.svg?branch=master)](https://travis-ci.com/PharmCat/ReplicateBE.jl).

### Installation and using

Installation:

```julia
using Pkg; Pkg.add("ReplicateBE")
```

Basic using:

```julia
using ReplicateBE, StatsBase
be = ReplicateBE.rbe!(df, dvar = :var, subject = :subject, formulation = :formulation, period = :period, sequence = :sequence);
ci = confint(be, 0.1)
```

Standard output:

```
Bioequivalence Linear Mixed Effect Model (status: converged)

-2REML: 329.257    REML: -164.629

Fixed effect:
───────────────────────────────────────────────────────────────────────────────────────────
Effect           Value         SE          F          DF        t           P|t|
───────────────────────────────────────────────────────────────────────────────────────────
(Intercept)      4.42158       0.119232    1375.21    68.9135   37.0838     2.90956e-47*
sequence: 2      0.360591      0.161776    4.96821    63.0139   2.22895     0.0293917*
period: 2        0.027051      0.0533388   0.257206   159.645   0.507155    0.612746
period: 3        -0.00625777   0.0561037   0.012441   199.483   -0.111539   0.911301
period: 4        0.036742      0.0561037   0.428886   199.483   0.654894    0.51329
formulation: 2   0.0643404     0.0415345   2.39966    207.651   1.54908     0.122884
───────────────────────────────────────────────────────────────────────────────────────────
Intra-individual variation:
formulation: 1   0.108629
formulation: 2   0.0783544
Inter-individual variation:
formulation: 1   0.377846
formulation: 2   0.421356
ρ:               0.980288   Cov: 0.391143
Confidence intervals(90%):
formulation: 1 / formulation: 2
99.5725 - 114.221 (%)
formulation: 2 / formulation: 1
87.5496 - 100.4293 (%)
```

### Results

ReplicateBE was developed to get mixed model solution to bioequivalence clinical trial. Package repository: [![GitHub version](https://badge.fury.io/gh/PharmCat%2FReplicateBE.jl.svg)](https://badge.fury.io/gh/PharmCat%2FReplicateBE.jl), Julia 1.0.5 or latest should be installed.

### Discussion

ReplicateBE not designed for modeling in a general purpose, but can be used in situation with similar structure. In part of datasets ReplicateBE showed better optimization result as SPSS. Also ReplicateBE based on direct inversing of variance-covarance matrix V, so computation of ``V^(-1)`` may be time expensive if size of matrix is big. This does not happen in bioequivalence study where size of ``V`` is no more 4 (4 periods). But in general this can be serious disadvantage. This situation can be avoided using sweep based transformations (Wolfinger et al., 1994). In ReplicateBE variance structure strictly denoted and can’t be changed, but it can be a target in package developing path. In ReplicateBE Satterthwaite degree of freedom (DF) not equal with SAS/SPSS DF estimate in part of datasets.

### Acknowledgments

D.Sc. in Physical and Mathematical Sciences Anastasia Shitova a.shitova@qayar.ru

### Literature Cited

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
