- v1.0.0

    * validation
    * test procedures
    * rho & df "problem" solved
    * fix text output
    * twostep deprecated
    * varlink() - optimization with logVar & sigmoid rho
    * rho initial value fix
    * simulation object
    * documuntation

- v0.2.0

    * bugfix
    * optimization
    * documentation
    * validation tests
    * memory caching
    * post optimization (can be done if Newton() failed)

- v0.1.4
    * StatsModels v0.6
    * change struct RBE (some fields renamed, new fixed effect struct)
    * new struct for tables
    * DF calculations
    * Type III Table
    * p values
    * methods to get estimates and contrasts and others
    * optimization, bugfix, cosmetics

- v0.1.3

    * Performance optimization
    * Random DataSet function randrbeds()
    * Increase coverage
    * Add test
    * Code cosmetics

- v0.1.2

    * change optimization algorithm
    * reml2(::RBE, Array{::Float64,1}) function for -REML2 calculation
    * contrast, lsm, emm, lmean utils
    * DF2
    * Optimizations
    * Changes in struct RBE
    * Additional options
    * Code redesign


- v0.1.1
    * change keyword var -> dvar
    * Step 0 variance calculation
    * Split REML β dependent and REML2 β independent
    * Hessian matrix now come from ForwardDiff
    *  g_tol, x_tol, f_tol keywords for Optim
    * Optimization
    * Confidence intervals: confint(::RBE, ::Float64)
    * Show result
    * Bugfix



- v0.1.0
  * Initial alpha version
