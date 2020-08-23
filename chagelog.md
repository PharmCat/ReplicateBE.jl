- v1.0.11
    * Optimizations
    * Cosmetics
    * PDMats exclude
    * Fix experimental sigmoid function

- v1.0.10
    * Cosmetics
    * StatsBase.coeftable
    * StatsBase.vcov
    * Minor optimizations
    * Some stable tests for 1.5

- v1.0.9
    * PDMats 0.10
    * Opimizations

- v1.0.8
    * Changing in RBE struct!
    * Optim 0.21
    * StatsBase 0.33
    * All output now in result: theta, reml, G, H, C, A, optim results
    * fix test with accordance with public API
    * change RBEResult

- v1.0.7
    * settings: initial
    * settings: initial
    * randrbeds "light" generation
    * optimizations
    * Distributions bump
    * "generalized" simulation
    * Optim 0.20
    * Julia 1.3

- v1.0.6
    * many optimizations
    * linear algebra custom functions
    * drop all interim data
    * improve code structure

- v1.0.5
    * inv(H) used where possible
    * DF for typeiii checks
    * Type optimization
    * Output fix

- v1.0.4
    * CSV test dataset
    * rho linking function fix
    * docs update

- v1.0.3
    * output fix
    * add "linking" functions
    * cosmetics

- v1.0.2 (hotfix 2)
    * cmat fix
    * cosmetics

- v1.0.1 (hotfix)
    * increase stability
    * simulation stability
    * change default options: g_tol = 1e-12
    * vlm = 1.0 in rbe() & rbe!()
    * initvar() changes
    * theta() - return final estimates

- v1.0.0 (unstable)
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
