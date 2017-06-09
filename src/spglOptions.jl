# GenSPGL
# spglOptions

export spglOptions, spgl_setparms

"""
GenSPGL input options composite type. spgl_setparms is provided as an outer constructor.

"""
type spglOptions
    fid::Integer
    verbosity::Integer
    iterations::Integer
    nPrevVals::Integer
    bpTol::Number
    lsTol::Number
    optTol::Number
    decTol::Number
    stepMin::Number
    stepMax::Number
    rootMethod::Integer
    activeSetIt::Number
    subspaceMin::Bool
    iscomplex::Bool
    maxMatvec::Number
    weights::Number
    quitPareto::Bool
    minPareto::Integer
    lineSrchIt::Integer
    feasSrchIt::Integer
    ignorePErr::Bool
    project::Function
    primal_norm::Function
    dual_norm::Function
    funPenalty::Function 
    proxy::Bool # Appropriate type?
    linear::Bool # ^
    restore::Bool # ^^
end


"""
A placeholder in the default spgl settings until all 
functionality is ported over.

"""
function dummyfun()
    println("This doesn't do anything")
end


"""
#########################################################################################\n
spglOptions Outer Constructor\n
Arguments not specified will be set to their default values. See the spglOptions type docs
for info regarding argument types.\n
#########################################################################################\n

    spgl_setparms(;fid = 1,
                    verbosity = 3,
                    iterations = 100000,
                    nPrevVals = 10,
                    bpTol = 1e-6,
                    lsTol = 1e-6,
                    optTol = 1e-4,
                    decTol = 1e-4,
                    stepMin = 1e-16,
                    stepMax = 1e+5,
                    rootMethod = 2,
                    activeSetIt = Inf,
                    subspaceMin = false,
                    iscomplex = false,
                    maxMatvec = Inf,
                    weights = 1,
                    quitPareto = false,
                    minPareto = 3,
                    lineSrchIt = 1,
                    feasSrchIt = 10000,
                    ignorePErr = false,
                    project = dummyfun,
                    primal_norm = dummyfun,
                    dual_norm = dummyfun,
                    funPenalty = dummyfun,
                    proxy = false,
                    linear = false,
                    restore = false)

# Examples 
```julia
# Default Options
julia> opts = spgl_setparms()

# Make optional changes
julia> opts = spgl_setparms(verbosity = 1, iterations = 20000)
```


Arguments:\n
    |___Field____||__Type__||_____________________Description____________________________|\n
    fid          |Integer|   ..... File ID for output\n
    verbosity    |Integer|   ..... Verbosity level\n
    iterations   |Integer|   ..... Max number of iterations\n
    nPrevVals    |Integer|   ..... Number previous func values for linesearch\n
    bpTol        |Number|    ..... Tolerance for basis pursuit solution \n
    lsTol        |Number|    ..... Least-squares optimality tolerance\n
    optTol       |Number|    ..... Optimality tolerance\n
    decTol       |Number|    ..... Reqd rel. change in primal obj. for Newton\n
    stepMin      |Number|    ..... Minimum spectral step\n
    stepMax      |Number|    ..... Maximum spectral step\n
    rootMethod   |Integer|   ..... Root finding method: 2=quad,1=linear (not used).\n
    activeSetIt, |Number|    ..... Exit with EXIT_ACTIVE_SET if nnz same for # its.\n
    subspaceMin, |Bool|      ..... Use subspace minimization\n
    iscomplex    |Bool|      ..... Flag set to indicate complex problem\n
    maxMatvec    |Number|    ..... Maximum matrix-vector multiplies allowed\n
    weights      |Number|    ..... Weights W in ||Wx||_1\n
    quitPareto   |Bool|      ..... Exits when pareto curve is reached\n
    minPareto    |Integer|   ..... If quitPareto is on, the minimum number of iterations
                                    before checking for quitPareto conditions\n
    lineSrchIt   |Integer|   ..... Maximum number of line search iterations for
                                    spgLineCurvy, originally 10\n
    feasSrchIt   |Integer|   ..... Maximum number of feasible direction line search
                                    iteraitons, originally 10\n
    ignorePErr   |Bool|      ..... Ignores projections error by issuing a warning instead
                                    of an error\n
    project      |Function|  ..... Projection function handle\n
    primal_norm  |Function|  ..... Primal Norm function handle\n
    dual_norm    |Function|  ..... Dual Norm function handle\n
    funPenalty   |Function|  ..... Penalty function handle\n
    proxy        |Bool|      ..... Advanced option that computes pareto curve in a user
                                    specified way. \n
    linear       |Bool|      ..... Advanced option that allows you to declare input
                                    functions to be linear\n
    restore      |Bool|      ..... Whether to restore best previous answer. for large
                                    problems, dont want to do this.\n
"""
function spgl_setparms(;fid = 1,
                        verbosity = 3,
                        iterations = 100000,
                        nPrevVals = 10,
                        bpTol = 1e-6,
                        lsTol = 1e-6,
                        optTol = 1e-4,
                        decTol = 1e-4,
                        stepMin = 1e-16,
                        stepMax = 1e+5,
                        rootMethod = 2,
                        activeSetIt = Inf,
                        subspaceMin = false,
                        iscomplex = false,
                        maxMatvec = Inf,
                        weights = 1,
                        quitPareto = false,
                        minPareto = 3,
                        lineSrchIt = 1,
                        feasSrchIt = 10000,
                        ignorePErr = false,
                        project = dummyfun,
                        primal_norm = dummyfun,
                        dual_norm = dummyfun,
                        funPenalty = dummyfun,
                        proxy = false,
                        linear = false,
                        restore = false)

    opts = spglOptions(fid,
                        verbosity,
                        iterations,
                        nPrevVals,
                        bpTol,
                        lsTol,
                        optTol,
                        decTol,
                        stepMin,
                        stepMax,
                        rootMethod,
                        activeSetIt,
                        subspaceMin,
                        iscomplex,
                        maxMatvec,
                        weights,
                        quitPareto,
                        minPareto,
                        lineSrchIt,
                        feasSrchIt,
                        ignorePErr,
                        project,
                        primal_norm,
                        dual_norm,
                        funPenalty,
                        proxy,
                        linear,
                        restore)



    return opts
end

