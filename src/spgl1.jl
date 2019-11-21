# SPGL

export spgl1, project, SpotFunForward, snr

"""
# INFO 
    Use:
        x, r, g, info = spgl1(A::AbstractArray, b::AbstractVector{ETb};
                            x::AbstractVector{ETx}=Array{ETb,1}(),
                            tau::AbstractFloat=NaN,
                            sigma::AbstractFloat=NaN,
                            options::spgOptions = spgOptions(),
                            params::Dict{String,Any} = Dict{String,Any}())\n
        Solve regularized composite programs, including:
            a) basis pursuit, basis pursuit denoise, and lasso
            b) non-linear versions of above problems, where forward model is non-linear

# Inputs 
    A:
        - An explicit 'm' x 'n' matrix
        - An implicit JOLI opperator (in development)
        - A non-linear function handle (in development)\n
    b:
        - An m-vector\n
    tau:
        - A non-negative scalar\n
    sigma:
        - If sigma != NaN then GenSPGL will launch into a root-finding mode to find the 
            tau above that solves. In this case it is STRONGLY recommended that tau = 0.\n
    x:
        - An n-vector estimate of the solution (possibly all zeros). If empty, then
            GenSPGL determines the legnth n via n = length(A'b) and sets x0 = zeros(n).\n
    options:
        - An instance of the composite tpe spgOptions. If not specified, default values
            are used. See the spgOptions documents for default values.

# Outputs 
    x:
        - A solution of the problem\n
    r: 
        - The residual, r = b - f(x)\n
    g:
        - The gradient, g = \u2207h(b - f(x))\n
    info:
        - An instance of the spgInfo composite type. See the spgInfo documents for
            more information.
# Author
    Keegan Lensink
        Seismic Laboratory for Imaging and Modeling
        The University of British Columbia
        keeganlensink@gmail.com

    This code is an adaptation of Michael P. Friedlander, Ewout van den Berg, 
    and Aleksandr Aravkin's MATLAB program SPGL1. 

"""
function spgl1(A::TA,
               b::AbstractVector{ETb};
               x::AbstractVector{ETx}=Array{ETb,1}(),
               tau::AbstractFloat=NaN,
               sigma::AbstractFloat=NaN,
               options::spgOptions = spgOptions(),
               params::Dict{String,Any} = Dict{String,Any}()) where
                 {TA<:Union{joAbstractLinearOperator,AbstractArray}, ETx<:Number, ETb<:Number}
    
    REVISION = "0.1"
    DATE = "June, 2017"

    #DEVNOTE# Could make Tau and Sigma nullable types? However, as long as
    # Tau and Sigma are always Float64 this wont be a problem
    (options.verbosity > 1) && println("Script made it to spgl1 for A::AbstractArray")


    # removed for 0.7+ #tic()
    m = length(b)

    #Add options proxy to params dict
    params["proxy"] = options.proxy

    # Check Tau and Sigma
    if isnan(tau) & isnan(sigma)
        tau = 0.
        sigma = 0.
        singleTau = false
    elseif isnan(sigma)
        singleTau = true
    else
        if isnan(tau)
            tau = 0.
        end
        singleTau = false
    end

    # Definitely dont do subspacemin in the non LS case
    #DEVNOTE# Make sure name matches once funLS is written
    string(options.funPenalty)=="GenSPGL.funLS" || (options.subspaceMin = 0) 
    
    # Threshold for signifigant Newton step
    pivTol = 1e-12

    # Account for precision of data format in line search
    if options.stepMin < 10*eps(real(ETx))
        options.stepMin = 10 * eps(real(ETx))
        println("WARNING: options.stepMin is below the precision of the data type. Setting to 10*eps(DT)")
    end

    ##--------------------------------------------------------------------------------
    # Initialize Local Variables
    ##--------------------------------------------------------------------------------  
    iter = 0; itnTotLSQR = 0#Total SPGL1 and LSQR iterations.
    nProdA = 0; nProdAt = 0
    lastFv = [-Inf for i=1:options.nPrevVals] # Last m functions values
    nLineTot = 0            # Total number of linesearch steps
    printTau = false
    nNewton = 0;
    bNorm, b_normalized = options.funPenalty(b, params)
    stat = false
    timeProject = zero(Float64)
    timeMatProd = zero(Float64)
    nnzIter = zero(Int64) # No. of Its with fixed pattern
    nnzIdx = BitArray{1}()  # Active set indicator
    subspace = false        # Flag if did subspace min in current itn
    stepG = 1.0               # Step length for projected gradient
    testUpdateTau = false

    ##-------------------------------------------------------------------------------
    # End Init
    ##-------------------------------------------------------------------------------

    #DEVNOTE# This could be a splitting point for multiple dispatch

    # Determine Initial x, vector length n, and check if complex
    # Explicit Method
    #DEVNOTE# Change name to JOLI  once things are working
    funForward = SpotFunForward
    options.linear = true

    if isempty(x)
        if eltype(A)<:Number
            n = size(A,2)
            realx = isreal(A) & isreal(b)
        else
            x = funForward(A, x, -b) #DEVNOTE# Check to make sure this is tmp
            n = length(x)
            realx = isreal(x) & isreal(b)
        end
        x = zeros(ETb, n)
    else
        n = length(x)
        realx = isreal(x) & isreal(A)
    end

    (eltype(A)<:Number) && (realx = realx & isreal(A))

    # Override if complex flag was used
    isnull(options.iscomplex) || (realx = ~get(options.iscomplex))
    
    # Check if all weights (if any) are strictly positive. In previous
    # versions we also checked if the number of weights was equal to
    # n. In the case of multiple measurement vectors, this no longer
    # needs to apply, so the check was removed.
    any(isinf.(options.weights)) && error("Weights must be finite")
    any(options.weights .> 0) || error("Weights must be strictly positive")
    
    # Quick exit if sigma >= ||b||. Set Tau = 0 to short circuit the loop
    if bNorm <= sigma
        println("W: sigma >= ||b||.  Exact solution is x = 0.")
        tau = zero(real(ETb))
        singleTau = true
    end

    # Don't do subspaceMin if x is complex
    if (~realx & options.subspaceMin)
        println("W: Subspace minimization disabled when variables are complex")
        subspaceMin = false
    end


    #DEVNOTE# Once these are being updated check to see if type should be inferred
            # from promotion rules
    # Pre-Allocate iteration info vectors
    xNorm1 = zeros(Float64, min(options.iterations,10000))
    rNorm2 = zeros(Float64, min(options.iterations,10000))
    lambda = zeros(Float64, min(options.iterations,10000))

    # Create ExitCondition with null trigger
    exit_status = spgExitCondition()

    #Prepare Log Header
    logheader_1 = """
    ================================================================================  
    GenSPGL.jl, Rev $(REVISION), $(DATE)
    ================================================================================\n

    No. Rows                :$(m)               
    No. Columns             :$(n)
    Initial Tau             :$(tau)             
    Penalty                 :$(options.funPenalty)
    Regularizer             :$(options.primal_norm)
    Penalty(b)              :$(bNorm)
    Optimality tol          :$(options.optTol)
    Basis Pursuit tol       :$(options.bpTol)
    Maximum Iterations      :$(options.iterations)
    """
   
    options.verbosity > 1 && println("singleTau: ", singleTau)

    if singleTau 
        (logheader_2 = """
    Target reg. Norm of x   :$(tau)\n

    Iter    Objective      Relative_Error  gNorm         stepG
    ------------------------------------------------------------------------
    """)

    else 
        (logheader_2 = """
    Target Objective        :$(sigma)\n

    Iter   Objective      Relative_Error  RelError  gNorm        stepG  tau
    ----------------------------------------------------------------------------------------
    """)
    end
    
    logheader = logheader_1*logheader_2
    
    # Decide what to do with log header based on verbosity setting
    (options.verbosity > 0) && println(logheader)
    


    # Project the stating point and evaluate function and gradient
    r = typeof(b)() 
   
    
    if isempty(x)#DEVNOTE# matlab Legacy, this can never be invoked. see 90-99
        
        #DEVNOTE# Why copy? Waste of mem and time. Check if b is used again
        r = deepcopy(b) 
        f,g,g2 = options.funCompositeR(A, x, r, funForward, options.funPenalty, timeMatProd, nProdAt)
    else
        x,itn = project(x,tau, timeProject, options, params)
        r = b - funForward(A, x, [], params)
        nProdA += 1
        f,g,g2 = options.funCompositeR(A, x,r, funForward, options.funPenalty, nProdAt, params)
        dx_tmp, itn_tmp = project(x-g, tau, timeProject, options, params)
        dx = dx_tmp - x
        itn += itn_tmp
    end
       
    dxNorm = norm(dx,Inf)
    if dxNorm < (1/options.stepMax)
        gStep = options.stepMax
    else
        gStep = min(options.stepMax, max(options.stepMin, 1/dxNorm))
    end

    # Required for non-monotone strategy
    lastFv[1] = copy(f)
    fBest = copy(f)
    xBest = copy(x)
    fOld = copy(f)

    (options.verbosity > 1) && println("Init Finished")

    # Wrap up initialized variables
    init = spgInit(A,
                    b,
                    x,
                    tau,
                    sigma,
                    g,
                    g2,
                    f,
                    r,
                    nnzIdx,
                    nnzIter,
                    options,
                    params,
                    timeProject,
                    timeMatProd,
                    exit_status,
                    singleTau,
                    bNorm,
                    fOld,
                    testUpdateTau,
                    iter,
                    nNewton,
                    printTau,
                    subspace,
                    stepG,
                    xNorm1,
                    rNorm2,
                    lambda,
                    lastFv,
                    gStep,
                    funForward,
                    nLineTot,
                    nProdAt,
                    nProdA,
                    fBest,
                    xBest)
    
    # Wrap main loop in a function to ease type stability
    init, rNorm, gNorm, rErr  = spglcore(init)

    # Prepare output
    info = spgInfo( init.tau,
                    rNorm,
                    gNorm,
                    rErr,
                    init.exit_status,
                    init.iter,
                    init.nProdA,
                    init.nProdAt,
                    init.nNewton,
                    init.timeProject,
                    init.timeMatProd,
                    init.options,
                    init.xNorm1[1:init.iter],
                    init.rNorm2[1:init.iter],
                    init.lambda[1:init.iter])

    if (options.verbosity > 0) 
        println("-----------------------------------------------------------------------\n")
        print(info.exit_status)
    end
    
    return init.x, init.r, init.g, info
end #func



"""
Use:    x = project(x::AbstractArray, tau::Number, timeProject::Float64
                    options::spgOptions, params::Dict)
"""
function project(x::Tx, tau::Number, timeProject::Float64,
                    options::spgOptions, params::Dict) where
                      {ETx<:Number, Tx<:AbstractVector{ETx}}
    
    (options.verbosity > 1) && println("Begin Project")

    x_out::Tx, itn::Int64 = options.project(x, tau, options.weights, params) 
 
    #DEVNOTE# Replace with @elapsed at call #timeProject[1] += (toc() - tStart)

    (options.verbosity > 1) && println("Finish Project")
    
    return x_out, itn

end



"""
Activated when an explicit or JOLI operator is passed in.
"""
function SpotFunForward(A::TA, x::AbstractArray, g::AbstractArray, params::Dict) where
                       {TA<:Union{joAbstractLinearOperator,AbstractArray}}
    #DEVNOTE# Double check type of g once in use

    isempty(g) && (f = A*x)
    isempty(g) || (f = A'*g)

    return f

end

"""
    Use:
        x, r, g, info = spgl1(A::Function, b::AbstractVector{ETb};
                            x::AbstractVector{ETx}=Array{ETb,1}(),
                            tau::AbstractFloat=NaN,
                            sigma::AbstractFloat=NaN,
                            options::spgOptions = spgOptions(),
                            params::Dict{String,Any} = Dict{String,Any}())\n
"""
function spgl1(A::Function, b::AbstractVector{ETb};
               x::AbstractVector{ETx}=Array{ETb,1}(),
               tau::AbstractFloat=NaN,
               sigma::AbstractFloat=NaN,
               options::spgOptions = spgOptions(),
               params::Dict{String,Any} = Dict{String,Any}()) where
                  {ETx<:Number, ETb<:Number}
    
    REVISION = "0.1"
    DATE = "June, 2017"

    #DEVNOTE# Could make Tau and Sigma nullable types? However, as long as
    # Tau and Sigma are always Float64 this wont be a problem
    (options.verbosity > 1) && println("Script made it to spgl1 for A::Function")


    # removed for 0.7+ #tic()
    m = length(b)

    #Add options proxy to params dict
    params["proxy"] = options.proxy

    # Check Tau and Sigma
    if isnan(tau) & isnan(sigma)
        tau = 0.
        sigma = 0.
        singleTau = false
    elseif isnan(sigma)
        singleTau = true
    else
        if isnan(tau)
            tau = 0.
        end
        singleTau = false
    end

    # Definitely dont do subspacemin in the non LS case
    #DEVNOTE# Make sure name matches once funLS is written
    string(options.funPenalty)=="GenSPGL.funLS" || (options.subspaceMin = 0) 
    
    # Threshold for signifigant Newton step
    pivTol = 1e-12

    # Account for precision of data format in line search
    if options.stepMin < eps(real(ETx))
        options.stepMin = 10 * eps(real(ETx))
        println("WARNING: options.stepMin is below the precision of the data type. Setting to 10*eps(DT)")
    end

    ##--------------------------------------------------------------------------------
    # Initialize Local Variables
    ##--------------------------------------------------------------------------------  
    iter = 0; itnTotLSQR = 0#Total SPGL1 and LSQR iterations.
    nProdA = 0; nProdAt = 0
    lastFv = [-Inf for i=1:options.nPrevVals] # Last m functions values
    nLineTot = 0            # Total number of linesearch steps
    printTau = false
    nNewton = 0;
    bNorm, b_normalized = options.funPenalty(b, params)
    stat = false
    timeProject = zero(Float64)
    timeMatProd = zero(Float64)
    nnzIter = zero(Int64) # No. of Its with fixed pattern
    nnzIdx = BitArray{1}()  # Active set indicator
    subspace = false        # Flag if did subspace min in current itn
    stepG = 1.0               # Step length for projected gradient
    testUpdateTau = false

    ##-------------------------------------------------------------------------------
    # End Init
    ##-------------------------------------------------------------------------------

    #DEVNOTE# This could be a splitting point for multiple dispatch

    # Determine Initial x, vector length n, and check if complex
    # Explicit Method
    funForward = A

    if isempty(x)
        x_tmp = funForward(A, x, -b, params) 
        n = length(x_tmp)
        realx = isreal(x_tmp) & isreal(b)

        x = zeros(n)
    else
        n = length(x)
        realx = isreal(x) & isreal(b)
    end

    (eltype(A)<:Number) && (realx = realx & isreal(A))

    # Override if complex flag was used
    isnull(options.iscomplex) || (realx = ~get(options.iscomplex))
    
    # Check if all weights (if any) are strictly positive. In previous
    # versions we also checked if the number of weights was equal to
    # n. In the case of multiple measurement vectors, this no longer
    # needs to apply, so the check was removed.
    any(isinf.(options.weights)) && error("Weights must be finite")
    any(options.weights .> 0) || error("Weights must be strictly positive")
    
    # Quick exit if sigma >= ||b||. Set Tau = 0 to short circuit the loop
    if bNorm <= sigma
        println("W: sigma >= ||b||.  Exact solution is x = 0.")
        tau = 0
        singleTau = true
    end

    # Don't do subspaceMin if x is complex
    if (~realx & options.subspaceMin)
        println("W: Subspace minimization disabled when variables are complex")
        subspaceMin = false
    end


    #DEVNOTE# Once these are being updated check to see if type should be inferred
            # from promotion rules
    # Pre-Allocate iteration info vectors
    xNorm1 = zeros(Float64, min(options.iterations,10000))
    rNorm2 = zeros(Float64, min(options.iterations,10000))
    lambda = zeros(Float64, min(options.iterations,10000))

    # Create ExitCondition with null trigger
    exit_status = spgExitCondition()

    #Prepare Log Header
    logheader_1 = """
    ================================================================================  
    GenSPGL.jl, Rev $(REVISION), $(DATE)
    ================================================================================\n

    No. Rows                :$(m)               
    No. Columns             :$(n)
    Initial Tau             :$(tau)             
    Penalty                 :$(options.funPenalty)
    Regularizer             :$(options.primal_norm)
    Penalty(b)              :$(bNorm)
    Optimality tol          :$(options.optTol)
    Basis Pursuit tol       :$(options.bpTol)
    Maximum Iterations      :$(options.iterations)
    """
   
    options.verbosity > 1 && println("singleTau: ", singleTau)

    if singleTau 
        (logheader_2 = """
    Target reg. Norm of x   :$(tau)\n

    Iter    Objective      Relative_Error  gNorm         stepG 
    -------------------------------------------------------------------
    """)

    else 
        (logheader_2 = """
    Target Objective        :$(sigma)\n

    Iter   Objective      Relative_Error  RelError  gNorm         stepG  tau
    -------------------------------------------------------------------------
    """)
    end
    
    logheader = logheader_1*logheader_2
    
    # Decide what to do with log header based on verbosity setting
    (options.verbosity > 0) && println(logheader)
    


    # Project the stating point and evaluate function and gradient
    r = typeof(b)() 
   
    
    if isempty(x)#DEVNOTE# matlab Legacy, this can never be invoked. see 90-99
        
        #DEVNOTE# Why copy? Waste of mem and time. Check if b is used again
        r = deepcopy(b) 
        f,g,g2 = options.funCompositeR(A, x, r, funForward, options.funPenalty, timeMatProd, nProdAt)
    else
        x,itn = project(x,tau, timeProject, options, params)
        r = b - funForward(A, x, Array{ETx,1}(), params)[1]
        nProdA += 1
        f,g,g2 = options.funCompositeR(A, x,r, funForward, options.funPenalty, nProdAt, params)
        dx_tmp, itn_tmp = project(x-g, tau, timeProject, options, params)
        dx = dx_tmp - x
        itn += itn_tmp
    end

    dxNorm = norm(dx,Inf)
    if dxNorm < (1/options.stepMax)
        gStep = options.stepMax
    else
        gStep = min(options.stepMax, max(options.stepMin, 1/dxNorm))
    end

    # Required for non-monotone strategy
    lastFv[1] = copy(f)
    fBest = copy(f)
    xBest = copy(x)
    fOld = copy(f)

    (options.verbosity > 1) && println("Init Finished")
    
    # Wrap up initialized variables
    init = spgInit(A,
                    b,
                    x,
                    tau,
                    sigma,
                    g,
                    g2,
                    f,
                    r,
                    nnzIdx,
                    nnzIter,
                    options,
                    params,
                    timeProject,
                    timeMatProd,
                    exit_status,
                    singleTau,
                    bNorm,
                    fOld,
                    testUpdateTau,
                    iter,
                    nNewton,
                    printTau,
                    subspace,
                    stepG,
                    xNorm1,
                    rNorm2,
                    lambda,
                    lastFv,
                    gStep,
                    funForward,
                    nLineTot,
                    nProdAt,
                    nProdA,
                    fBest,
                    xBest)
    
    # Wrap main loop in a function to ease type stability
    init, rNorm, gNorm, rErr  = spglcore(init)
    # Prepare output
    info = spgInfo(  init.tau,
                    rNorm,
                    gNorm,
                    rErr,
                    init.exit_status,
                    init.iter,
                    init.nProdA,
                    init.nProdAt,
                    init.nNewton,
                    init.timeProject,
                    init.timeMatProd,
                    init.options,
                    init.xNorm1[1:init.iter],
                    init.rNorm2[1:init.iter],
                    init.lambda[1:init.iter])

    if (options.verbosity > 0) 
        println("---------------------------------------------------------------------------\n")
        print(info.exit_status)
    end
    
    return init.x, init.r, init.g, info
end #func


snr(raw,interp) = -20log10(norm(interp-raw)/norm(raw))
