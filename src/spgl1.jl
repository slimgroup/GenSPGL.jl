# SPGL

export spgl1, project

"""
This will contain info on use of spgl1

EXPLICIT METHOD

When implementing JOLI support, provide new method. e.g A::joOp....
"""
function spgl1{Tx<:AbstractFloat, Tb<:Number}(A::AbstractArray, b::AbstractVector{Tb};
                    x::AbstractVector{Tx}=Array{Tb,1}(),
                    tau::AbstractFloat=NaN,
                    sigma::AbstractFloat=NaN,
                    options::spgOptions = spgOptions(),
                    params::Dict{String,Number} = Dict{String,Number}())
    
    REVISION = "0.1"
    DATE = "June, 2017"

    #DEVNOTE# Could make Tau and Sigma nullable types? However, as long as
    # Tau and Sigma are always Float64 this wont be a problem
    println("Script made it to spgl1")

    tic()
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



    ##--------------------------------------------------------------------------------
    # Initialize Local Variables
    ##--------------------------------------------------------------------------------  
    iter = 0; itnTotLSQR = 0#Total SPGL1 and LSQR iterations.
    nProdA = [0]; nProdAt = [0]
    lastFv = [-Inf for i=1:options.nPrevVals] # Last m functions values
    nLineTot = 0            # Total number of linesearch steps
    pintTau = false
    nNewton = 0;
    bNorm, b_normalized = options.funPenalty(b, params)
    stat = false
    timeProject = Float64[0]
    timeMatProd = Float64[0]
    nnzIter = 0             # No. of Its with fixed pattern
    nnzIdx = []             # Active set indicator
    subspace = false        # Flag if did subspace min in current itn
    stepG = 1               # Step length for projected gradient
    testUpdateTau = 0

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
        x = zeros(n)
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
    any(isinf(options.weights)) && error("Weights must be finite")
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
    
    singleTau && (logheader_2 = """
    Target reg. Norm of x   :$(tau)\n

    Iter    Objective   Relative_Error  gNorm   stepG   nnzX    nnzG
    --------------------------------------------------------------------------------\n
    """)

    singleTau || (logheader_2 = """
    Target Objective        :$(sigma)\n

    Iter    Objective   Relative_Error  Rel_Error   gNorm   stepG   nnzX    nnzG    tau
    -----------------------------------------------------------------------------------\n
    """)
    logheader = logheader_1*logheader_2
    
    # Decide what to do with log header based on verbosity setting
    (options.verbosity == 1) && println(logheader)
    


    # Project the stating point and evaluate function and gradient
    r = typeof(b)() 
   
    
    if isempty(x)#DEVNOTE# matlab Legacy, this can never be invoked. see 90-99
        
        #DEVNOTE# Why copy? Waste of mem and time. Check if b is used again
        r = deepcopy(b) 
        f,g,g2 = funCompositeR(A, x, r, funForward, options.funPenalty, timeMatProd, nProdAt)
    else
        x,itn = project(x,tau, timeProject, options, params)
        r = b - funForward(A, x, [], params)
        nProdA += 1
        f,g,g2 = funCompositeR(A, x,r, funForward, options.funPenalty, nProdAt, params)
        dx_tmp, itn_tmp = project(x-g, tau, timeProject, options, params)
        dx = dx_tmp - x
        itn += itn_tmp
    end
       
    dxNorm = norm(dx,Inf)
    if dxNorm < (1/options.stepMax)
        gStep = stepMax
    else
        gStep = min(options.stepMax, max(options.stepMin, 1/dxNorm))
    end

    # Required for non-monotone strategy
    lastFv[1] = f
    fBest = f
    xBest = x
    fOld = f

    println("Init Finished")

    # Wrap up initialized variables
    init = spgInit(x,
                    tau,
                    sigma,
                    g,
                    g2,
                    f,
                    nnzIdx,
                    options,
                    params,
                    timeProject)
                    
    # Wrap main loop in a function to ease type stability
    @code_warntype spglcore(init)

    return init
end #func



"""
Use:    x = project(x::AbstractArray, tau::Number, timeProject::AbstractArray,
                    options::spgOptions, params::Dict)
"""
function project{Tx<:Number}(x::AbstractArray{Tx}, tau::Number, timeProject::AbstractArray,
                    options::spgOptions, params::Dict)
    
    (options.verbosity == 1) && println("Begin Project")

    x_out::typeof(x), itn::Int = options.project(x, tau, options.weights, params) 
 
    #DEVNOTE# Replace with @elapsed at call #timeProject[1] += (toc() - tStart)

    (options.verbosity == 1) && println("Finish Project")
    
    return x_out,itn

end



"""
GenSPGL

Use:    f,g1,g2 = funCompositeR(A, r, funForward, funPenalty, params, nProdAt, params)
"""
function funCompositeR(A::AbstractArray,x::AbstractArray,r::AbstractArray,
                        funForward::Function, funPenalty::Function, 
                        nProdAt::AbstractArray,
                        params::Dict{String,Number})

    nProdAt[1] += one(eltype(nProdAt))
    f,v = funPenalty(r, params)
    
    if ~(params["proxy"])
        g1 = funForward(A, x, -v, params)
        g2 = zero(eltype(g1))
    else
        g1,g2 = funForward(A, x, -v, params)
    end

    return f,g1,g2
end



"""
Activated when an explicit or JOLI operator is passed in.
"""
function SpotFunForward(A::AbstractArray, x::AbstractArray, g::AbstractArray, params::Dict)
    #DEVNOTE# Double check type of g once in use

    isempty(g) && (f = A*x)
    isempty(g) || (f = A'*g)

    return f

end
