# SPGL

export spgl1

"""
This will contain info on use of spgl1

EXPLICIT METHOD

When implementing JOLI support, provide new method. e.g A::joOp....
"""
function spgl1{xT<:AbstractFloat}(A::AbstractArray, b::AbstractArray;
                    x::AbstractArray{xT}=Float64[],
                    tau::AbstractFloat=NaN,
                    sigma::AbstractFloat=NaN,
                    options::spgOptions = spgOptions(),
                    params::Dict{String,Number} = Dict{String,Number}())

    #DEVNOTE# Could make Tau and Sigma nullable types? However, as long as
    # Tau and Sigma are always Float64 this wont be a problem
    println("Script made it to spgl1")

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
    nProdA = 0; nProdAt = 0
    lastFv = [-Inf for i=1:options.nPrevVals] # Last m functions values
    nLineTot = 0            # Total number of linesearch steps
    pintTau = false
    nNewton = 0;
    
    #DEVNOTE# Consider making a composite params type for each funPenalty instead of
            # the current Dict solution
    bNorm, b_normalized = options.funPenalty(b, params)

    stat = false
    timeProject = 0
    timeMatProd = 0
    nnzIter = 0             # No. of Its with fixed pattern
    nnzIdx = []             # Active set indicator
    subspace = false        # Flag if did subspace min in current itn
    stepG = 1               # Step length for projected gradient
    testUpdateTau = 0

    ##-------------------------------------------------------------------------------
    # End Init
    ##-------------------------------------------------------------------------------



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
            x = funForward(x, -b) #DEVNOTE# Check to make sure this is tmp
            n = length(x)
            realx = isreal(x) & isreal(b)
        end
        x = zeros(n,1)
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
    isinf(options.weights) && error("Weights must be finite")
    options.weights > 0 || error("Weights must be strictly positive")
    
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

    # Exit Condition Dict

    
end #func


"""
Activated when an explicit or JOLI operator is passed in.
"""
function SpotFunForward(x::AbstractArray, g::AbstractArray, params::Dict)
    #DEVNOTE# Double check type of g once in use

    isempty(g) && (f = A*x)
    isempty(g) || (f = A'*g)

    return f

end
