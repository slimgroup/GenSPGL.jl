export spglcore

"""
Use: spglcore()

Main loop of spgl1.jl
"""
function spglcore{Txg<:Number, Titn<:Number}(init::spgInit{Txg,Titn})

    #DEVNOTE# Create spgInit type to hold all these initialized paramaters
    
    println("script has entered spglcore\n 
            Begin Main Loop.")

    # Pull out options type for easier use
    options = init.options
    params = init.params

    #Main Loop
    while true

        # Test Exit Conditions
        # ================================================================================ 
        gNorm::Float64 = zero(Float64)

        if (options.proxy)
            gNorm = options.dual_norm(init.g2, options.weights, init.params)
        else
            gNorm = options.dual_norm(init.g,  options.weights, init.params)
        end
        
        # rNorm and f are the same thing
        rNorm = init.f

        tmp_proj::Array{Txg,1},tmp_itn::Titn = project(init.x - init.g,
                                        init.tau, init.timeProject, options, params)
        Err::Txg = norm(init.x - tmp_proj)
   
        aError1 = rNorm - init.sigma
        aError2 = rNorm^2 - init.sigma^2
        rError1 = abs(aError1) / max(1,rNorm)
        rError2 = abs(aError2) / max(1,init.f)

        # Count number of consecutive iterations with identical support
        nnzOld = init.nnzIdx
    
        nnzX,nnzG,nnzIdx,nnzDiff = activevars(init.x, init.g, init.nnzIdx, options, params)
        
        break #DEVNOTE# Remove this when done main loop
    
    end #Main Loop

end

"""
Use: nnzX,nnzG,nnzIdx,nnzDiff = activevars(x,g,nnzIdx,options, params)

Find the current active set.
nnzX    is the number of nonzero x.
nnzG    is the number of elements in nnzIdx.
nnzIdx  is a vector of primal/dual indicators.
nnzDiff is the no. of elements that changed in the support.
"""
function activevars{ETxg<:Number, Txg<:AbstractVector{ETxg}}(x::Txg,g::Txg,nnzIdx,
                                                    options::spgOptions,params::Dict)

    xTol = min(.1,10*options.optTol)
    gTol = min(.1,10*options.optTol)

    gNorm::ETxg = options.dual_norm(g,options.weights, params)

    nnzOld = nnzIdx

    #Reduced costs for postive and negative parts of x
    z1 = gNorm + g
    z2 = gNorm - g

    #Primal/dual based indicators
    if(~options.proxy)
        xPos = (x .>  xTol) .& (z1 .< gTol)
        xNeg = (x .< -xTol) .& (z2 .< gTol)
        nnzIdx = xPos .| xNeg
    end


end
