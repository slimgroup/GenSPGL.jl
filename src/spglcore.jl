export spglcore

"""
Use: spglcore()

Main loop of spgl1.jl
"""
function spglcore{Txg<:Number}(init::spgInit{Txg})

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
       
        # Declare type of gNorm
        gNorm::Array{Txg,1} = Txg[]

        if (options.proxy)
            gNorm = options.dual_norm(init.g2,options.weights,init.params)
        else
            gNorm = [options.dual_norm(init.g, options.weights)]
        end

        # rNorm and f are the same thing
        rNorm = init.f

        tmp_proj,tmp_itr = project(init.x - init.g, init.tau, init.timeProject, options, params)
        Err = norm(init.x - tmp_proj)
   
        aError1 = rNorm - init.sigma
        aError2 = rNorm^2 - init.sigma^2
        rError1 = abs(aError1) / max(1,rNorm)
        rError2 = abs(aError2) / max(1,init.f)

        # Count number of consecutive iterations with identical support
        nnzOld = init.nnzIdx
        
        break #DEVNOTE# Remove this when done main loop
    end #Main Loop

end
