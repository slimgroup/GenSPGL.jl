export spglcore

"""
Use: spglcore()

Main loop of spgl1.jl
"""
function spglcore{Tx<:Number, Tg<:Number}(x::AbstractVector{Tx}, tau::Number, sigma::Number,
                            g::AbstractArray{Tg}, g2::Tg, f::Number, timeProject,
                            options::spgOptions, params::Dict)

    #DEVNOTE# Create spgInit type to hold all these initialized paramaters
    
    println("script has entered spglcore\n 
            Begin Main Loop.")

    #Main Loop
    while true

        # Test Exit Conditions
        # ================================================================================ 
       
        # Declare type of gNorm
        gNorm::Array{Tg,1} = Tg[]

        if (options.proxy)
            gNorm = options.dual_norm(g2,options.weights,params)
        else
            gNorm = [options.dual_norm(g, options.weights)]
        end

        # rNorm and f are the same thing
        rNorm = f

        tmp_proj,tmp_itr = project(x - g, tau, timeProject, options, params)
        Err = norm(x - tmp_proj)
   
        aError1 = rNorm - sigma
        aError2 = rNorm^2 - sigma^2
        rError1 = abs(aError1) / max(1,rNorm)
        rError2 = abs(aError2) / max(1,f)

        # Count number of consecutive iterations with identical support
        nnzOld = nnIdx
        
        break #DEVNOTE# Remove this when done main loop
    end #Main Loop

end
