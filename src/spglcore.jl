export spglcore

"""
Use: spglcore()

Main loop of spgl1.jl
"""
function spglcore(g::AbstractArray{Tg}, g2::Tg, options::spgOptions, params::Dict)

    println("script has entered spglcore\n 
            Begin Main Loop.")

    #Main Loop
    while 1

        # Test Exit Conditions
        # ================================================================================ 
        
        if (options.proxy)
            gNorm = options.dual_norm(g2,options.weights,params)
        else
            gNorm = options.dual_norm(g, options.weights)
        end

    end #Main Loop

end
