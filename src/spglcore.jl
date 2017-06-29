export spglcore

"""
Use: spglcore()

Main loop of spgl1.jl
"""
function spglcore{Txg<:Number, Tidx<:BitArray}(init::spgInit{Txg, Tidx})

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
        gNorm::Txg = zero(Txg)

        if (options.proxy)
            gNorm = options.dual_norm(init.g2, options.weights, init.params)
        else
            gNorm = options.dual_norm(init.g,  options.weights, init.params)
        end
        
        # rNorm and f are the same thing
        rNorm = init.f

        tmp_proj::Array{Txg,1},tmp_itn::Int64 = project(init.x - init.g,
                                        init.tau, init.timeProject, options, params)
        Err::Txg = norm(init.x - tmp_proj)
        rErr::Txg = Err/max(1,init.f)
   
        aError1 = rNorm - init.sigma
        aError2 = rNorm^2 - init.sigma^2
        rError1 = abs(aError1) / max(1,rNorm)
        rError2 = abs(aError2) / max(1,init.f)

        # Count number of consecutive iterations with identical support
        
        nnzOld::Tidx = deepcopy(init.nnzIdx) #DEVNOTE# Not stable, but shouldnt be a big deal
        nnzX,nnzG,init.nnzIdx,nnzDiff = activevars(init.x, init.g, init.nnzIdx, options, params)
      
        println("""
        nnzX:       $nnzX\n
        nnzG:       $nnzG\n
        nnzIdx:     $(init.nnzIdx)\n
        nnzDiff:    $nnzDiff\n
        """)
        
        if (nnzDiff == -1)
            init.nnzIter = 0
        end

        options.verbosity == 1 && println("fin CompConditions")
        
        if isempty(init.x)
            println("DEVNOTE x is empty in spglcore, check spgl1 90-99")
        else
            init.nnzIter += 1
            (init.nnzIter >= options.activeSetIt) && (init.exit_status.triggered = 9)
            nnzX = sum(abs.(init.x) .>= min(0.1,10*options.optTol))::Int64
        end

        # SingleTau, check if optimal
        # Second condition guards against large tau
        if init.singleTau
            if ((rErr .<= options.optTol) | (rNorm .< options.optTol*init.bNorm))
                init.exit_status.triggered = 4 
            end
        else
            # Test for LS solution found
            if gNorm <= options.lsTol
                init.exit_status.triggered = 3
            end

            if (rErr <= max(options.optTol, rError2)) | rError1 <= optTol

                # Problem is nearly optimal for current tau
                # Check optimailty of current root
                test1 = (rNorm <= options.bpTol*init.bNorm)::Bool
                test3 = (rError1 <= options.optTol)::Bool
                test4 = (rNorm <= init.sigma)::Bool
                
                test1 && (init.exit_status.triggered = 7)
                test3 && (init.exit_status.triggered = 1)
                test4 && (init.exit_status.triggered = 2)

            end

            testRelChange1 = (abs(init.f - init.fOld) <= options.decTol * init.f)::Bool
            testRelChange2 = (abs(init.f - init.fOld) <= 1e-1*init.f*(abs(rNorm - init.sigma)))
            init.testUpdateTau = ((testRelChange1 & rNorm >  2*init.sigma) |
                            (testRelChange2 & rNorm <= 2*init.sigma)) &
                            isnull(init.exit_status.triggered) &
                            ~init.testUpdateTau

            if testUpdateTau

                if (options.quitPareto & iter >= minPareto) 
                    init.exit_status.triggered = 10
                end

                tauOld = copy(init.tau)
                init.tau = max(zero(typeof(init.tau)), tau + (aError1)/ gNorm)::typeof(tauOld)
                init.nNewton += one(typeof(init.nNewton))
                
                #DEVNOTE# Unstable, but may remove later
                printTau = (abs(tauOld - init.tau) >= 1e-6 * init.tau)::Bool 

                if tau < tauOld
                    init.x, tmp_itn = project(init.x, init.tau, init.timeProject, options, params)
                end

            end

        end

        if (isnull(init.exit_status.triggered) & init.iter >= options.iterations)
            init.exit_status.triggered = 5 
        end

        (options.verbosity == 1) && println("fin CheckConverge")


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
function activevars{Ti<:BitArray{1}, ETxg<:Number, Txg<:AbstractVector{ETxg}}(x::Txg,g::Txg,
                        nnzIdx::Ti, options::spgOptions,params::Dict{String,Number})

    xTol::ETxg = min(.1,10*options.optTol)
    gTol::ETxg = min(.1,10*options.optTol)

    gNorm::ETxg = options.dual_norm(g,options.weights, params)
    
    if isnull(nnzIdx)
        nnzOld_inner = Ti()
    else
        nnzOld_inner = copy(nnzIdx)
    end

    #Reduced costs for postive and negative parts of x
    z1 = gNorm + g
    z2 = gNorm - g

    #Primal/dual based indicators
    xPos = BitArray{1}()
    xNeg = BitArray{1}()

    if(~options.proxy)
        xPos = (x .>  xTol) .& (z1 .< gTol)
        xNeg = (x .< -xTol) .& (z2 .< gTol)
        nnzIdx = xPos .| xNeg
    end

    nnzX = sum(abs.(x) .>= xTol)::Int64
    nnzG = sum(nnzIdx)::Int64

    if isempty(nnzOld_inner)
        #DEVNOTE# Use -1 instead of Inf? Inf creates a type stability
                # Need to document what this means though
        nnzDiff = -one(Int64)
    else
        nnzDiff = sum(nnzIdx .!== nnzOld_inner)::Int64
    end

    return nnzX, nnzG, nnzIdx, nnzDiff

end
