export spglcore

"""
Use: spglcore()

Main loop of spgl1.jl
"""
function spglcore{ETxg<:Number, Txg<:AbstractVector{ETxg}, Tidx<:BitArray}(init::spgInit{ETxg, Txg, Tidx})

    #DEVNOTE# Create spgInit type to hold all these initialized paramaters
    
    (init.options.verbosity > 1) && println("script has entered spglcore\n 
            Begin Main Loop.")

    # Pull out options type for easier use
    options = init.options
    params = init.params

    #Main Loop
    while true

        # ================================================================================
        # Test Exit Conditions
        # ================================================================================ 
        gNorm::ETxg = zero(ETxg)

        if (options.proxy)
            gNorm = options.dual_norm(init.g2, options.weights, init.params)
        else
            gNorm = options.dual_norm(init.g,  options.weights, init.params)
        end
        
        # rNorm and f are the same thing
        rNorm = init.f

        tmp_proj::Txg,tmp_itn::Int64 = project(init.x - init.g,
                                        init.tau, init.timeProject, options, params)
        Err::ETxg = norm(init.x - tmp_proj)
        rErr::ETxg = Err/max(1,init.f)
   
        aError1 = rNorm - init.sigma
        aError2 = rNorm^2 - init.sigma^2


        # Count number of consecutive iterations with identical support
        
        nnzOld::Tidx = deepcopy(init.nnzIdx) #DEVNOTE# Not stable, but shouldnt be a big deal
        nnzX,nnzG,init.nnzIdx,nnzDiff = activevars(init.x, init.g, init.nnzIdx, options, params)
      
        if (nnzDiff == -1)
            init.nnzIter = 0
        end

        options.verbosity > 1 && println("fin CompConditions")
        
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

                if (options.quitPareto & init.iter >= minPareto) 
                    init.exit_status.triggered = 10
                end

                tauOld = copy(init.tau)
                init.tau = max(zero(typeof(init.tau)), tau + (aError1)/ gNorm)::typeof(tauOld)
                init.nNewton += one(typeof(init.nNewton))
                
                #DEVNOTE# Unstable, but may remove later
                init.printTau = (abs(tauOld - init.tau) >= 1e-6 * init.tau)::Bool 

                if tau < tauOld
                    init.x, tmp_itn = project(init.x, init.tau, init.timeProject, options, params)
                end

            end

        end

        if (isnull(init.exit_status.triggered) & init.iter >= options.iterations)
            init.exit_status.triggered = 5 
        end

        (options.verbosity > 1) && println("fin CheckConverge")

        # ===============================================================================
        # Print log, update history, and check exit conditions
        # ===============================================================================

        if (options.verbosity > 0) | init.singleTau | init.printTau |
                                 (init.iter == 0) | ~isnull(init.exit_status.triggered)

            if init.singleTau
                
                s = @sprintf "%5i   %13.7e  %13.7e  %9.2e   %6.1f   %6i     %6i" init.iter rNorm rErr rNorm log10(init.stepG) nnzX nnzG
                (options.verbosity > 0) && println(s)

                if init.subspace
                    println("$itnLSQR")
                end
            else
                
                #DEVNOTE# Check ML ver. This line should be different
                s = @sprintf "%5i   %13.7e  %13.7e  %9.2e   %6.1f   %6i     %6i" init.iter rNorm rErr rNorm log10(init.stepG) nnzX nnzG
                (options.verbosity > 0) && println(s)

                if init.printTau | init.subspace
                    println("$Tau   $itnLSQR\n")
                end
            end

        end

        init.printTau = false
        init.subspace = false

        #Update History
        if isempty(init.x)
            init.xNorm1[init.iter+1] = 0.
        else
            init.xNorm1[init.iter+1] = options.primal_norm(init.x, options.weights, params)
        end
        init.rNorm2[init.iter+1] = copy(rNorm)
        init.lambda[init.iter+1] = copy(gNorm)
        
        # ================================================================================
        # Begin Iterations
        # ================================================================================
        init.iter += 1
        xOld = copy(init.x)
        fOld = copy(init.f)
        rOld = copy(init.r)
       
        try
            # ================================================================================
            # Projected gradient step and line search
            # ================================================================================

            (options.verbosity > 1) && println("begin LineSearch")

            init.f, init.x, init.r, nLine, stepG, lnErr, localProdA = spglinecurvy(init.A,
                                                                    init.x, 
                                                                    init.gStep*init.g,
                                                                    maximum(init.lastFv),
                                                                    init.funForward,
                                                                    options.funPenalty,
                                                                    init.b,
                                                                    project, 
                                                                    init.timeProject,
                                                                    init.tau,
                                                                    options,
                                                                    params)            
            
            (options.verbosity > 1) && println("fin LineSearch")
            init.nLineTot += nLine
            
            if lnErr == -1
                warn("Line Search Error not set in call to spglinecurvy")
            end

            # If an error was triggered in the line search
            if (lnErr !== 0)
                #DEVNOTE# Finish this if statement
                throw(error("SPGLine Error in development"))

                (options.verbosity > 1) && println("begin FeasLineSearch")

                # Projected backtrack failed. Retry with feasible dir'n line search
                init.x = xOld

                # In-place scaling of gradient and updating of x
                if ~isempty(xOld)
                    dx = project(xOld - init.gStep.*init.g, init.tau, init.timeProject,
                                                            options, params)[1] - xOld
                else
                    throw(error("Empty X")) # This should never be invoked 
                end

                gtd = dot(g,dx)
            end
            
            # Failed again, Revert to Previous iterates and damp max BB step
            if (lnErr !== 0)
                options.stepMax /= 10
                warn("Line Search Failed") #DEVNOTE# Include more info
            end


            doSubspaceMin = false
            if options.subspaceMin
                throw(error("subspaceMin is in development"))
            end

            primNorm_x = options.primal_norm(init.x,options.weights,params)
            targetNorm = init.tau + options.optTol

            if options.ignorePErr

                if primNorm_x > targetNorm
                    warn("Primal norm of projected x is larger than expected")
                end
                
            end

            (options.verbosity > 1) && println("fin UpdateX")

            gOld = copy(init.g) 

            init.f, init.g, init.g2 = funCompositeR(init.A, init.x, init.r, init.funForward,
            options.funPenalty, init.nProdAt, params)

            # xOld plays the role of s
            xOld = init.x - xOld

            y = init.g - gOld
            sts = dot(xOld,xOld)
            sty = dot(xOld,y)

            if sty <= 0
                gStep = options.stepMax
            else
                gStep = min(options.stepMax,max(options.stepMin, sts/sty))
            end

            (options.verbosity > 1) && println("fin CompScaling")


        catch exc
            
            #DEVNOTE# Add in MaxMatVec error catch
            println("""
            ===============================CAUGHT_ERROR===================================== 
            """)
            throw(exc)

            println("""
            ================================================================================ 
            """)
        end

        # ================================================================================
        # Update Function History
        # ================================================================================

        # Dont update if superoptimal 
        if init.singleTau | (init.f > init.sigma)
            init.lastFv[mod(init.iter, options.nPrevVals)+1] = init.f
            if init.fBest > init.f
                fBest = copy(init.f)
                xBest = copy(init.x)
            end
        end
        
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
