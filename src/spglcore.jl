export spglcore

"""
# Info
    Use:
        spglcore{ETxg<:Number, Txg<:AbstractVector{ETxg}, Tidx<:BitArray}(init::spgInit{ETxg, Txg, Tidx})
# Input
    init:
        - An instance of the composite type spgInit which contains all initalized variables
            created in the preceding call to spgl1.\n

# Outputs
    init: 
        - An instance of the composite type spgInit that contains all updates to variables
            initialized in the preceding call to spgl1\n
    rNorm:
        - Value of the objective\n

# Author
    Keegan Lensink
        Seismic Laboratory for Imaging and Modeling
        The University of British Columbia
        keeganlensink@gmail.com
 
    This code is an adaptation of Michael P. Friedlander, Ewout van den Berg, 
    and Aleksandr Aravkin's MATLAB program SPGL1. 
"""
function spglcore(init::spgInit{TA, ETb, ETx, ETg, ETr, Tidx}) where
            {TA<:Union{joAbstractLinearOperator, AbstractArray, Function},
             ETb<:Number, ETx<:Number, ETg<:Number, ETr<:Number, Tidx<:BitArray}

    #DEVNOTE# Create spgInit type to hold all these initialized paramaters
    
    (init.options.verbosity > 1) && println("script has entered spglcore\n 
            Begin Main Loop.")

    # Pull out options type for easier use
    options = init.options
    params = init.params
    
    # Initialise rErr
    rErr = zero(ETr)

    #Main Loop
    while true

        # ================================================================================
        # Test Exit Conditions
        # ================================================================================ 
        gNorm = zero(ETg)

        if (options.proxy)
            #@code_warntype options.dual_norm(init.g2, options.weights, init.params)
            # prev prod: 1829
            gNorm = options.dual_norm(init.g2, options.weights, init.params)
        else
            gNorm = options.dual_norm(init.g,  options.weights, init.params)
        end
        
        # rNorm and f are the same thing
        rNorm = copy(init.f)
        #println("rNorm = $rNorm")

        tmp_proj, tmp_itn = project(init.x - init.g,
                                        init.tau, init.timeProject, options, params)
        Err = norm(init.x - tmp_proj)
        rErr = Err/max(1,init.f)
   
        aError1 = rNorm - init.sigma
        aError2 = rNorm^2 - init.sigma^2
        rError1 = abs(aError1) / max(1,rNorm)
        rError2 = abs(aError2) / max(1,init.f)

        #= Count number of consecutive iterations with identical support
        
        nnzOld = deepcopy(init.nnzIdx)
        
        #DEVNOTE# -Performance: Expensive function call 
        if ~(options.proxy) 
            nnzX,nnzG,init.nnzIdx,nnzDiff = activevars(init.x, init.g, init.nnzIdx, options, params)
        else
            nnzX,nnzG,init.nnzIdx,nnzDiff = activevars(init.x, init.g2, init.nnzIdx, options, params)
        end

        if (nnzDiff == -1)
            init.nnzIter = 0
        end

        options.verbosity > 1 && println("fin CompConditions")
        
        if isempty(init.x)
            println("DEVNOTE x is empty in spglcore, check spgl1 90-99")
        else
            init.nnzIter += 1
            (init.nnzIter >= options.activeSetIt) && (init.exit_status.triggered = 9)
            tmp_nnzX1 = min(0.1,10*options.optTol)::Float64
            tmp_nnzX2::BitArray{1} = (abs.(init.x) .>= tmp_nnzX1)
            nnzX = sum(tmp_nnzX2)
        end
        =#

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

            if (rErr <= max(options.optTol, rError2)) | (rError1 <= options.optTol)

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

            (false) && println("""
            decTol: $(options.decTol)
            f: $(init.f)
            fOld: $(init.fOld)
            $(abs(init.f - init.fOld)) <= $(options.decTol * init.f)
            testRelChange1: $testRelChange1
            testRelChange2: $testRelChange2
            """)

            init.testUpdateTau = ((testRelChange1 & (rNorm >  2*init.sigma)) |
                            (testRelChange2 & (rNorm <= 2*init.sigma))) &
                            isnull(init.exit_status.triggered) &
                            ~init.testUpdateTau

            if init.testUpdateTau
                
                (options.verbosity > 1) && println("WARNING: Update Tau")
                if (options.quitPareto & (init.iter >= options.minPareto)) 
                    init.exit_status.triggered = 10
                end

                tauOld = copy(init.tau)
                init.tau = max(zero(typeof(init.tau)), init.tau + (aError1)/ gNorm)::typeof(tauOld)
                init.nNewton += one(typeof(init.nNewton))
                
                init.printTau = (abs(tauOld - init.tau) >= (1e-6 * init.tau)) 

                if init.tau < tauOld
                    (options.verbosity > 1) && println("WARNING: Tau Decreasing")
                    init.x, tmp_itn = project(init.x, init.tau, init.timeProject, options, params)

                end

            end

        end

        if (isnull(init.exit_status.triggered) & ((init.iter+1) >= options.iterations))
            init.exit_status.triggered = 5 
        end

        (options.verbosity > 1) && println("fin CheckConverge")

        # ===============================================================================
        # Print log, update history, and check exit conditions
        # ===============================================================================
        

        if (options.verbosity > 0) | init.singleTau | init.printTau |
                                 (init.iter == 0) | ~isnull(init.exit_status.triggered)

            if init.singleTau
                
                s = @sprintf "%5i   %13.7e  %13.7e  %9.2e   %6.1f" init.iter rNorm rErr gNorm log10(init.stepG) #nnzX nnzG
                (options.verbosity > 0) && println(s)

                if init.subspace
                    #println("$itnLSQR")
                end
            else
                
                if init.printTau | init.subspace
                    s = @sprintf "%5i  %13.7e  %13.7e  %9.2e  %9.3e  %6.1f  %13.7e" init.iter rNorm rErr rError1 gNorm log10(init.stepG) #=nnzX nnzG=# init.tau
                else
                    s = @sprintf "%5i  %13.7e  %13.7e  %9.2e  %9.3e  %6.1f" init.iter rNorm rErr rError1 gNorm log10(init.stepG) #nnzX nnzG
                end
                    (options.verbosity > 0) && println(s)
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
        
        # Act on exit conditions
        ~(isnull(init.exit_status.triggered)) && break
        # ================================================================================
        # Begin Iterations
        # ================================================================================
        init.iter += 1
        xOld = copy(init.x)
        init.fOld = copy(init.f)
        rOld = copy(init.r)
       
        try
            # ================================================================================
            # Projected gradient step and line search
            # ================================================================================

            (options.verbosity > 1) && println("begin LineSearch")
            (options.verbosity > 1) && println("Tau: $(init.tau)") 

            init.f, init.x, init.r, nLine, init.stepG, lnErr, localProdA = spglinecurvy(init.A,
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
            init.nProdA += localProdA
            
            if lnErr == -1
                println("WARNING: Line Search Error not set in call to spglinecurvy")
            end

            # If an error was triggered in the line search
            if (lnErr !== 0)

                (options.verbosity > 1) && println("Line Error: $lnErr \n begin FeasLineSearch")

                # Projected backtrack failed. Retry with feasible dir'n line search
                init.x = copy(xOld)

                # In-place scaling of gradient and updating of x
                if ~isempty(xOld)

                    dx = project(xOld - init.gStep.*init.g, init.tau, init.timeProject,
                                                            options, params)[1] - xOld
                else
                    throw(error("Empty X")) # This should never be invoked 
                end

                gtd = dot(init.g,dx)

                if options.linear

                    init.f, step, init.r, nLine, lnErr, localProdA = spgline(init.A,
                                                                init.f,
                                                                dx,
                                                                gtd,
                                                                rOld,
                                                                maximum(init.lastFv),
                                                                init.funForward,
                                                                options.funPenalty,
                                                                params,
                                                                init.b,
                                                                options.feasSrchIt,
                                                                options.linear,
                                                                options,
                                                                init.timeProject)
                else

                    init.f, step, init.r, nLine, lnErr, localProdA = spgline(init.A,
                                                                init.f,
                                                                dx,
                                                                gtd,
                                                                init.x,
                                                                maximum(init.lastFv),
                                                                init.funForward,
                                                                options.funPenalty,
                                                                params,
                                                                init.b,
                                                                options.feasSrchIt,
                                                                options.linear,
                                                                options,
                                                                init.timeProject)
                end

                init.nProdA += localProdA
                
                (options.verbosity > 1) && println("Fin FeasLineSearch")

                if isempty(xOld)
                    init.x = step*dx
                else
                    init.x = xOld + step*dx
                end

                init.x, tmp_itn = project(init.x, init.tau, init.timeProject, options, params)

            end
            
            # Failed again, Revert to Previous iterates and damp max BB step
            if (lnErr !== 0)
                options.stepMax /= 10
                println("WARNING Line Search Failed. Line Error: $(lnErr)") #DEVNOTE# Include more info
            end


            doSubspaceMin = false
            if options.subspaceMin
                throw(error("subspaceMin is in development"))
            end

            primNorm_x = options.primal_norm(init.x,options.weights,params)
            targetNorm = init.tau + options.optTol

            if options.ignorePErr

                if primNorm_x > targetNorm
                    println("WARNING: Primal norm of projected x is larger than expected")
                end
                
            end

            (options.verbosity > 1) && println("fin UpdateX")

            gOld = copy(init.g) 


            init.f, init.g, init.g2 = options.funCompositeR(init.A, init.x, init.r, init.funForward,
            options.funPenalty, init.nProdAt, params)

            # xOld plays the role of s
            xOlds = init.x - xOld

            y = init.g - gOld
            sts = dot(xOlds,xOlds)
            sty = dot(xOlds,y)

            #DEVNOTE# Double check that it is okay to only compare the real part of sty
            if real(sty) <= 0
                init.gStep = options.stepMax
            else
                init.gStep = min(options.stepMax,max(options.stepMin, real(sts/sty)))
            end
            
            (options.verbosity > 1) && println("fin CompScaling")

        catch exc
            
            #DEVNOTE# Add in MaxMatVec error catch
            println("""
            ===============================CAUGHT_ERROR===================================== 
            """)
            println("ERROR: ", exc) 

            display(stacktrace(catch_backtrace()))

            println("""
            ================================================================================ 
            """)
            break
        end

        # ================================================================================
        # Update Function History
        # ================================================================================

        # Dont update if superoptimal 
        if init.singleTau | (init.f > init.sigma)
            init.lastFv[mod(init.iter, options.nPrevVals)+1] = init.f
            if init.fBest > init.f
                init.fBest = copy(init.f)
                init.xBest = copy(init.x)
            end
        end
        
    end #Main Loop
    
    # Restore best solution(only if solving single problem)
    if (init.singleTau) & (init.f >init.fBest)
        if options.restore
            rNorm = init.fBest
            (options.verbosity > 1) && println("Restoring best iterate to objective: $(rNorm)")
            init.x = init.xBest
            init.r = init.b - init.funForward(init.A, init.x, Array{ETx,1}(), params)
            init.f, init. g, init.g2 = options.funCompositeR(init.A, init.x, init.r, init.funForward,
                                            options.funPenalty, init.nProdAt, params)
            if options.proxy
                gNorm = options.dual_norm(init.g2, options.weights)
            else
                gNorm = options.dual_norm(init.g, options.weights)
            end
            rNorm = init.f
        else
            (options.verbosity > 0) && println("""
                Note: Solution not actually optimal. Best objective value is $(init.fBest)
                """)
        end
    end
    
    return init, rNorm, gNorm, rErr
end



"""
Use: nnzX,nnzG,nnzIdx,nnzDiff = activevars(x,g,nnzIdx,options, params)

Find the current active set.
nnzX    is the number of nonzero x.
nnzG    is the number of elements in nnzIdx.
nnzIdx  is a vector of primal/dual indicators.
nnzDiff is the no. of elements that changed in the support.
"""
function activevars(x::Txg,g::Txg, nnzIdx::Ti, options::spgOptions,params::Dict{String,Any}) where
            {Ti<:BitArray{1}, ETxg<:Number, Txg<:AbstractVector{ETxg}}

    xTol = min(.1,10*options.optTol)
    gTol = min(.1,10*options.optTol)

    #DEVNOTE# -Performance: Expensive Line
    gNorm = options.dual_norm(g,options.weights, params)
    
    if isnull(nnzIdx)
        nnzOld_inner = Ti()
    else
        nnzOld_inner = copy(nnzIdx)
    end

    #Reduced costs for postive and negative parts of x
    z1 = gNorm .+ g
    z2 = gNorm .- g

    #Primal/dual based indicators
    xPos = BitArray{1}()
    xNeg = BitArray{1}()

    if(~options.proxy)
        # Unwrap temp vars for performance
        a = x .> xTol
        b = z1 .< gTol
        xPos = a .& b
        
        c = x .< -xTol
        d = z2 .< gTol
        xNeg = c .& d
        
        nnzIdx = xPos .| xNeg
    end
    
    nnzX = sum(abs.(x) .>= xTol)::Int64
    nnzG = sum(nnzIdx)::Int64

    if isempty(nnzOld_inner)
        #DEVNOTE# Use -1 instead of Inf? Inf creates a type stability
                # Need to document what this means though
        nnzDiff = -one(Int64)
    else
        e = nnzIdx .!== nnzOld_inner
        nnzDiff = sum(e)::Int64
    end

    return nnzX, nnzG, nnzIdx, nnzDiff
end

mtypeof(c...) = [typeof(a) for a in c]
