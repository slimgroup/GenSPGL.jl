export spglinecurvy

"""
Use: fNew,xNew,rNew,iter,step,err, nProd = spglinecurvy(x,g,fMax,funForward, funPenalty, b,project,tau, params)\n
Projected backtracking linesearch. On entry, g is the (possibly scaled) steepest
descent direction. 
"""
function spglinecurvy{TA<:Union{joAbstractLinearOperator, AbstractArray},
                                                        Tf<:Number,
                                                        ETx<:Number,
                                                        ETg<:Number,
                                                        ETb<:Number}(
                                    A::TA, 
                                    x::AbstractArray{ETx}, 
                                    g::AbstractArray{ETg}, 
                                    fMax::Tf, 
                                    funForward::Function,
                                    funPenalty::Function,
                                    b::AbstractVector{ETb}, 
                                    project::Function,
                                    timeProject::Float64,
                                    tau::Float64, 
                                    options::spgOptions,
                                    params::Dict{String,Any})

    (options.verbosity > 1) && println("Script entered spglinecurvy for A::AbstractArray")

    nProd = 0
    nProdAt = 0
    EXIT_CONVERGED = 0
    EXIT_ITERATIONS = 1
    EXIT_NODESCENT = 2
    gamma = 1e-4
    maxIts = 10
    step = one(real(ETg))
    sNorm = zero(real(ETx))
    scale = one(real(ETg))
    nSafe = 0
    iter = 0
    debug = false
    n = length(x)

    # Define outputs outside of while loop scope
    xNew = similar(x)
    rNew = Array{ETb,1}()
    fNew = zero(Tf)
    err = -1
    while true

        xNew, tmp_itr = project(x - step*scale*g, tau, timeProject, options, params)
        rNew = b - funForward(A, xNew, [], params)::Array{ETb,1}
        nProd += 1
        fNew, dummy_g = funPenalty(rNew, params)

        s = xNew - x
        gts = scale * real(g'*s)

        if gts >= 0
            err = EXIT_NODESCENT
            break
        end

        if fNew < fMax + gamma*step*gts
            err = EXIT_CONVERGED
            break
        elseif iter >= maxIts
            err = EXIT_ITERATIONS
            break
        end

        #New linesearch iteration
        iter += 1
        step /= 2

        # Safeguard: If stepMax is huge, then even damped search
        # directions can give exactly the same point after projection.  If
        # we observe this in adjacent iterations, we drastically damp the
        # next search direction.

        sNormOld = sNorm
        sNorm = norm(s) / sqrt(n)
        
        if abs(sNorm - sNormOld) <= 1e-6*sNorm
            gNorm = norm(g) / sqrt(n)
            scale = sNorm/gNorm/(2^nSafe)
            nSafe+=1
        end

    end

    return fNew, xNew, rNew, iter, step, err, nProd 
end

function spglinecurvy{TA<:Function, Tf<:Number, ETx<:Number,
                                                        ETg<:Number,
                                                        ETb<:Number}(
                                    A::TA, 
                                    x::AbstractArray{ETx}, 
                                    g::AbstractArray{ETg}, 
                                    fMax::Tf, 
                                    funForward::Function,
                                    funPenalty::Function,
                                    b::AbstractVector{ETb}, 
                                    project::Function,
                                    timeProject::Float64,
                                    tau::Float64, 
                                    options::spgOptions,
                                    params::Dict{String,Any})

    (options.verbosity > 1) && println("Script entered spglinecurvy for A::Function")

    nProd = 0
    nProdAt = 0
    EXIT_CONVERGED = 0
    EXIT_ITERATIONS = 1
    EXIT_NODESCENT = 2
    gamma = 1e-4
    maxIts = 10
    step = one(real(ETg))
    sNorm = zero(real(ETx))
    scale = one(real(ETg))
    nSafe = 0
    iter = 0
    debug = false
    n = length(x)

    # Define outputs outside of while loop scope
    xNew = similar(x)
    rNew = Array{ETb,1}()
    fNew = zero(Tf)
    err = -1
    while true

        xNew, tmp_itr = project(x - step*scale*g, tau, timeProject, options, params)
        tmp1::Array{ETb,1} = funForward(A, xNew, [], params)[1]
        rNew::Array{ETb,1} = b - tmp1
        nProd += 1

        #DEVNOTE# Expensive Line
        fNew, dummy_g = funPenalty(rNew, params)

        s = xNew - x
        gts = scale * real(g'*s)

        if gts >= 0
            err = EXIT_NODESCENT
            break
        end

        if fNew < fMax + gamma*step*gts
            err = EXIT_CONVERGED
            break
        elseif iter >= maxIts
            err = EXIT_ITERATIONS
            break
        end

        #New linesearch iteration
        iter += 1
        step /= 2

        # Safeguard: If stepMax is huge, then even damped search
        # directions can give exactly the same point after projection.  If
        # we observe this in adjacent iterations, we drastically damp the
        # next search direction.

        sNormOld = sNorm
        sNorm = norm(s) / sqrt(n)
        
        if abs(sNorm - sNormOld) <= 1e-6*sNorm
            gNorm = norm(g) / sqrt(n)
            scale = sNorm/gNorm/(2^nSafe)
            nSafe+=1
        end

    end

    return fNew, xNew, rNew, iter, step, err, nProd 
end

