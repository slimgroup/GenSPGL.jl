import Base.sortperm

export oneProjector, oneProjectorMex, sortperm_col #DEVNOTE# Don't need to expose ...Mex after debugging

"""
GenSPGL version of oneProjector.m\n

See MATLAB version documents below for info on use.\n

================================================================================

ONEPROJECTOR  Projects b onto the weighted one-norm ball of radius tau

    [X,ITN] = ONEPROJECTOR(B,TAU) returns the orthogonal projection
    of the vector b onto the one-norm ball of radius tau. The return
    vector X which solves the problem

            minimize  ||b-x||_2  st  ||x||_1 <= tau.
               x

    [X,ITN] = ONEPROJECTOR(B,D,TAU) returns the orthogonal
    projection of the vector b onto the weighted one-norm ball of
    radius tau, which solves the problem

            minimize  ||b-x||_2  st  || Dx ||_1 <= tau.
               x

    If D is empty, all weights are set to one, i.e., D = I.

    In both cases, the return value ITN given the number of elements
    of B that were thresholded.

 See also spgl1.

   oneProjector.m
   oneProjector.m 1200 2008-11-21 19:58:28Z mpf 
================================================================================
"""
#DEVNOTE# Eventually use multiple dispatch to support scalar weights
function oneProjector(b::AbstractArray, d::AbstractArray, tau::AbstractFloat)

    println("Script made it into oneProjector") 

    ~(length(d)==1) && ~(length(b) == length(d)) && println("""
    Vectors 'b' and 'd' must be the same length
    """)
    
    # Quick return for the easy case
    if (length(d)==1) & (d[1] == 0) 
        x = b
        itn = 0
        return x, itn
    end

    # Get sign of b and set to absolute values
    b_abs = abs(b)

    # Perform projection
    if length(d)==1
        x,itn = oneProjectorMex(b_abs, d[1], tau)
        return x,itn
    else
        d_abs = abs(d)
        idx = find(d .> eps())
        x = deepcopy(b_abs) #DEVNOTE# Double check avoiding referencing is necessary
        x[idx],itn = oneProjectorMex(b_abs[idx], d[idx], tau)
        return x,itn
    end
end



"""
oneProjectorMex for scalar weight
oneProjectorMex_I clone
"""
function oneProjectorMex{T<:Number}(b::AbstractArray{T}, d::Number, tau::AbstractFloat)

    println("Script made it into oneProjectorMex for scalar weight")
    
    tau = tau/abs(d)

    #Initialization
    n = length(b)
    x = zeros(T,n,1)
    bNorm = norm(b,1)

    #Check for quick exit
    (tau >= bNorm) && (x=b; itn=0; return x,itn)
    (tau < eps()) && (itn = 0; return itn)

    # Preprocessing (b is assumed to be >= 0)
    idx = sortperm(b, rev=true)
    b_sort = b[idx]

    csb = -tau
    alphaPrev = 0

    for j = 1:n
        csb += b[j]
        alpha = csb/j

        # Finish as soon as constraint can be satisfied w/o exceeding current min val of b
        (alpha >= b[j]) && break

        alphaprev = alpha

    end

    # Set the solution by apply soft-thresholding with previous value of alpha
    x[idx] = max(0, b .- alphaPrev)

    # Set number of iterations
    itn = j

    return x, itn


end



"""
#DEVNOTE# Don't need this right now, maybe not ever
This function is part of GenSPGL

Use: sortperm_rev(A::AbstractMatrix)

Returns the indicies for sorted columns of A
"""
function sortperm_col(A::AbstractMatrix; rev::Bool = false)

    n,m = size(A)

    # Init idx
    idx = zeros(Int64,n,m)
   
    # Loop over columns
    for i = 1:m
        idx[:,i] = sortperm(A[:,i], rev = rev)
    end

    return idx
end
            
