export oneprojector, oneprojectormex, sortperm_col #DEVNOTE# Don't need to expose ...Mex after debugging

"""
GenSPGL version of oneprojector.m\n

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

   oneprojector.m
   oneprojector.m 1200 2008-11-21 19:58:28Z mpf 
================================================================================
"""
function oneprojector(b::AbstractArray, d, tau::AbstractFloat)

    
    len_d = length(d)
    len_b = length(b)

    ~(len_d==1) && ~(len_b == len_d) && println("""
    Vectors 'b' and 'd' must be the same length
    Length b: $(len_b)
    Length d: $(len_d)
    
    """)
   
    # Declare x for stability
    x::typeof(b) = similar(b)
    itn::Int = zero(Int)
    #DEVNOTE# Not necessary if there is a scalar method 
    # Quick return for the easy case
    if (len_d==1) & (d[1] == 0) 
        x = b
        itn = 0
        return x, itn
    end

    # Get sign of b and set to absolute values
    s = sign.(b)
    b_abs = abs.(b)

    # Perform projection
    if len_d==1
        x,itn = oneprojectormex(b_abs, d[1], tau)

    else
        
        d_abs = abs.(d)
        idx = find(d .> eps())
        x = deepcopy(b_abs) 
        x[idx],itn = oneprojectormex(b_abs[idx], d[idx], tau)
    end
   
    # Restore signs of x
    x .*= s
    return x,itn
end



"""
oneprojectormex for scalar weight
oneprojectormex_I clone
"""
function oneprojectormex(b::AbstractVector{T}, d::Number, tau::Number) where {T<:Number}

    tau = tau/abs.(d)
    len_b = length(b)
    
    #Initialization
    n = len_b
    x = zeros(T,n)
    bNorm = norm(b,1)

    #Check for quick exit
    (tau >= bNorm) && (x=b; itn=0; return x,itn)
    (tau < eps()) && (itn = zero(Int64); return x,itn)
    
    # Preprocessing (b is assumed to be >= 0)
    idx = sortperm_col(b, rev=true)
    b_sort = b[idx]

    csb = -tau
    alphaprev = zero(T)

    j_out = 1
    for j = 1:n
        csb += b_sort[j]
        alpha = csb/j

        # Finish as soon as constraint can be satisfied w/o exceeding current min val of b
        (alpha >= b_sort[j]) && break

        alphaprev = alpha

        j_out = j
    end
    
    # Set the solution by apply soft-thresholding with previous value of alpha
    x[idx] = max.(0, b_sort .- alphaprev)

    # Set number of iterations
    itn = j_out

    return x, itn

end



"""
Use: x,itn = oneprojectormex(b::Abstractvector, d::AbstractVector, tau::Number)
"""
function oneprojectormex(b::AbstractVector{<:Number}, d::AbstractVector{<:Number}, tau::Number)
    

    #Get type of b.*d
    Tdb = promote_type(eltype(b), eltype(d))

    len_d = length(d)
    len_b = length(b)
    
    #Check for quick exit
    (tau >= norm(d.*b,1)) && (x=b; itn= 0; return x,itn)
    (tau < eps()) && (itn = 0; return x,itn)

    n = len_b
    x = zeros(Tdb,n,1)
   
    # Preprocessing
    bd = b./d
    idx = sortperm_col(bd, rev = true)
    b_sort = b[idx]
    d_sort = d[idx]
    bd_sort = bd[idx]
    
    # Optimize
    csdb = zero(Tdb)
    csd2 = zero(Tdb)
    soft = zero(Tdb)
    alpha1 = zero(Tdb)
    soft = zero(Tdb)
    i = 1

    while i <= n
        csdb += d_sort[i].*b_sort[i]
        csd2 = csd2 + d[i].*d[i]

        alpha1 = (csdb - tau)/ csd2
        alpha2 = bd_sort[i]

        (alpha1 >= alpha2) && break

        soft = alpha1
        i += 1
    end

    x[idx[1:i-1]] = b_sort[1:i-1] - d_sort[1:i-1]*max(0,soft)
    itn = i

    return x,itn
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
            
function sortperm_col(A::AbstractVector; rev::Bool = false)

    idx = sortperm(A, rev = rev)

    return idx
end
            
