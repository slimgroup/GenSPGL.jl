export oneProjector

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
function oneProjector(b::AbstractArray, d::AbstractArray, tau::AbstractFloat)

    println("Script made it into oneProjector") 

    ~(length(d)==1) && ~(length(b) == length(d)) && println("""
    Vectors 'b' and 'd' must be the same length
    """)

end
