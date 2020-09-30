export funCompR2

"""
GenSPGL (this function is used for parallel_weighted spgl1, it is a traditional spgl2 function)

Use:    f,g1,g2 = funCompositeR(A, r, funForward, funPenalty, params)
"""
function funCompR2(A::TA,
                       x::AbstractArray,
                       r::AbstractArray,
                       funForward::Function, funPenalty::Function, 
                       nProdAt::Int64,
                       params::Dict{String,Any}) where {TA<:InType}

    nProdAt += one(Int64)
    f,v = funPenalty(r, params)
    
    if ~(params["proxy"])
        g1 = funForward(A, x, -v, params)
        g2 = g1[:]
    else
        g1 = funForward(A, x, -v, params)
        g2 = g1[:]
    end
    return f,g1,g2
end
