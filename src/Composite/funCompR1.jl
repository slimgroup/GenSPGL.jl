export funCompR1

"""
GenSPGL (this function is used for traditional spgl1)

Use:    f,g1,g2 = funCompositeR(A, r, funForward, funPenalty, params)
"""
function funCompR1(A::TA,
                       x::AbstractArray,
                       r::AbstractArray,
                       funForward::Function, funPenalty::Function, 
                       nProdAt::Int64,
                       params::Dict{String,Any}) where
                         {TA<:Union{joAbstractLinearOperator,AbstractArray,Function}}

    nProdAt += one(Int64)
    f,v = funPenalty(r, params)
    
    if ~(params["proxy"])
        g1 = funForward(A, x, -v, params)
        g2 = [zero(eltype(g1))]
    else
        g1,g2 = funForward(A, x, -v, params)
    end
    
    return f,g1,g2
end
