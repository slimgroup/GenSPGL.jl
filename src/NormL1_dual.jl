export NormL1_dual

"""
Use: NormL1_dual(x::AbstractArray, weights::Number)
"""
function NormL1_dual(x::AbstractArray, weights::Number)

    d = norm(x/weights, Inf)

end

"""
     NormL1_dual(x::AbstractArray, weights::AbstractArray) 
"""
function NormL1_dual(x::AbstractArray, weights::AbstractArray)

    d = norm(x./weights, Inf)

end

