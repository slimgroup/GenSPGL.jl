export NormL1_dual

"""
Use: NormL1_dual(x::AbstractArray, weights::Number, params::Dict)
"""
function NormL1_dual(x::AbstractArray, weights::Number, params::Dict)

    d = norm(x/weights, Inf)

end

"""
     NormL1_dual(x::AbstractArray, weights::AbstractArray, params::Dict)
"""
function NormL1_dual(x::AbstractArray, weights::AbstractArray, params::Dict)

    d = norm(x./weights, Inf)

end

