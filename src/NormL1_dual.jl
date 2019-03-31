export NormL1_dual

"""
     NormL1_dual(x::AbstractArray, weights::AbstractArray, params::Dict)
"""
function NormL1_dual(x::AbstractVector{ETx},
                     weights::Tw,
                     params::Dict{String, Any}) where
                       {Tw<:Union{AbstractVector, Number}, ETx<:Number}

    d = norm(x./weights, Inf)

end

