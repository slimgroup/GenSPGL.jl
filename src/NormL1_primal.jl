export NormL1_primal

"""
     NormL1_primal( x::AbstractArray{<:Number}, 
                    weights<:Union{AbstractVector, Number}, 
                    params::Dict{String, Any})
"""
function NormL1_primal(x::AbstractVector{ETx},
                       weights::Tw,
                       params::Dict{String, Any}) where
                         {Tw<:Union{AbstractVector, Number}, ETx<:Number}

    d = norm(x.*weights, 1)

end

