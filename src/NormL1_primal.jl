export NormL1_primal

"""
     NormL1_primal( x::AbstractArray{<:Number}, 
                    weights<:Union{AbstractVector, Number}, 
                    params::Dict{String, Any})
"""
function NormL1_primal{Tw<:Union{AbstractVector, Number}, ETx<:Number}(
                                    x::AbstractVector{ETx},
                                    weights::Tw,
                                    params::Dict{String, Any})

    d = norm(x.*weights, 1)

end

