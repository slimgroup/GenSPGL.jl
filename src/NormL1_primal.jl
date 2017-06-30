export NormL1_primal

"""
Use: NormL1_primal(x::AbstractArray, weights::Number, params::Dict)
"""
function NormL1_primal{ETx<:AbstractFloat}(x::AbstractVector{ETx}, weights::Number, params::Dict)::ETx

    d = norm(x*weights, Inf)

end

"""
     NormL1_primal(x::AbstractArray, weights::AbstractArray, params::Dict)
"""
function NormL1_primal{ETx<:AbstractFloat}(x::AbstractVector{ETx},
                                    weights::AbstractVector, params::Dict)::ETx

    d = norm(x.*weights, Inf)

end

