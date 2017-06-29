export NormL1_dual

"""
Use: NormL1_dual(x::AbstractArray, weights::Number, params::Dict)
"""
function NormL1_dual{ETx<:AbstractFloat}(x::AbstractVector{ETx}, weights::Number, params::Dict)::ETx

    d = norm(x/weights, Inf)

end

"""
     NormL1_dual(x::AbstractArray, weights::AbstractArray, params::Dict)
"""
function NormL1_dual{ETx<:AbstractFloat}(x::AbstractVector{ETx},
                                    weights::AbstractVector, params::Dict)::ETx

    d = norm(x./weights, Inf)

end

