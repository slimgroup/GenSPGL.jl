export TraceNorm_primal

"""
Force the rows of L and R to have norm at most B
"""
function TraceNorm_primal(x::Tx, weights, params::Dict{String,Any}) where {ETx<:Number, Tx<:AbstractVector{ETx}}

p = 0.5 * norm(x.*weights)^2

return p

end
