export TraceNorm_dual

function TraceNorm_dual{ETx<:Number, Tx<:AbstractVector{ETx}}(x::Tx, weights, params::Dict{String,Any})

E = reshape(x, params["numr"], params["numc"])

tmp = svds(E; nsv = 1)[1]

d = tmp.S[1]

return d

end
