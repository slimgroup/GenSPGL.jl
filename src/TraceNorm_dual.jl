export TraceNorm_dual

function TraceNorm_dual{ETx<:Number, Tx<:AbstractVector{ETx},
                        ETw<:Number, Tw<:AbstractArray{ETw}}(
                        x::Tx, weights::Tw, params::Dict{String,Any})

numr::Int = params["numr"]
numc::Int= params["numc"]
nr::Int = params["nr"]
E = reshape(x, numr, numc)

#DEVNOTE# -Performance: Very expensive line
tmp = svds(E; nsv = 1, ritzvec = false)[1]

d = tmp.S[1]

return d

end
