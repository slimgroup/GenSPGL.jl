export TraceNorm_project

"""
Force the rows of L and R to have norm at most B
"""
function TraceNorm_project{ETx<:Number, Tx<:AbstractVector{ETx}}(x::Tx, B, weights, params::Dict{String,Any})
numr::Int = params["numr"]
numc::Int= params["numc"]
nr::Int = params["nr"]
e = numr * nr

L_v::Vector{ETx} = x[1:e]
R_v::Vector{ETx} = x[(e+1):end]

L::Array{ETx,2} = reshape(L_v, numr, nr)
R::Array{ETx,2} = reshape(R_v, numc, nr)

c = sqrt(B/(0.5*norm(x).^2))

Lout = min(1,c)*vec(L)
Rout = min(1,c)*vec(R)

xout = append!(vec(Lout), vec(Rout))

#Dummy
itn = 1

return xout, itn

end
