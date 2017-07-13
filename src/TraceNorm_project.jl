export TraceNorm_project

"""
Force the rows of L and R to have norm at most B
"""
function TraceNorm_project{ETx<:Number, Tx<:AbstractVector{ETx}}(x::Tx, B, weights, params::Dict{String,Any})

e = params["numr"] * params["nr"]

L = x[1:e]
R = x[(e+1):end]

L = reshape(L, params["numr"], params["nr"])
R = reshape(R, params["numc"], params["nr"])

c = sqrt(B/(0.5*norm(x).^2))
Lout = min(1,c)*vec(L)
Rout = min(1,c)*vec(R)

xout = append!(vec(Lout), vec(Rout))

#Dummy
itn = 1

return xout, itn

end
