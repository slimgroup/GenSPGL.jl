export TraceNorm_project

"""
Force the rows of L and R to have norm at most B
"""
function TraceNorm_project{ETx<:Number, Tx<:AbstractVector{ETx}}(x::Tx, B, weights, params::Dict{String,Any})

if isapprox(norm(x), 0.0)
    warn("norm(x) cannot be 0 in TraceNorm_project")
end
  
c = sqrt(B/(0.5*norm(x).^2))
xout = min(1,c).*x

#Dummy
itn = 1

return xout, itn

end
