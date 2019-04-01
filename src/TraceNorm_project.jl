export TraceNorm_project

"""
Force the rows of L and R to have norm at most B
"""
function TraceNorm_project(x::Tx, B, weights, params::Dict{String,Any}) where {ETx<:Number, Tx<:AbstractVector{ETx}}

if isapprox(norm(x), 0.0)
    println("WARNING: norm(x) cannot be 0 in TraceNorm_project")
end
  
c = sqrt(B/(0.5*norm(x).^2))
xout = min(1,c).*x

#Dummy
itn = 1

return xout, itn

end
