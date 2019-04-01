export afun

function afun(x::Tx, params::Dict{String, Any}) where {ETx<:Number, Tx<:AbstractArray{ETx}}

if params["logical"] == 1
    xv = vec(x)
    Ind = params["Ind"]
    xv[.~Ind] = zero(ETx)
else
    xv = vec(x)
    Ind = params["Ind"]
    xv[Ind] .= zero(ETx)
end

return xv
end
