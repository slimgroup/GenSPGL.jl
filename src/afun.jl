export afun

function afun{ETx<:Number, Tx<:AbstractArray{ETx}}(x::Tx, params::Dict{String, Any})

if params["logical"] == 1
    xv = vec(x)
    Ind = params["Ind"]
    xv[.~Ind] = zero(ETx)
else
    xv = vec(x)
    Ind = params["Ind"]
    xv[Ind] = zero(ETx)
end

return xv
end
