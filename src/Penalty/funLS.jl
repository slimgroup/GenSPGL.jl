export funLS

"""
f,g = funLS(r::AbstractArray; params::Dict)\n

Normalize r using the 2 norm. This function uses no params\n

Returns: f: 2-Norm of r\n
         g: Normalized r\n
"""
function funLS{ETr<:Real, Tr<:Union{AbstractArray{ETr}, AbstractArray{Complex{ETr}}}}(
                                                                r::Tr, 
                                                                params::Dict{String,Any})

    f = norm(r, 2)
    g = r/f

    return f, g
end
