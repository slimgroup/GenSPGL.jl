# GenSPGL
# spgInit

export spgInit

"""
GenSPGL Initialized local vars composite type. \n
Passed to spglcore.jl

"""
type spgInit{Txg<:Number}

    x::AbstractVector{Txg}
    tau::Number
    sigma::Number
    g::AbstractVector{Txg}
    g2::Number
    f::Number
    nnzIdx::AbstractArray
    options::spgOptions
    params::Dict
    timeProject::AbstractArray
    
end
