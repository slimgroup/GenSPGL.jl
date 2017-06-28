# GenSPGL
# spgInit

export spgInit

"""
GenSPGL Initialized local vars composite type. \n
Passed to spglcore.jl

"""
type spgInit{Txg<:Float64, T<:Int64}

    x::AbstractVector{Txg}
    tau::Txg
    sigma::Txg
    g::AbstractVector{Txg}
    g2::Txg
    f::Txg
    nnzIdx::T
    options::spgOptions
    params::Dict
    timeProject::Txg
    
end
