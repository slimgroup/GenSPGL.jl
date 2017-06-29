# GenSPGL
# spgInit

export spgInit

"""
GenSPGL Initialized local vars composite type. \n
Passed to spglcore.jl

"""
type spgInit{Txg<:Float64,Tidx<:BitArray}

    x::AbstractVector{Txg}
    tau::Txg
    sigma::Txg
    g::AbstractVector{Txg}
    g2::Txg
    f::Txg
    nnzIdx::Tidx
    nnzIter::Int64
    options::spgOptions
    params::Dict{String,Number}
    timeProject::Txg
    exit_status::spgExitCondition  
    singleTau::Bool
    bNorm::Txg
    fOld::Txg
    testUpdateTau::Bool
    iter::Int64
    nNewton::Int64
end
