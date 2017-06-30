# GenSPGL
# spgInit

export spgInit

"""
GenSPGL Initialized local vars composite type. \n
Passed to spglcore.jl

"""
type spgInit{ETxg<:Float64, Txg<:AbstractVector{ETxg}, Tidx<:BitArray}

    x::Txg
    tau::ETxg
    sigma::ETxg
    g::Txg
    g2::ETxg
    f::ETxg
    r::Txg
    nnzIdx::Tidx
    nnzIter::Int64
    options::spgOptions
    params::Dict{String,Number}
    timeProject::ETxg
    exit_status::spgExitCondition  
    singleTau::Bool
    bNorm::ETxg
    fOld::ETxg
    testUpdateTau::Bool
    iter::Int64
    nNewton::Int64
    printTau::Bool
    subspace::Bool
    stepG::ETxg
    xNorm1::Vector{ETxg}
    rNorm2::Vector{ETxg}
    lambda::Vector{ETxg}
end
