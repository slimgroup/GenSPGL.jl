# GenSPGL
# spgInit

export spgInit

"""
GenSPGL Initialized local vars composite type. \n
Passed to spglcore.jl

"""
type spgInit{ETxg<:Float64, Txg<:AbstractVector{ETxg}, Tidx<:BitArray}
    
    A::AbstractArray
    b::AbstractVector
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
    timeProject::Float64
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
    lastFv::Vector{Float64} # Has to be Float64 to support Inf
    gStep::ETxg
    funForward::Function
    nLineTot::Int64
    nProdAt::AbstractArray
end
