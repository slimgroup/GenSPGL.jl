# GenSPGL
# spgInit

export spgInit

"""
GenSPGL Initialized local vars composite type. \n
Passed to spglcore.jl

"""
type spgInit{TA<:Union{joAbstractLinearOperator,AbstractArray, Function},
                                                ETb<:Number,
                                                ETx<:Number,
                                                ETg<:Number,
                                                ETr<:Number,
                                                Tidx<:BitArray}
    
    A::TA
    b::AbstractVector{ETb}
    x::AbstractVector{ETx}
    tau::Float64
    sigma::Number
    g::AbstractVector{ETg}
    g2::AbstractVector{ETg}
    f::Float64
    r::AbstractVector{ETr}
    nnzIdx::Tidx
    nnzIter::Int64
    options::spgOptions
    params::Dict{String,Any}
    timeProject::Float64
    timeMatProd::Float64
    exit_status::spgExitCondition  
    singleTau::Bool
    bNorm::AbstractFloat
    fOld::Float64
    testUpdateTau::Bool
    iter::Int64
    nNewton::Int64
    printTau::Bool
    subspace::Bool
    stepG::Float64
    xNorm1::Vector{Float64}
    rNorm2::Vector{Float64}
    lambda::Vector{Float64}
    lastFv::Vector{Float64} # Has to be Float64 to support Inf
    gStep::Float64
    funForward::Function
    nLineTot::Int64
    nProdAt::Int64
    nProdA::Int64
    fBest::Float64
    xBest::AbstractVector{ETx}

end
