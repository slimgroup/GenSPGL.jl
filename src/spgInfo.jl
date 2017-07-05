# GenSPGL
# spgInfo

export spgInfo

"""
GenSPGL output vars composite type. \n
Passed to spgl1.jl from spglcore.jl

"""
type spgInfo{ETxg<: AbstractFloat, Txg<:AbstractVector{ETxg}}
   
    tau::ETxg
    rNorm::ETxg
    gNorm::ETxg
    rErr::ETxg
    exit_status::spgExitCondition
    iter::Int64
    nProdA::Int64
    nProdAt::Int64
    nNewton::Int64
    timeProject::ETxg
    timeMatProd::ETxg
    options::spgOptions
    xNorm1::Txg
    rNorm2::Txg
    lambda::Txg

end
