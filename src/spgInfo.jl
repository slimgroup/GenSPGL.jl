# GenSPGL
# spgInfo

export spgInfo

"""
GenSPGL output vars composite type. \n

Tip: Use print(ec::spgExitCondition) to display an exit status.
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
    timeProject::Float64
    timeMatProd::Float64
    options::spgOptions
    xNorm1::Txg
    rNorm2::Txg
    lambda::Txg

end
