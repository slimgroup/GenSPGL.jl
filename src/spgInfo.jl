# GenSPGL
# spgInfo

export spgInfo

"""
GenSPGL output vars composite type. \n

Tip: Use print(ec::spgExitCondition) to display an exit status.
"""
mutable struct spgInfo
   
    tau::Number
    rNorm::Number
    gNorm::Number
    rErr::Number
    exit_status::spgExitCondition
    iter::Int64
    nProdA::Int64
    nProdAt::Int64
    nNewton::Int64
    timeProject::Float64
    timeMatProd::Float64
    options::spgOptions
    xNorm1::Vector{Number}
    rNorm2::Vector{Number}
    lambda::Vector{Number}

end
