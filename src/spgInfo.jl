# GenSPGL
# spgInfo

export spgInfo

"""
GenSPGL output vars composite type. \n
Passed to spgl1.jl from spglcore.jl

"""
type spgInfo
   
    tau
    rNorm
    gNorm
    rErr
    exit_status
    iter
    nProdA
    nProdAt
    nNewton
    timeProject
    timeMatProd
    options
    xNorm1
    rNorm2
    lambda
end
