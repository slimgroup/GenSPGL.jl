#GenSPGL

export NormL1_project

"""
Use: x = NormL1_project(x,weights,tau)

"""
function NormL1_project(x::AbstractArray, tau::Number; weights::AbstractArray = [1])

    println("Script made it into NORML1_project")

    if isreal(x)
        x,itn = oneProjector(x, weights, tau)
    else
        xa = abs(x)
        idx = find(xa .< eps())
        xc,itn = oneProjector(xa, weights, tau)
        xc = xc ./ xa
        xc[idx] = 0
        x = x.*xc 
    end
    

    println("""
    ================================================================================
    NormL!_project results:
    x = \n
    $x\n

    itn = \n
    $itn

    """)
    return x,itn
end
