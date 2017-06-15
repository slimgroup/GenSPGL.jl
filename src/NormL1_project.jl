#GenSPGL

export NormL1_project

"""
Use: x = NormL1_project(x,weights,tau)

"""
function NormL1_project(x::AbstractArray, tau::Number; weights::AbstractArray = [1])

    println("Script made it into NORML1_project")

    if isreal(x)
        x = oneProjector(x,weights,tau)
    else #DEVNOTE# Not done here
        println("NormL1_project doesnt support non-real x yet")
    end

end
