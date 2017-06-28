#GenSPGL

export NormL1_project


#DEVNOTE# Backtrace weights so only oneprojectormex is type tight
"""
Use: x,itn = NormL1_project(x,weights,tau)

"""
function NormL1_project{ETx<:Number,Tx<:AbstractVector{ETx}}(x::Tx, tau::Number, weights, params)

    println("Script made it into NORML1_project for real x")
    
    x_out::Tx,itn::Int64 = oneprojector(x, weights, tau)
    
    return x_out,itn
end



"""
For Complex
"""
function NormL1_project{Tx<:Complex}(x::AbstractVector{Tx}, tau::Number, weights, params)

    println("Script made it into NORML1_project for complex x")
    
    xa = abs(x)
    idx = find(xa .< eps())
    xc,itn = oneprojector(xa, weights, tau)
    xc = xc ./ xa
    xc[idx] = 0
    x_out = x.*xc 
    
    return x_out,itn
end
