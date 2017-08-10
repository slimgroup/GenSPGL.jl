#GenSPGL

export NormL1_project


#DEVNOTE# Backtrace weights so only oneprojectormex is type tight
"""
Use: x,itn = NormL1_project(x,weights,tau)

"""
function NormL1_project(x::AbstractVector{<:Real}, tau::Number, weights, params)
    
    x_out,itn = oneprojector(x, weights, tau)
    
    return x_out,itn
end



"""
For Complex
"""
function NormL1_project(x::AbstractVector{<:Complex}, tau::Number, weights, params)

    xa = abs.(x)
    idx = find(xa .< eps())
    xc,itn = oneprojector(xa, weights, tau)
    xc = xc ./ xa
    xc[idx] = 0
    x_out = x.*xc 
    
    return x_out,itn
end
