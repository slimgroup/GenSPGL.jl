## JL-ML SPGL Comparison
export jl_cs

"""
An example of using implicit JOLI operators to solve a system using GenSPGL

Use: x, r, g, info, SNR = jl_cs()
"""
function jl_cs()
    n = 1024
    k = 200

    # Create a system
    p = randperm(n)
    x0 = zeros(n)
    x0[p[1:k]] = sign.(randn(k))

    # Create a modelling op
    F = joDCT(n)
    
    # Create a restriction op
    ind = randperm(n)
    ind = ind[1:Int(floor(0.6*n))]
    R = joRestriction(n,ind, DDT = Float64, RDT = Float64)

    # Create data
    b = R*F*x0

    # Solve
    opts = spgOptions(optTol = 1e-4,
                         verbosity = 1)

    @time x, r, g, info = spgl1(R*F, vec(b), tau = 0., sigma = 1e-3, options = opts) 

    # Calc SNR of recovery
    SNR = -20*log10(norm(x0-x)/norm(x0));
   
    return x, r, g, info, SNR
end
