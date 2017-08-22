## JL-ML SPGL Comparison
export scale_cs

"""
An example of using implicit JOLI operators to solve a system using GenSPGL

Use: x, r, g, info, SNR = ex_cs()
"""
function scale_cs(n,k,i)

    #Create a sparse vector
    #inds = randperm(n)
    #x0 = zeros(n)
    #keep = Int(floor(n/1000))
    #x0[inds[1:keep]] = randn(keep)
    #r_inds = inds[1:k]

    #= Load sparse vector x and restriction mask
    fid = matopen(Pkg.dir("GenSPGL")*"/src/Examples/scaletest/x0.mat")
    x0 = vec(read(fid, "x0"))
    r_inds = Int.(vec(read(fid, "r_inds")))
    n = length(vec(x0))
    close(fid)
    =#

    # Load sparse vector x and restriction mask
    fid = matopen(Pkg.dir("GenSPGL")*"/src/Examples/scaletest/x0.mat")
    x0 = vec(read(fid, "x0"))
    r_inds = Int.(vec(read(fid, "r_inds")))
    n = length(vec(x0))
    close(fid)

    # Sparsity Promoting Transform
    S = joDCT(n)
    
    # Restriction Operator
    R = joRestriction(n, r_inds)
    
    # Modelling op
    A = R*S

    # Create data
    b = A*x0

    # Set BLAS Threads
    BLAS.set_num_threads(8)

    # Solve
    opts = spgOptions(verbosity = 1, iterations = 100)
    t = [@elapsed spgl1(A, b, options = opts) for ii = 1:i]
    return t

end
