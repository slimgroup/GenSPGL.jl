## JL-ML SPGL Comparison
export ex_cs

"""
An example of using implicit JOLI operators to solve a system using GenSPGL

Use: x, r, g, info, SNR = ex_cs()
"""
function ex_cs()
    #n = 1024
    #k = 50

    # Create a sparse vector
    #inds = randperm(n)
    #x0 = zeros(n)
    #x0[inds[1:10]] = randn(10)

    # Load sparse vector x and restriction mask
    fid = matopen(Pkg.dir("GenSPGL")*"/src/Examples/pres/x0.mat")
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
    opts = spgOptions(verbosity = 1)
    x, r, g, info = spgl1(A, b, options = opts) 

    # Solve again silently for timing
    opts = spgOptions(verbosity = 0)
    t_jl = @elapsed spgl1(A, b, options = opts)

    # Calc SNR of recovery
    SNR = -20*log10(norm(x0-x)/norm(x0));

    # Write results
    fw = matopen("x_jl.mat", "w")
    write(fw, "x_jl", x)
    write(fw, "t_jl", t_jl)
    close(fw)
   
    return x, r, g, info, SNR
end
