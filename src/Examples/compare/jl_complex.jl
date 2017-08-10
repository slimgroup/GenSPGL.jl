# Load a common workspace and solve a lasso problem
export jl_complex

function jl_complex()
    exampledir = Pkg.dir("GenSPGL")*"/src/Examples/compare/"
    file = MAT.matopen(exampledir*"vars_complex.mat")

    b = read(file,"b")
    tau = read(file,"tau")
    test = read(file,"test")
    xinit = read(file,"xinit")
    sigma = read(file,"sigma")
    close(file)

    # Choose options for GenSPGL
    opts = spgOptions(  optTol = 1e-5,
                        bpTol = 1e-5,
                        decTol = 1e-4,
                        project = TraceNorm_project,
                        primal_norm = TraceNorm_primal,
                        dual_norm = TraceNorm_dual,
                        proxy = true,
                        ignorePErr = true,
                        iterations = 150,
                        verbosity = 1)

    # Create Params Dict
    afunT(x) = reshape(x,355,709)
    params = Dict("nr"=> 40,
                    "Ind"=> vec(b) .== 0,
                    "numr"=> 355,
                    "numc"=> 709,
                    "funForward"=> NLfunForward,
                    "afunT"=> afunT,
                    "afun"=> afun,
                    "mode"=> 1,
                    "ls"=> 1,
                    "logical"=> 0,
                    "funPenalty"=> funLS)

    # Set threads and solve problems
    #BL.set_num_threads(8)
    @time xLS_jl, r, g, info = spgl1(NLfunForward, vec(b), x = vec(xinit),
                                                           tau = tau,
                                                           sigma = sigma, 
                                                           options = opts,
                                                           params = params)

    # Write output for MATLAB comparison
    fw = MAT.matopen(exampledir*"xLS_jl.mat", "w")
    write(fw, "xLS_jl", xLS_jl)
    close(fw)
    
    return xLS_jl, r, g, info
end
