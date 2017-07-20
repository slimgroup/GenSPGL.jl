# Load a common workspace and solve a lasso problem
using GenSPGL
file = MAT.matopen("vars_complex.mat")

b = read(file,"b")
tau = read(file,"tau")
xinit = read(file,"xinit")
sigma = read(file,"sigma")

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

# Avoid anon func
afunT(x) = reshape(x,355,709)
params = Dict{String, Any}([("nr", 40)
                                ("Ind", vec(b) .== 0)
                                ("numr", 355)
                                ("numc", 709)
                                ("funForward", NLfunForward)
                                ("afunT", afunT)
                                ("afun", afun)
                                ("mode", 1)
                                ("ls", 1)
                                ("logical", 0)
                                ("funPenalty", funLS)])

@time x_LS, r, g, info = spgl1(NLfunForward, vec(b), x = vec(xinit),
                                                       tau = tau,
                                                       sigma = sigma, 
                                                       options = opts,
                                                       params = params)
fid2 = MAT.matopen("mat_complex_sol.mat")
mat_xLS = vec(read(fid2, "xLS"))

