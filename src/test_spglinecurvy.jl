# Script to test spglinecurvy

include("spglinecurvy.jl")

using GenSPGL
using MAT
var = matopen("/home/slim/klensink/.julia/v0.6/GenSPGL/compare/testspglc.mat")

gStep= read(var,"gStep")
tau = read(var,"tau")
lastFv = vec(read(var,"lastFv"))
b = vec(read(var,"b"))
x = vec(read(var,"x"))
g = vec(read(var,"g"))
timeProject = zero(Float64)
A = NLfunForward
opts = spgOptions(  optTol = 1e-5,
                    bpTol = 1e-5,
                    decTol = 1e-5,
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
f, x, r, nLine, stepG, lnErr, localProdA = spglinecurvy(A, x, gStep*g, maximum(lastFv),
                                A, funLS, b, project, timeProject, tau, opts, params)
