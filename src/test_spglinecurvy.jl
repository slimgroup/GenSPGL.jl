# Script to test spglinecurvy

include("spglinecurvy.jl")

using GenSPGL
using MAT
var = matread("/scratch/slim/klensink/test/spglinecurvy_vars.mat")

gStep= var["gStep"]
tau = var["tau"]
lastFv = vec(var["lastFv"])
b = vec(var["b"])
x = vec(var["x"])
g = vec(var["g"])
params = Dict{String, Int64}()
timeProject = zero(Float64)
options = spgOptions()

# Create explicit A
m = 50; n = 128; k = 14
A,Rtmp = qr(randn(n,m))
A = A'

fNew, xNew, rNew, iter, step, err, nProd = spglinecurvy(A, x, gStep*g, maximum(lastFv), SpotFunForward, funLS, b, project, timeProject, tau, options, params)
