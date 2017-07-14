# Script to test spglinecurvy

using GenSPGL
using MAT
fid = matopen("../compare/spgLine_vars.mat")

f = read(fid, "f")
dx = vec(read(fid, "dx"))
gtd = read(fid, "gtd")
lastFv = vec(read(fid, "lastFv"))
funForward = NLfunForward
funPenalty = funLS
b = vec(read(fid, "b"))
x = vec(read(fid, "x"))
feasSrchIt = convert(Int64, read(fid,"feasSrchIt"))
timeProject = zero(Float64)
linear_in = read(fid,"linear")

if linear_in == 0
    linear = false
else
    linear = true
end
A = funForward
opts = spgOptions(  optTol = 1e-5,
                    bpTol = 1e-5,
                    decTol = 1e-5,
                    project = TraceNorm_project,
                    primal_norm = TraceNorm_primal,
                    dual_norm = TraceNorm_dual,
                    proxy = true,
                    ignorePErr = true,
                    iterations = 150,
                    verbosity = 2)

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

fNew, step, rNew, iter, err, localProdA = spgline(A,
        f,
        dx,
        gtd,
        x,
        maximum(lastFv),
        funForward,
        funPenalty,
        params,
        b,
        feasSrchIt,
        linear,
        opts,
        timeProject)
