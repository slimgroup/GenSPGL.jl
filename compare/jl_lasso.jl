# Load a common workspace and solve a lasso problem
using MAT
vars = matread("vars_lasso.mat")

A = vars["A"]
b = vec(vars["b"])
tau = vars["tau"]

using GenSPGL
opts = spgOptions(verbosity = 1)

println("""
================================================================================
spg_lasso\n
""")
x, r, g, info = spg_lasso(A,b; tau = tau, options = opts)

println("""
================================================================================\n
Basis Pursuit\n
""")
sigma = 0.
tau = 0.

x, r, g, info = spgl1(A, b; tau = tau, sigma = sigma, options = opts)

