# Load a common workspace and solve a lasso problem
using MAT
vars = matread("vars_lasso.mat")

A = vars["A"]
b = vec(vars["b"])
b_noise = vec(vars["b_noise"])
tau = vars["tau"]
x0 = vec(vars["x0"])

using GenSPGL
opts = spgOptions(verbosity = 1)

println("""
================================================================================
spg_lasso\n
""")
x_lasso, r, g, info = spg_lasso(A,b; tau = tau, options = opts)

println("""
================================================================================\n
Basis Pursuit\n
""")
sigma = 0.
tau = 0.

x_bp, r, g, info = spgl1(A, b; tau = tau, sigma = sigma, options = opts)

println("""
================================================================================\n
Basis Pursuit Denoise\n
""")
sigma = 0.10
tau = 0.

x_bpdn, r, g, info = spgl1(A, b_noise; tau = tau, sigma = sigma, options = opts)



