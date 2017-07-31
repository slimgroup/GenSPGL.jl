## JL-ML SPGL Comparison

# Create Restriction Operator
using JOLI
n = 512
k = 20

p = randperm(n)
x0 = zeros(n)
x0[p[1:k]] = sign.(randn(k))

A = joMatrix(randn(n,n))
ind = randperm(n)
ind = ind[1:Int(floor(0.6*n))]
R = joRestriction(n,ind)

# Create data
b = R*A*x0

# Solve
using GenSPGL
opts = spgOptions(optTol = 1e-4,
                     verbosity = 1)

@time x, r, g, info = spgl1(R*A, vec(b), tau = 0., sigma = 1e-3, options = opts) 

using PyPlot
plot(1:n,x0, 1:n,x)
SNR = -20*log10(norm(x0-x)/norm(x0));
title("GenSPGL Solution - SNR: $(SNR)")

