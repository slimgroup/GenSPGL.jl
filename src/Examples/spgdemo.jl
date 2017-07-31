export spgdemo_ex, spgdemo_im
using JOLI

function spgdemo_ex()

# Create random m-by-n encoding matric and sparse vector
m = 50; n = 128; k = 14
A,Rtmp = qr(randn(n,m))
A = A'
p = randperm(n)
p = p[1:k]
x0 = zeros(n,1)
x0[p] = randn(k,1);

b = A*x0
tau=convert(Float64,pi)

# LASSO
sigma = []
x0=[]
opts = spgOptions(verbosity = 1)

out = spg_lasso(A,vec(b); tau = tau, options = opts)

end

function spgdemo_im()

# Create random m-by-n encoding matric and sparse vector
m = 50; n = 128; k = 14
A,Rtmp = qr(randn(n,m))
A = A'
p = randperm(n)
p = p[1:k]
x0 = zeros(n,1)
x0[p] = randn(k,1);

b = A*x0
tau=convert(Float64,pi)

# LASSO
sigma = []
x0=[]
opts = spgOptions(verbosity = 1)

A_im = joMatrix(A; name="A")
out = spg_lasso(A_im, vec(b); tau = tau, options = opts)

end
