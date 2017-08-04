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
    mytau=convert(Float64,pi)

    # Optionally set the verbosity to level 1
    opts = spgOptions(verbosity = 1)

    x, r, g, info = spgl1(A, vec(b); tau = mytau, options = opts)

    # Calc the SNR of the recovery
    SNR = -20*log10(norm(x0-x)/norm(x0));

    return x, r, g, info, SNR

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
    mytau=convert(Float64,pi)

    # LASSO
    sigma = []
    x0=[]

    # Optionally set the verbosity to level 1
    opts = spgOptions(verbosity = 1)

    A_im = joMatrix(A; name="A")
    x, r, g, info = spgl1(A_im, vec(b); tau = mytau, options = opts)

    # Calc the SNR of the recovery
    SNR = -20*log10(norm(x0-x)/norm(x0));

    return x, r, g, info, SNR
end
