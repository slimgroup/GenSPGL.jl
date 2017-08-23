using PyPlot, JLD, MAT, opesciSLIM.TimeModeling
include("/home/slim/klensink/GIT/WG_3D_Coil/julia/NLfunForward_local.jl")

function data_load()
    cd("/scratch/slim/rkumar/CoilSyn")
    nsrcx  = 102;
    nsrcy  = 102;
    nrecx  = 202;
    nrecy  = 202;
    nt     = 751;
    nsrc   = 102^2;
    xsrc   = 5000.f0:100.f0:15100.f0;  # 500 m space from boundary
    ysrc   = 5000.f0:100.f0:15100.f0;

    #= loop over all sources to extract any frequency slice
    # initialize frequency slice
    Df     = zeros(Complex{Float64}, nrecx,nrecy,nsrcx,nsrcy)
    index  = 80;

    tic();
    for i in 1:2# length(xsrc)
        for j in 1:2#length(ysrc)
            D_1 = load("shot_$(xsrc[i])_$(ysrc[j]).jld");
            D_2 = D_1["d"];
            D_3 = D_2.data[1];
            D_tmp = reshape(fft(D_3),nt,nrecx,nrecy);
            test = D_tmp[index,:,:];
            Df[:,:,i,j] = test;
        end
    end
    toc()
    =#
    
    Df = load("/scratch/slim/shared/rajiv_keegan/slice_80.jld")["slice"]

    # load matlab Mask
    file = matopen("/scratch/slim/shared/rajiv_keegan/MASK_Coil.mat");
    Mask = BitArray(read(file,"Mask"));
    close(file)

    Df = reshape(permutedims(Df,[1, 3, 2, 4]), 202*102, 202*102)
    D_sub = copy(Df)
    D_sub[.~Mask] = zero(eltype(D_sub))
    
    return D_sub, Mask
end

function scale_test(D_sub, Mask)
    # Choose options for GenSPGL
    opts = GenSPGL.spgOptions(  optTol = 1e-5,
                        bpTol = 1e-5,
                        decTol = 1e-4,
                        project = GenSPGL.TraceNorm_project,
                        primal_norm = GenSPGL.TraceNorm_primal,
                        dual_norm = GenSPGL.TraceNorm_dual,
                        proxy = true,
                        ignorePErr = true,
                        iterations = 50,
                        verbosity = 1)

    
    nr, nc = size(D_sub)
    # Create Params Dict
    params = Dict("nr"=> 100,
                    "numr"=> nr,
                    "numc"=> nc,
                    "funForward"=> NLfunForward_local,
                    "mode"=> 1,
                    "ls"=> 1,
                    "mask" => Mask,
                    "logical"=> 0,
                    "funPenalty"=> GenSPGL.funLS)

    # Set threads and solve problems
    Linit = complex(rand(params["numr"], params["nr"]))
    Rinit = complex(rand(params["numr"], params["nr"]))
    xinit = [vec(Linit); vec(Rinit)]

    BLAS.set_num_threads(8)

    t = GenSPGL.spgl1(NLfunForward_local, vec(D_sub), x = vec(xinit),
                                                           options = opts,
                                                           params = params,
                                                           tau = 3.14e1)
    return t
end
