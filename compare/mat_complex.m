clear;
addpath(genpath('/home/Software/slim/SLIM-release-apps/tools/solvers/GenSPGL1/'));
load('vars_complex.mat')
sigma = 1e-2*norm(vec(b),'fro');
tic;
[xLS,r1,g1,info1] = spgl1(@NLfunForward,b(:),tau,sigma,xinit,opts,params);
toc

% Reconstruct MATLAB solution and calc SNR
e   = params.numr*params.nr;
L  = xLS(1:e);
R  = xLS(e+1:end);
L  = reshape(L,params.numr,params.nr);
R  = reshape(R,params.numc,params.nr);
Drecf = reshape(MH'*vec(L*R'),nr,nc);
SNR = -20*log10(norm(test-Drecf,'fro')/norm(test,'fro'));

% Reconstruct Julia solution and calc SNR
load('/home/slim/klensink/.julia/v0.6/GenSPGL/compare/xLS_jl.mat')
L_jl  = xLS_jl(1:e);
R_jl  = xLS_jl(e+1:end);
L_jl  = reshape(L_jl,params.numr,params.nr);
R_jl  = reshape(R_jl,params.numc,params.nr);
Drecf_jl = reshape(MH'*vec(L_jl*R_jl'),nr,nc);
SNR_jl = -20*log10(norm(test-Drecf_jl,'fro')/norm(test,'fro'));
