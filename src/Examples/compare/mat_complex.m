clear;

% Prepend this path appropriately 
mypath = '/home/Software/slim/';
addpath(genpath([mypath 'SLIM-release-apps/tools/solvers/GenSPGL1/']));
load('vars_complex.mat')
sigma = 1e-2*norm(vec(b),'fro');
tic;
[xLS,r1,g1,info1] = spgl1(@NLfunForward,b(:),tau,sigma,xinit,opts,params);
toc

% Recreate MH
perc=0.2;
inds =randperm(nr);
ind =inds(1:floor(nc*0.2));
R1  = opRestriction(nc,ind);
R   = opKron(R1,opDirac(nr));
MH  = opMH(nr,nc);

% Reconstruct MATLAB solution and calc SNR
e   = params.numr*params.nr;
L  = xLS(1:e);
R  = xLS(e+1:end);
L  = reshape(L,params.numr,params.nr);
R  = reshape(R,params.numc,params.nr);
Drecf = reshape(MH'*vec(L*R'),nr,nc);
SNR = -20*log10(norm(test-Drecf,'fro')/norm(test,'fro'));

% Reconstruct Julia solution and calc SNR
load('xLS_jl.mat')
L_jl  = xLS_jl(1:e);
R_jl  = xLS_jl(e+1:end);
L_jl  = reshape(L_jl,params.numr,params.nr);
R_jl  = reshape(R_jl,params.numc,params.nr);
Drecf_jl = reshape(MH'*vec(L_jl*R_jl'),nr,nc);
SNR_jl = -20*log10(norm(test-Drecf_jl,'fro')/norm(test,'fro'));
