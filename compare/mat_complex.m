clear;
addpath(genpath('/home/Software/slim/SLIM-release-apps/tools/solvers/GenSPGL1/'));
load('vars_complex.mat')
sigma = 1e-2*norm(vec(b),'fro');
tic;
[xLS,r1,g1,info1] = spgl1(@NLfunForward,b(:),tau,sigma,xinit,opts,params);
toc
