% Load a common workspace and solve a lasso problem
clear;

addpath(genpath('/home/Software/slim/SLIM-release-apps/tools/solvers/GenSPGL1/'));
load('vars_lasso.mat')

tic;
x_lasso = spg_lasso(A, b, tau);
toc;

display('------------- Basis Pursuit --------------')

opts = spgSetParms('verbosity',2);
x_bp = spg_bp(A, b, opts);
