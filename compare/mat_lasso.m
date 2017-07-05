% Load a common workspace and solve a lasso problem
addpath(genpath('/home/Software/slim/SLIM-release-apps/tools/solvers/GenSPGL1/'));
load('vars_lasso.mat')

spg_lasso(A, b, tau);
