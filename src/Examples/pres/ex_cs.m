%% ML-JL SPGL Comparison
clear

mypath = '/home/Software/slim/';
addpath(genpath([mypath 'SLIM-release-apps/tools/solvers/GenSPGL1/']));

% Define Size of problem
%n  = 1024;
%k  = 50; 

% Sparse Vector x0
%x0 = zeros(n,1);
%x0(randperm(n, 10)) = randn(10,1);
%r_inds = randperm(n, k);

% Load sparse vector and restriction mask
load('x0.mat')
n = length(x0);

% Sparsifying Transform
S  = opDCT(n);

% Restriction Operator
R  = opRestriction(n, r_inds);

% Modelling Operator
A = R*S

% Create Data
b  = A*x0;

% Solve
tic
[x,r,g,info] = spgl1(A, b);
t_ml = toc;
SNR = -20*log10(norm(x0-x)/norm(x0));

%Compare
load('x_jl.mat')
figure(1)
subplot 121
    plot(1:length(x), x, 'k-', 'LineWidth', 2); hold on
    plot(1:length(x_jl), x_jl,'r--', 'LineWidth', 2)
    legend('MATLAB', 'JULIA')
    title(['Solution, Julia Speedup:' num2str(round(t_ml/t_jl, 2))])
subplot 122
    plot(x-x_jl)
    title('Residual')
