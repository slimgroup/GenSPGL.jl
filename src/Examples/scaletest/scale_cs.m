%% ML-JL SPGL Comparison
function [t]=scale_cs(n,k, i)
    mypath = '/home/Software/slim/';
    addpath(genpath([mypath 'SLIM-release-apps/tools/solvers/GenSPGL1/']));

    % Sparse Vector x0
    x0 = zeros(n,1);
    keep = floor(n/1000);
    x0(randperm(n, keep)) = randn(keep,1);
    r_inds = randperm(n, k);
    
    save('x0.mat', 'x0', 'r_inds', '-v7.3')
    % Load sparse vector and restriction mask
    n = length(x0);

    % Sparsifying Transform
    S  = opDCT(n);

    % Restriction Operator
    R  = opRestriction(n, r_inds);

    % Modelling Operator
    A = R*S;

    % Create Data
    b  = A*x0;

    % Solve
    parms = spgSetParms('verbosity', 1, 'iterations', 100);
    for ii = 1:i
        tic;
        [x,r,g,info] = spgl1(A, b,[], [], [], parms);
        t(ii) = toc;
    end
end
