%% ML-JL SPGL Comparison
clear

% Define Size of problem
n  = 5000;
k  = 200; 

% Create data from sparse vector
p  = randperm(n);
x0 = zeros(n,1);
x0(p(1:k)) = sign(randn(k,1));
A  = opMatrix(randn(n));
ind = randperm(n);
ind = ind(1:floor(0.6*length(ind)));
R  = opRestriction(n,ind);
b  = R*A*x0;

% Solve
opts = spgSetParms('optTol',1e-4);
tic
[x,r,g,info] = spgl1(R*A, b, 0, 1e-3, [], opts); % Find BP sol'n.]
toc
SNR = -20*log10(norm(x0-x)/norm(x0));
