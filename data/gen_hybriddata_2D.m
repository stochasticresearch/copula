function [ X, U ] = gen_hybriddata_2D( n )
%GEN_HYBRIDDATA_2D Generates 2-D mixed joint distribution
% Inputs:
%  n - the number of RV's to generate
%
% Outputs:
%  y - the random variables, has size = [n x 2]

rho = .7;
Z = mvnrnd([0 0], [1 rho; rho 1], n);
U = normcdf(Z);

alpha = 5;
beta = 1;
betaDist = makedist('Beta', alpha, beta);
lambda = 0.5;
poissonDist = makedist('Poisson', lambda);

X = [icdf(poissonDist, U(:,1)) icdf(betaDist, U(:,2))];

end

