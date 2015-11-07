function [ X, U ] = gen_hybriddata_3D( n )
%GEN_HYBRIDDATA_3D Generates 2-D mixed joint distribution
% Inputs:
%  n - the number of RV's to generate
%
% Outputs:
%  y - the random variables, has size = [n x 2]

Rho = [1 .4 .2; .4 1 -.8; .2 -.8 1];
Z = mvnrnd([0 0 0], Rho, n);
U = normcdf(Z,0,1);

alpha = 5;
beta = 1;
betaDist = makedist('Beta', alpha, beta);
lambda = 0.5;
poissonDist = makedist('Poisson', lambda);

X = [icdf(poissonDist, U(:,1)) icdf(betaDist, U(:,2)), tinv(U(:,3),5)];

end

