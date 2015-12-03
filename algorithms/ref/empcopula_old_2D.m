function [ C, U, c ] = empcopula_old_2D( X, n )
%EMPCOPULA Calculates the empiricial copula in a unit hypercube
% Inputs:
%  X - a [M x D] matrix, where M is the number of samples, and D is the
%      dimensionality of the data
%  n - the number of evenly spaced points in the unit hypercube to
%      calculate the empirical copula over
%
% Outputs:
%  C - the empirical copula, which is a [n x D] matrix
%  U - the points over which the copula was calculated
%  c - the empirical copula density, same dimensions as C

R = tiedrank(X);
M = size(X,1);
D = size(X,2);

if(D<2)
    error('D must be >= 2!');
end

bounds = linspace(0,1,n);
bounds = bounds';

C = zeros(length(bounds),length(bounds));
U = zeros(length(bounds),length(bounds),D);

% scale the ranks to create pseudoobservations
u_i = R(:,1)/(M+1); v_i = R(:,2)/(M+1);

% TODO: improve this to become n-dimensional, for now we handle 2
for ii=2:length(bounds)     % start from 2, b/c of copula grounded property
    for jj=2:length(bounds)
        u = bounds(ii);
        v = bounds(jj);
                
        total = sum( (u_i<=u).*(v_i<=v) );
        C(ii,jj) = total/M;
        
        U(ii,jj,1) = u;
        U(ii,jj,2) = v;
    end
end

% calculate c from C, just take the numerical derivative
[cc] = gradient(C);
[~,c] = gradient(cc);

end
