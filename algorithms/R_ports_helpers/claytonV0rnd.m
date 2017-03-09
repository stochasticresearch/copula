function [y] = claytonV0rnd(n, theta)
%CLAYTONV0 - LS^{-1}[exp(-V_0psi_0^{-1}(psi_1(t)))]
%
% Inputs:
%  n - the number of V0's to generate
%  theta - dependency parameter
%
% Outputs:
%  y - a column vector of n random variables which follows the distribution 
%      of the inverse Laplace Stiles integral given a dependency parameter.

y = gamrnd(1.0/theta, 1, n, 1);

end