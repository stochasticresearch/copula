function [y] = gumbelV0rnd(n, theta)
%GUMBELV0 - LS^{-1}[exp(-t)]
% Sample from S(alpha,1,(cos(alpha*pi/2))^(1/alpha),0;1)
% with Laplace-Stieltjes transform exp(-t^alpha)
%
% Inputs:
%  n - the number of V0's to generate
%  theta - dependency parameter
%
% Outputs:
%  y - a column vector of n random variables which follows the distribution 
%      of the inverse Laplace Stiles integral given a dependency parameter.

if(theta==1)
    y = ones(n,1);
else
    alpha = 1/theta;
    gamma = cos(pi/2.0*alpha)^(1/alpha);
    delta = 0;
    y = stable1rnd(n, alpha, gamma, delta);
end

end