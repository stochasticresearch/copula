function [ y ] = frankcopulapdf( u, alpha, varargin )
%GUMBELCOPULAPDF Computes the PDF of the Frank Copula for N>=2
% Inputs:
%  u - an M x N matrix of all the points over which to compute the Frank
%      copula PDF, M is the number of points in a unit-hypercube of
%      dimension N
%  alpha - the dependency parameter of the Frank copula.  
%      Please NOTE!! There is a difference in terminology between Mathworks
%      and R, in R (and seemingly in academic literature), theta is used as
%      the dependency parameter, and alpha = 1/theta.  Mathworks seems to
%      like to define the dependency parameter as alpha.  In this code, we
%      stick to the (unfortunate) Mathworks convention of using alpha as
%      the dependency parameter!!, 
%      Consequently we assign 1_div_alpha = 1/alpha
% Optional Inputs:
%  varargin{1} - if 0, then compute pdf directly, 
%                else, compute log of pdf at specified value
%
% Outputs:
%  y - the value of the Clayton copula density at the specified value in
%      the unit hypercube
%
% Acknowledgements:
% This code modeled after the function definitions in the paper:
% Estimators for Archimedean copulas in high dimensions, 
% by Marius Hofert1, Martin Mächler, and Alexander J. McNeil AND
% the R code in cop_objects.R in the "copula" package in R, found at:
% https://cran.r-project.org/web/packages/copula/

nVarargs = length(varargin);
if(nVarargs==0)
    wantLog = 0;
else
    if(isnumeric(varargin{1}))
        wantLog = varargin{1};
    else
        warning('Invalid varargin{1} for frankcopulapdf, defaulting to LOG=FALSE');
        wantLog = 0;
    end
end

[~,N] = size(u);

u_sum = sum(u,2);
lp  = log1mexp(alpha);    % log(1 - exp(-alpha))
lpu = log1mexp(alpha*u); % log(1 - exp(-alpha * u))
lu  = sum(lpu,2);

Li_arg = lp + sum(lpu-lp,2);
Li = log(polylog(-(N-1), exp(Li_arg)));
y = (N-1)*log(alpha) + Li - alpha*u_sum - lu;

if(~wantLog)
    y = exp(y);
end

end