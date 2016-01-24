function [ y ] = claytoncopulapdf( u, alpha, varargin )
%GUMBELCOPULAPDF Computes the PDF of the Clayton Copula for N>=2
% Inputs:
%  u - an M x N matrix of all the points over which to compute the Clayton
%      copula PDF, M is the number of points in a unit-hypercube of
%      dimension N
%  alpha - the dependency parameter of the Clayton copula.  
%      Please NOTE!! There is a difference in terminology between Mathworks
%      and R, in R (and seemingly in academic literature), theta is used as
%      the dependency parameter, and alpha = 1/theta.  Mathworks seems to
%      like to define the dependency parameter as alpha.  In this code, we
%      stick to the (unfortunate) Mathworks convention of using alpha as
%      the dependency parameter!!, 
%      Consequently we assign 1_div_alpha = 1/alpha
%
% Outputs:
%  y - the value of the Clayton copula density at the specified value in
%      the unit hypercube
%
% Acknowledgements:
% This code modeled after the function definitions in the paper:
% Estimators for Archimedean copulas in high dimensions, 
% by Marius Hofert1, Martin Mächler, and Alexander J. McNeil AND
% the R code in cop_objects.R in the copula package in R.

nVarargs = length(varargin);
if(nVarargs==0)
    wantLog = 0;
else
    if(isnumeric(varargin{1}))
        wantLog = varargin{1};
    else
        warning('Invalid varargin{1} for claytoncopulapdf, defaulting to LOG=FALSE');
        wantLog = 0;
    end
end


[~,N] = size(u);
if(N==2)
    y = copulapdf('Clayton', u, alpha);
else
    % compute log of the pdf
    lu = sum(log(u),2);
    t = sum(iPsi_clayton(u, alpha), 2);
    y = sum(log1p(alpha*(0:N-1))) - (1+alpha)*lu - (N + 1.0/alpha)*log1p(t);
    % exponentiate the result if the desired value is not the log version
    if(~wantLog)
        y = exp(y);
    end
end

end

function [ y ] = iPsi_clayton( u, alpha, varargin )

nVarargs = length(varargin);
if(nVarargs==0)
    useLog = 0;
else
    if(isnumeric(varargin{1}))
        useLog = varargin{1};
    else
        warning('Invalid varargin{1} for iPsi_clayton, defaulting to LOG=FALSE');
        useLog = 0;
    end
end

y = u.^(-alpha) - 1;
if(useLog)
    y = log(y);
end

end