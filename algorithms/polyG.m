function [ y ] = polyG( lx, one_div_alpha, d, varargin )
%POLYG - Compute the coefficients a_{dk}(theta) involved in the generator derivatives
%and the copula density of a Gumbel copula
% Inputs:
%  lx = log(x); where x: evaluation point (vector);
%        e.g., for copGumbel@dacopula, lx = alpha*log(rowSums(iPsi(u)))
%        where u = (u_1,..,u_d) is the evaluation point of the density of Joe's copula)
%  one_div_alpha - alpha parameter (1/theta) in (0,1];
%      Please NOTE!! There is a difference in terminology between Mathworks
%      and R, in R (and seemingly in academic literature), theta is used as
%      the dependency parameter, and alpha = 1/theta.  Mathworks seems to
%      like to define the dependency parameter as alpha.  In this code, we
%      stick to the (unfortunate) Mathworks convention of using alpha as
%      the dependency parameter!  Therefore, in here we have one_div_alpha,
%      meaning that the input to this function should be 1/theta, where
%      theta is the dependency parameter of the Gumbel copula 
%  d number of summands, >= 1
% Outputs:
%  y - the output value
%
% Acknowledgements:
% This code modeled after the function definition in the paper:
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
        warning('Invalid varargin{1} for polyG, defaulting to LOG=FALSE');
        wantLog = 0;
    end
end

% TODO:
%stopifnot(length(alpha)==1, 0 < alpha, alpha <= 1,
%              d == as.integer(d), d >= 1)

M = size(lx,1);
x = exp(lx);

k = 1:d;
s = zeros(1,d);
for ss=1:d
    s(ss) = stirling1(d,ss);
end
S = cell(1,d);
for ss=1:d
    S_tmp = zeros(1,ss);
    for sss=1:ss
        S_tmp(sss) = stirling2(ss,sss);
    end
    S{ss} = S_tmp;
end


% slow :( TODO: vectorize this
y = zeros(size(lx));
for mm=1:M
    lst = zeros(1,length(k));
    for kk=k
        lst(kk) = (-1)^(d-1)*x(mm)*one_div_alpha^kk*s(kk)*polyval(fliplr(S{kk}),-x(mm));
    end
    zz = sum(lst);
    y(mm,:) = zz;
end

if(wantLog)
    y = log(y);
end

end

