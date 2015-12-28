function [mte_info] = estMteDensity(X)
%ESTMTEDENSITYUNIVARIATE - estimates the density of a univariate dataset 
% with a mixture of truncated exponentials
% Inputs
%  x - a MxD matrix of data-points from the source density (unknown)
% Outputs:
%  mte_params - a rvEmpiricalInfo object w/ the domain and the density
%
% Code inspired by the paper: "Estimating Mixtures of Truncated
% Exponentials in Hybrid Bayesian Networks" - by Rumi, Salmeron, and Moral
% Requires Matlab's Curve Fitting Toolbox

D = size(X,2);
mte_params_1 = estMteDensityUnivariate(X(:,1));
if(D==1)
    % reconstruct the pdf from the MTE
    mte_pdf_estimate = [];
    xi = [];
    for ii=1:length(mte_params_1)
        x = mte_params_1{ii}.xi_subset;
        a = mte_params_1{ii}.a; b = mte_params_1{ii}.b;
        c = mte_params_1{ii}.c; d = mte_params_1{ii}.d;
        mte_subset_reconstruct = a*exp(b*x)+c*exp(d*x);
        mte_pdf_estimate = [mte_pdf_estimate; mte_subset_reconstruct];
        xi = [xi; x];
    end
    mte_info = rvEmpiricalInfo(xi, mte_pdf_estimate, []);
else
    error('Currently unsupported!\n');
end

end

function [mte_params] = estMteDensityUnivariate(x)
%ESTMTEDENSITYUNIVARIATE - estimates the density of a univariate dataset 
% with a mixture of truncated exponentials

% perform KDE on the given dataset
[f,xi] = ksdensity(x,'width',10);       % high bandwidth to avoid too many sub-intervals
f = f(:); xi = xi(:);

[intervals, numSubsets] = genIntervals(f,xi);

% for each subset, estimate MTE parameters of the form:
%  y = a*exp(b*x) + h*c*exp(d*x) + k
for subset=1:numSubsets
    f_subset = f(intervals(subset,1):intervals(subset,2));
    xi_subset = xi(intervals(subset,1):intervals(subset,2));
    
    % fit to 2 term exponential using matlab's fit functionality
    f2 = fit(xi_subset,f_subset,'exp2');
    
    tmp = struct;
    tmp.subset = subset;
    tmp.f_subset = f_subset;
    tmp.xi_subset = xi_subset;
    tmp.interval = [intervals(subset,1) intervals(subset,2)];
    tmp.a = f2.a;
    tmp.b = f2.b;
    tmp.c = f2.c;
    tmp.d = f2.d;
    
    mte_params{subset} = tmp;
end

end