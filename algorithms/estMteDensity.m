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
[f,xi] = ksdensity(x);
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

function [intervals, numSubsets] = genIntervals(f, xi)

% determine inflection points
[peaks, valleys] = findinflections(0,f);
allPts = sort([peaks(:); valleys(:)]');
numSubsets = length(allPts)+1;
allPtsAugment = [0 allPts length(xi)];

mte_params = cell(1,numSubsets);

% make a matrix of intervals over which to calculate the MTE parameters
intervals = zeros(numSubsets,2);
intervals(1,1) = allPtsAugment(1)+1; intervals(1,2) = allPtsAugment(2);    
minIntervalSize = 4;
for ii=2:numSubsets
    % if any intervals are less than a configurable amount of points,
    % extend the intervals so that curve fitting can work
    prevIntervalSize = intervals(ii-1,2)-intervals(ii-1,1);
    if(prevIntervalSize<minIntervalSize)
        intervals(ii-1,2) = intervals(ii-1,2) + (minIntervalSize-prevIntervalSize);
    end

    intervals(ii,1) = intervals(ii-1,2)+1;
    intervals(ii,2) = allPtsAugment(ii+1);
end
% if the last interval is less than the minimum interval size, merge it
% into the second to last interval
lastIntervalSize = intervals(end,2)-intervals(end,1);
if(lastIntervalSize<minIntervalSize)
    intervals(end-1,2) = length(xi);
    intervals(end,:) = [];
    numSubsets = numSubsets - 1;
end

end