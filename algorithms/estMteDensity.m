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
%
%**************************************************************************
%* 
%* Copyright (C) 2016  Kiran Karra <kiran.karra@gmail.com>
%*
%* This program is free software: you can redistribute it and/or modify
%* it under the terms of the GNU General Public License as published by
%* the Free Software Foundation, either version 3 of the License, or
%* (at your option) any later version.
%*
%* This program is distributed in the hope that it will be useful,
%* but WITHOUT ANY WARRANTY; without even the implied warranty of
%* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%* GNU General Public License for more details.
%*
%* You should have received a copy of the GNU General Public License
%* along with this program.  If not, see <http://www.gnu.org/licenses/>.
%* 
%**************************************************************************

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
        term1_exp = exp(b*x);
        term2_exp = exp(d*x);
        if(any(isinf(term1_exp)))
            if(a==0)
                term1_exp = zeros(length(x),1);
            else
                term1_exp = 10000*ones(length(x),1);
            end
        end
        if(any(isinf(term2_exp)))
            if(b==0)
                term2_exp = zeros(length(x),1);
            else
                term2_exp = 10000*ones(length(x),1);
            end
        end
        mte_subset_reconstruct = a*term1_exp+c*term2_exp;
        mte_pdf_estimate = [mte_pdf_estimate; mte_subset_reconstruct];
        xi = [xi; x];
    end
    % eliminate any negative probabilities from curve-fitting
    if(min(mte_pdf_estimate)<0)
        mte_pdf_estimate = mte_pdf_estimate + abs(min(mte_pdf_estimate));
    end
    % normalize to integrate to 1
    isdiscrete = 0;
    mte_pdf_estimate = mte_pdf_estimate/trapz(xi,mte_pdf_estimate);
    mte_info = rvEmpiricalInfo(xi, mte_pdf_estimate, [], isdiscrete);
else
    error('Currently unsupported!\n');
end

end

function [mte_params] = estMteDensityUnivariate(x)
%ESTMTEDENSITYUNIVARIATE - estimates the density of a univariate dataset 
% with a mixture of truncated exponentials

% perform KDE on the given dataset
% [f,xi] = ksdensity(x,'width',10);       % high bandwidth to avoid too many sub-intervals
[f,xi] = ksdensity(x);
f = f(:); xi = xi(:);

[intervals, numSubsets] = genIntervals(f,xi);

% for each subset, estimate MTE parameters of the form:
%  y = a*exp(b*x) + h*c*exp(d*x) + k
for subset=1:numSubsets
    f_subset = f(intervals(subset,1):intervals(subset,2));
    xi_subset = xi(intervals(subset,1):intervals(subset,2));
    
    % fit to 2 term exponential using matlab's fit functionality
    warning('error', 'curvefit:fit:invalidStartPoint');
    try
        f2 = fit(xi_subset,f_subset,'exp2');
    catch
        f2 = struct;
        f2.a = 1;
        f2.b = 1;
        f2.c = 1;
        f2.d = 1;
    end
    
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