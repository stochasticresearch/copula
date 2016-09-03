%**************************************************************************
%*                                                                        *
%* Copyright (C) 2016  Kiran Karra <kiran.karra@gmail.com>                *
%*                                                                        *
%* This program is free software: you can redistribute it and/or modify   *
%* it under the terms of the GNU General Public License as published by   *
%* the Free Software Foundation, either version 3 of the License, or      *
%* (at your option) any later version.                                    *
%*                                                                        *
%* This program is distributed in the hope that it will be useful,        *
%* but WITHOUT ANY WARRANTY; without even the implied warranty of         *
%* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          *
%* GNU General Public License for more details.                           *
%*                                                                        *
%* You should have received a copy of the GNU General Public License      *
%* along with this program.  If not, see <http://www.gnu.org/licenses/>.  *
%*                                                                        *
%**************************************************************************

% A simple test to see how discretization affects likelihood

clear;
clc;

rng(12345);

M = 1000;
X = randn(M*2,1);
X_train = X(1:M);
X_test = X(M+1:end);

[muhat, sigmahat] = normfit(X);
nlogL = normlike([muhat,sigmahat],X_test);

% discretize the RV and compute the likelihood
[X_train_discretized, edges] = discretizeRv(X_train, 10);
% [testDiscretizedIdxs] = discretize(X_test, edges);
[~,~,testDiscretizedIdxs] = histcounts(X_test, edges);
% deal w/ NaN's :(
nanIdxs = [find(isnan(testDiscretizedIdxs)) find(testDiscretizedIdxs==0)];
belowMinIdxs = X_test(nanIdxs)<min(edges);
aboveMaxIdxs = X_test(nanIdxs)>max(edges);
testDiscretizedIdxs(nanIdxs(belowMinIdxs)) = 1;
testDiscretizedIdxs(nanIdxs(aboveMaxIdxs)) = edges(end);
testDiscretizedIdxs = floor(testDiscretizedIdxs);
X_test_discretized = edges(testDiscretizedIdxs);

% make X_test_discretized into a rvEmpiricalInfo object
isdiscrete = 1;
[f,xi] = emppdf(X_test_discretized', isdiscrete);
F = empcdf(X_test_discretized, isdiscrete);
discreteEmpInfo = rvEmpiricalInfo(xi, f, F, isdiscrete);

% compute negative log-likelihood manually for continuous and discretized
% versions
nLLVal_continuous = 0;
nLLVal_discrete = 0;
LOG_CUTOFF = 1e-5;
for mm=1:M
    continuousProb = normpdf(X_test(mm), muhat, sigmahat);
    if(continuousProb<LOG_CUTOFF)
        continuousProb = LOG_CUTOFF;
    end
    
    discreteProb = discreteEmpInfo.pdf(X_test_discretized(mm));
    if(discreteProb<LOG_CUTOFF)
        discreteProb = LOG_CUTOFF;
    end
    
    fprintf('X=%f \t X_disc=%f \t continuousProb=%f \t discreteProb=%f\n', ...
        X_test(mm), X_test_discretized(mm), continuousProb, discreteProb);
    
    nLLVal_continuous = nLLVal_continuous + log(normpdf(X_test(mm), muhat, sigmahat));
    nLLVal_discrete = nLLVal_discrete + log(discreteProb);
end

fprintf('Matlab=%f Continuous=%f Discretized=%f\n', -1*nlogL, nLLVal_continuous, nLLVal_discrete);

%% See the shape discretization produces

clear;
clc;

M = 1000;
X = randn(M*2,1);
[X_discretized, edges] = discretizeRv(X, 10);