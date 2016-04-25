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

% A Script to test Log-Likelihood calculation

%% test multi-modal MTE vs Gaussian
clear;
clc;

M = 500;
x = [normrnd(-2,0.3,M,1); normrnd(2,0.5,M,1)];
x = x(randperm(M*2));

x_train = x(1:M);
x_test = x(M+1:end);

% estimate w gaussian
[Mean, Covariance] = ecmnmle(x_train);

% estimate MTE
mte_info = estMteDensity(x_train);

% calculate likelihoods
ll_clg = 0;
ll_mte = 0;
minProbVal = 10e-5;
for ii=1:size(x_test,1)
    clg_likelihood = normpdf(x_test(ii), Mean, Covariance);
    if(clg_likelihood<minProbVal)
        clg_likelihood = minProbVal;
    end
    
    mte_likelihood = mte_info.pdf(x_test(ii));
    if(mte_likelihood<minProbVal)
        mte_likelihood = minProbVal;
    end
    
    ll_clg = ll_clg + log(clg_likelihood);
    ll_mte = ll_mte + log(mte_likelihood);
end

[zz,xzz]=ksdensity(x_train);

plot(xzz,zz,mte_info.domain, normpdf(mte_info.domain, Mean, Covariance), mte_info.domain, mte_info.density);
grid on; legend('True','Gaussian Approx', 'MTE Approx');
title(sprintf('CLG=%0.02f MTE=%0.02f',ll_clg,ll_mte));