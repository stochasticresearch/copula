%**************************************************************************
%* 
%* Copyright (C) 2017  Kiran Karra <kiran.karra@gmail.com>
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

%% Test 1 - entire multivariate distribution
% Generate a multivariate normal distribution w/ Gaussian Copula & Normal
% marginals, then ensure that the probability matches what is reported by
% copfamprob and by mvnpdf

clear;
clc;

D = 3;
R = [1 0 .2; 0 1 -.8; .2 -.8 1];
MVec = 100:100:1000;

% define the corresponding copula
copModel = 'gaussian';
copParams = R;

for M=MVec
    U = copularnd(copModel, copParams, M);
    data = U;

    isdiscrete = 0;

    % define the marginals -- we define them as rvEmpiricalInfo objects to
    % conform to the API
    rvEmpInfoObjs = cell(1,D);
    for ii=1:D
        data(:,ii) = norminv(U(:,ii),0,1);

        [f, xi] = emppdf(data(:,ii), isdiscrete);
        F = empcdf(data(:,ii), isdiscrete);
        obj = rvEmpiricalInfo(xi, f, F, isdiscrete);
        rvEmpInfoObjs{ii} = obj;
    end

    K = 50;
    idxsToTest = randperm(M); idxsToTest = idxsToTest(1:K);

    mse = 0;
    for ii=1:K
        p1 = copfamprob(data(ii,:), rvEmpInfoObjs, copModel, copParams);
        p2 = mvnpdf(data(ii,:),[0,0,0],R);
        mse = mse + (p1-p2).^2;
    end
    ave_mse = mse/K;
    fprintf('M=%d ave_mse=%0.05f\n', M, ave_mse);
end

%% Test 2 - marginals of the multivariate distribution