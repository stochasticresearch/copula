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

% Test the maxcopfamprob optimization routine

%% Start the test w/ a Gaussian
clear;
clc;

M = 5;
D = 3;
R = [1 0 .2; 0 1 -.8; .2 -.8 1];

% generate data from a multivariate normal that has this Rho matrix, and
% normal marginals w/ mu=0,std=1
data = mvnrnd(zeros(1,D),R,M);

% define the corresponding copula
copModel = 'gaussian';
copParams = R;

% define the marginals -- we define them as rvEmpiricalInfo objects to
% conform to the API
rvEmpInfoObjs = cell(1,D);
for ii=1:D
    tmpData = randn(1000,1);
    [f, xi] = emppdf(tmpData, 0);
    F = empcdf(tmpData, 0);
    isdiscrete = 0;
    obj = rvEmpiricalInfo(xi, f, F, isdiscrete);
    rvEmpInfoObjs{ii} = obj;
end

for ii=1:M
    % randomely determine which of the 3 to drop.  For now we are only
    % dropping 1
    missingIdxs = randi(D);
    xIn = data(ii,:);
    % for now, we test only continuous optimization
    intIdxs = [];
    
    optimOut = maxcopfamprob(copModel,copParams,xIn,rvEmpInfoObjs,missingIdxs,intIdxs);
    
end