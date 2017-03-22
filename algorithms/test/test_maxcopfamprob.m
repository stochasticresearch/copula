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

M = 1000;
D = 3;
R = [1 0 .2; 0 1 -.8; .2 -.8 1];

% define the corresponding copula
copModel = 'gaussian';
copParams = R;
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

% Test 1-D Continuous Optimization
f = figure;
for ii=1:5
    % randomely determine which of the 3 to drop.  For now we are only
    % dropping 1
    missingIdxs = randi(D);
    xIn = data(ii,:);
    intIdxs = [];
    
    [optimOut,fval] = maxcopfamprob(copModel,copParams,xIn,rvEmpInfoObjs,missingIdxs,intIdxs);
    
    % sanity check to see if the optimizer is working correctly
    lb = min(rvEmpInfoObjs{missingIdxs}.domain);
    ub = max(rvEmpInfoObjs{missingIdxs}.domain);
    valVec = linspace(lb,ub,100);
    Z = zeros(1,100);
    for jj=1:100
        val = valVec(jj);
        xx = xIn;
        xx(missingIdxs) = val;
        Z(jj) = copfamprob(xx,rvEmpInfoObjs,copModel,copParams);
    end
    
    clf(f);
    plot(valVec,Z); grid on; hold on;
    plot(optimOut,-1*fval,'r*','MarkerSize',14);
    pause;
end

% Test 2-D Continuous Optimization
for ii=1:5
    idxs = randperm(3);
    missingIdxs = sort(idxs(1:2));
    xIn = data(ii,:);
    intIdxs = [];
    
    [optimOut,fval] = maxcopfamprob(copModel,copParams,xIn,rvEmpInfoObjs,missingIdxs,intIdxs);
    
    % sanity check to see if the optimizer is working correctly
    lb = zeros(1,length(missingIdxs)); ub = zeros(1,length(missingIdxs));
    for jj=1:length(missingIdxs)
        lb(jj) = min(rvEmpInfoObjs{missingIdxs(jj)}.domain);
        ub(jj) = max(rvEmpInfoObjs{missingIdxs(jj)}.domain);
    end
    
    gridX = linspace(lb(1),ub(1),100);
    gridY = linspace(lb(2),ub(2),100);
    [X,Y] = ndgrid(gridX,gridY);
    Z = zeros(length(X),length(Y));
    for jj=1:length(X)
        for kk=1:length(Y)
            xx = xIn;
            xx(missingIdxs(1)) = X(jj);
            xx(missingIdxs(2)) = Y(kk);
            Z(jj,kk) = copfamprob(xx,rvEmpInfoObjs,copModel,copParams);
        end
    end
    
    clf(f);
    surf(X,Y,Z); grid on; hold on;
    plot(optimOut,-1*fval,'r*','MarkerSize',14);
    pause;
end
close(f);