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

%% another test

clear;
clc;
continuousType = 'other';
discreteType = {};
nodeA = [0.4 0.3 0.2 0.1]; discreteType{1} = nodeA;

a_dist = makedist('Multinomial','Probabilities',discreteType{1});

M = 2000;
% generate the copula random variables
D = 2; alpha = 10;
U = claytoncopularnd(M, D, alpha);
% U = copularnd('Gaussian', [1 0.8; 0.8 1], M);
X = zeros(M,2);

X(:,1) = a_dist.icdf(U(:,1));
% perform the inverse transform to generate the X data matrix
if(strcmpi(continuousType,'Gaussian'))
    rhoD = 0.6;
    X(:,2) = norminv(U(:,2),0,rhoD);
else
    % make it bimodal
%     x = [normrnd(-2,0.3,M/2,1); normrnd(2,0.5,M/2,1)];
%     [f,xi] = emppdf(x,0);
%     F = empcdf(x,0);
%     myObj = rvEmpiricalInfo(xi,f,F,0);
%     for ii=1:M
%         X(ii,2) = myObj.icdf(U(ii,2));
%     end
    X(:,2) = unifinv(U(:,2),-2,2);
end

X = X(randperm(M),:);

X_xform = X; X_xform(:,1) = continueRv(X(:,1));
U_xform = pobs(X_xform);

% now look at the pseudo-observations that are created to see the best way
% to create them
subplot(2,2,1); scatter(U(:,1),U(:,2)); grid on; title('U Actual');
subplot(2,2,2); scatter(U_xform(:,1),U_xform(:,2)); grid on; title('U Estimated');
subplot(2,2,3); scatter(X(:,1),X(:,2)); grid on; title('X Actual');
subplot(2,2,4); scatter(X_xform(:,1),X_xform(:,2)); grid on; title('X Continued');


%% test how HCBN is calculating LL values
clear;
clc;

% setup global parameters
D = 2;

%       A-->B
aa = 1; bb = 2;
dag = zeros(D,D);
dag(aa,bb) = 1;
discreteNodes = [aa];
nodeNames = {'A', 'B'};
bntPath = '../bnt';
discreteNodeNames = {'A'};
trainVecSize = 100:100:500; testSize = 500;
totalNumMC = 20;

continuousType = 'other';
discreteType = {};
nodeA = [0.4 0.3 0.2 0.1]; discreteType{1} = nodeA;
       
llMat = zeros(4,length(trainVecSize),totalNumMC);

for mcSim=1:totalNumMC
    
    idx = 1;
    for M=trainVecSize
        % Generate the synthetic data set
        [X, trueLLVec] = genSynthData3(discreteType, continuousType, M+testSize);
        X_train = X(1:M, :); X_test = X(M+1:end,:);
        trueLLTest = sum(trueLLVec(M+1:end));

        hcbnObj = hcbn(bntPath, X_train, nodeNames, discreteNodeNames, dag); hcbnLL = hcbnObj.dataLogLikelihood(X_test);
        clgObj = clgbn(X_train, discreteNodes, dag); clgLL = clgObj.dataLogLikelihood(X_test);
        mteObj = mtebn(X_train, discreteNodes, dag); mteLL = mteObj.dataLogLikelihood(X_test);

        fprintf('M=%d trueLL=%0.02f hcbnLL=%0.02f clgLL=%0.02f mteLL=%0.02f\n', ...
            M, trueLLTest, hcbnLL, clgLL, mteLL);

        llMat(1,idx,mcSim) = trueLLTest; 
        llMat(2,idx,mcSim) = hcbnLL;
        llMat(3,idx,mcSim) = clgLL; 
        llMat(4,idx,mcSim) = mteLL;
        
        idx = idx + 1;
    end
end

llAvg = mean(llMat,3);

fprintf('M=100 AVG(trueLL)=%0.02f AVG(hcbnLL)=%0.02f AVG(clgLL)=%0.02f AVG(mteLL)=%0.02f\n', ...
    llAvg(1,1), llAvg(2,1), llAvg(3,1), llAvg(4,1));
fprintf('M=200 AVG(trueLL)=%0.02f AVG(hcbnLL)=%0.02f AVG(clgLL)=%0.02f AVG(mteLL)=%0.02f\n', ...
    llAvg(1,2), llAvg(2,2), llAvg(3,2), llAvg(4,2));
fprintf('M=300 AVG(trueLL)=%0.02f AVG(hcbnLL)=%0.02f AVG(clgLL)=%0.02f AVG(mteLL)=%0.02f\n', ...
    llAvg(1,3), llAvg(2,3), llAvg(3,3), llAvg(4,3));
fprintf('M=400 AVG(trueLL)=%0.02f AVG(hcbnLL)=%0.02f AVG(clgLL)=%0.02f AVG(mteLL)=%0.02f\n', ...
    llAvg(1,4), llAvg(2,4), llAvg(3,4), llAvg(4,4));
fprintf('M=500 AVG(trueLL)=%0.02f AVG(hcbnLL)=%0.02f AVG(clgLL)=%0.02f AVG(mteLL)=%0.02f\n', ...
    llAvg(1,5), llAvg(2,5), llAvg(3,5), llAvg(4,5));