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

% A script to help us investigate why MTE has very seemingly poor
% performance!

%% Test the LL calculation for MTE and Gaussian w/ only one node that has a multimodal distribution
clear;
clc;

numMCSims = 20;
numTrain = 250:250:2000;

dag = 0;
discreteNodes = [];

mteLLVec = zeros(1, numMCSims); clgLLVec = zeros(1, numMCSims);
for trainVecSize=numTrain
    for mcSimNumber=1:numMCSims

        M = trainVecSize;
        x = [normrnd(-2,0.3,M/2,1); normrnd(2,0.5,M/2,1)];
        x = x(randperm(M));

        X_train = x(1:M/2);
        X_test = x(M/2+1:end);

        mteObj = mtebn(X_train, discreteNodes, dag); mteLLVal_Other = mteObj.dataLogLikelihood(X_test);
        clgObj = clgbn(X_train, discreteNodes, dag); clgLLVal_Other = clgObj.dataLogLikelihood(X_test);

        mteLLVec(mcSimNumber) = mteLLVal_Other;
        clgLLVec(mcSimNumber) = clgLLVal_Other;
    end

    mean_mteLL = mean(mteLLVec); mean_clgLL = mean(clgLLVec);
    percentageDiff = ((mean_clgLL-mean_mteLL)/mean_clgLL) * 100;
    
    progressStr = sprintf('trainVecSize=%d AVG(MTE_Other_LL)=%f AVG(CLG_Other_LL)=%f DIFF=%0.02f\n', ...
            trainVecSize, mean(mteLLVec), mean(clgLLVec), percentageDiff);
    fprintf(progressStr);
end

%%
clear;
clc;

% generate a small synthetic dataset
% setup global parameters
D = 2;

%       A-->B      

aa = 1; bb = 2;
dag = zeros(D,D);
dag(aa,bb) = 1;
discreteNodes = [aa];
nodeNames = {'A', 'B'};
bntPath = '../bnt'; addpath(genpath(bntPath));
discreteNodeNames = {'A'};

discreteType = {};
nodeA = [0.4 0.3 0.2 0.1]; discreteType{1} = nodeA;

numMCSims = 10; numTest = 1000;
numTrainVec = 500:500:1000;
M = max(numTrainVec)+numTest; 

llValMat = zeros(4,length(numTrainVec),numMCSims);  
for mcSimNumber=1:numMCSims
    idx = 1;
    for numTrain=numTrainVec

        continuousType = 'Gaussian';
        X = genSynthData3(discreteType, continuousType, M);
        X_train_full = X(1:max(numTrain),:);
        X_test = X(max(numTrain)+1:end,:);
        X_train = X_train_full(1:numTrain,:);
        mteObj = mtebn(X_train, discreteNodes, dag); mteLLVal_Gaussian = mteObj.dataLogLikelihood(X_test);
        clgObj = clgbn(X_train, discreteNodes, dag); clgLLVal_Gaussian = clgObj.dataLogLikelihood(X_test);
        
        continuousType = 'other';
        X = genSynthData3(discreteType, continuousType, M);
        X_train_full = X(1:max(numTrain),:);
        X_test = X(max(numTrain)+1:end,:);
        X_train = X_train_full(1:numTrain,:);
        mteObj = mtebn(X_train, discreteNodes, dag); mteLLVal_Other = mteObj.dataLogLikelihood(X_test);
        clgObj = clgbn(X_train, discreteNodes, dag); clgLLVal_Other = clgObj.dataLogLikelihood(X_test);
        
        progressStr = sprintf('numTrain=%d MTE_Gaussian_LL=%f CLG_Gaussian_LL=%f MTE_Other_LL=%f CLG_Other_LL=%f \n', ...
                numTrain, mteLLVal_Gaussian, clgLLVal_Gaussian, mteLLVal_Other, clgLLVal_Other);
        fprintf(progressStr);

        llValMat(1, idx, mcSimNumber) = clgLLVal_Gaussian;
        llValMat(2, idx, mcSimNumber) = mteLLVal_Gaussian;
        llValMat(3, idx, mcSimNumber) = clgLLVal_Other;
        llValMat(4, idx, mcSimNumber) = mteLLVal_Other;

        idx = idx + 1;
    end

end
llValsAvg = mean(llValMat,3);
fprintf('numTrain=%d CLG_Gaussian_LL_Avg=%f MTE_Gaussian_LL_Avg=%f CLG_Other_LL_Avg=%f MTE_Other_LL_Avg=%f\n', ...
    500, llValsAvg(1,1), llValsAvg(2,1), llValsAvg(3,1), llValsAvg(4,1));
fprintf('numTrain=%d CLG_Gaussian_LL_Avg=%f MTE_Gaussian_LL_Avg=%f CLG_Other_LL_Avg=%f MTE_Other_LL_Avg=%f\n', ...
    1000, llValsAvg(1,2), llValsAvg(2,2), llValsAvg(3,2), llValsAvg(4,2));


     
%%
clear;
clc;

% generate a small synthetic dataset
% setup global parameters
D = 3;

%       A   B      
%        \ /
%         C
aa = 1; bb = 2; cc = 3;
dag = zeros(D,D);
dag(aa,cc) = 1;
dag(bb,cc) = 1;
discreteNodes = [aa bb];
nodeNames = {'A', 'B', 'C'};
bntPath = '../bnt'; addpath(genpath(bntPath));
discreteNodeNames = {'A','B'};

discreteType = {};
nodeA = [0.4 0.3 0.2 0.1]; discreteType{1} = nodeA;
nodeB = [0.6 0.1 0.05 0.25]; discreteType{2} = nodeB;

numMCSims = 10; numTest = 1000;
numTrain = 500;
M = max(numTrain)+numTest; 

llValMat = zeros(4,length(numTrain),numMCSims);   

% dispstat('','init'); % One time only initialization
% dispstat(sprintf('Begining the simulation...'),'keepthis','timestamp');
progressIdx = 1;
for mcSimNumber=1:numMCSims
    idx = 1;
    for numTrain=numTrain

        numTotalLoops = numMCSims*length(numTrain);
        progress = (progressIdx/numTotalLoops)*100;
%         progressStr = sprintf('MTE -- ABC Progress: Training Size=%d %0.02f%%',numTrain, progress);
%         dispstat(progressStr,'timestamp');

%         continuousType = 'Gaussian';
%         X = genSynthData2(discreteType, continuousType, M);
%         X_train_full = X(1:max(trainVecSize),:);
%         X_test = X(max(trainVecSize)+1:end,:);
%         X_train = X_train_full(1:numTrain,:);
%         mteObj = mtebn(X_train, discreteNodes, dag); mteLLVal_Gaussian = mteObj.dataLogLikelihood(X_test);
%         clgObj = clgbn(X_train, discreteNodes, dag); clgLLVal_Gaussian = clgObj.dataLogLikelihood(X_test);
        
        continuousType = 'other';
        X = genSynthData2(discreteType, continuousType, M);
        
        X_train_full = X(1:max(numTrain),:);
        X_test = X(max(numTrain)+1:end,:);
        X_train = X_train_full(1:numTrain,:);
        mteObj = mtebn(X_train, discreteNodes, dag); mteLLVal_Other = mteObj.dataLogLikelihood(X_test);
        clgObj = clgbn(X_train, discreteNodes, dag); clgLLVal_Other = clgObj.dataLogLikelihood(X_test);
        
%         progressStr = sprintf('numTrain=%d MTE_Gaussian_LL=%f CLG_Gaussian_LL=%f MTE_Other_LL=%f CLG_Other_LL=%f \n', ...
%                 numTrain, mteLLVal_Gaussian, clgLLVal_Gaussian, mteLLVal_Other, clgLLVal_Other);
        progressStr = sprintf('numTrain=%d MTE_Other_LL=%f CLG_Other_LL=%f \n', ...
                numTrain, mteLLVal_Other, clgLLVal_Other);
        fprintf(progressStr);
%         dispstat(progressStr,'timestamp');

%         llValMat(1, idx, mcSimNumber) = clgLLVal_Gaussian;
%         llValMat(2, idx, mcSimNumber) = mteLLVal_Gaussian;
        llValMat(3, idx, mcSimNumber) = clgLLVal_Other;
        llValMat(4, idx, mcSimNumber) = mteLLVal_Other;

        idx = idx + 1;
        progressIdx = progressIdx + 1;
    end

end
% dispstat('Finished.','keepprev');

llValsAvg = mean(llValMat,3);
fprintf('MTE_Other_LL_Avg=%f CLG_Other_LL_Avg=%f\n', llValsAvg(4), llValsAvg(3));