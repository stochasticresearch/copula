%******************************************************************************
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

% synthetic data simulation for HCBN, MTE and CLG
clear;
clc;

tic

% setup global parameters
D = 5;

%       A   B      
%      / \ / \
%     C   D   E
aa = 1; bb = 2; cc = 3; dd = 4; ee = 5;
dag = zeros(D,D);
dag(aa,[cc dd]) = 1;
dag(bb,[dd ee]) = 1;
discreteNodes = [aa bb];
nodeNames = {'A', 'B', 'C', 'D', 'E'};
bntPath = '../bnt';
discreteNodeNames = {'A','B'};

% Generate the synthetic data set

discreteType = {};
nodeA = [0.4 0.3 0.2 0.1]; discreteType{1} = nodeA;
nodeB = [0.6 0.1 0.05 0.25]; discreteType{2} = nodeB;

% perform CLG/HCBN/MTE modeling, parametric to train/test size

% instantiate the CLG object
numMCSims = 100; numTest = 1000;
trainVecSize = 500:250:2000;
M = max(trainVecSize)+numTest; 
% (1,:) -> CLG-Gaussian, (2,:) -> MTE-Gaussian, (3,:) -> HCBN-Gaussian
% (4,:) -> CLG-Other     (5,:) -> MTE-Other     (6,:) -> HCBN-Other
llValMat = zeros(6,length(trainVecSize),numMCSims);   
for mcSimNumber=1:numMCSims
    idx = 1;
    for numTrain=trainVecSize

        fprintf('Processing training size=%d -- MC Sim #=%d\n', numTrain, mcSimNumber);

        continuousType = 'Gaussian';
        X = genSynthData(discreteType, continuousType, M);
        X_train_full = X(1:max(trainVecSize),:);
        X_test = X(max(trainVecSize)+1:end,:);
        X_train = X_train_full(1:numTrain,:);

        clgObj = clg(X_train, discreteNodes, dag); clgLLVal_Gaussian = clgObj.dataLogLikelihood(X_test);
        mteObj = mte(X_train, discreteNodes, dag); mteLLVal_Gaussian = mteObj.dataLogLikelihood(X_test);
        hcbnObj = hcbn(bntPath, X_train, nodeNames, discreteNodeNames, dag); hcbnLLVal_Guassian = hcbnObj.hcbnLogLikelihood(X_test);

        continuousType = 'other';
        X = genSynthData(discreteType, continuousType, M);
        X_train_full = X(1:max(trainVecSize),:);
        X_test = X(max(trainVecSize)+1:end,:);
        X_train = X_train_full(1:numTrain,:);

        clgObj = clg(X_train, discreteNodes, dag); clgLLVal_Other = clgObj.dataLogLikelihood(X_test);
        mteObj = mte(X_train, discreteNodes, dag); mteLLVal_Other = mteObj.dataLogLikelihood(X_test);
        hcbnObj = hcbn(bntPath, X_train, nodeNames, discreteNodeNames, dag); hcbnLLVal_Other = hcbnObj.hcbnLogLikelihood(X_test);

        llValMat(1, idx, mcSimNumber) = clgLLVal_Gaussian;
        llValMat(2, idx, mcSimNumber) = mteLLVal_Gaussian;
        llValMat(3, idx, mcSimNumber) = hcbnLLVal_Guassian;
        llValMat(4, idx, mcSimNumber) = clgLLVal_Other;
        llValMat(5, idx, mcSimNumber) = mteLLVal_Other;
        llValMat(6, idx, mcSimNumber) = hcbnLLVal_Other;
                
        idx = idx + 1;
    end
    
    % save off results after every MC sim
    fnameToSave = sprintf('/home/kiran/ownCloud/PhD/sim_results/llValMat_%d.mat', mcSimNumber);
    save(fnameToSave, 'llValMat');

end

%% 
clear;
clc;

load('/home/kiran/ownCloud/PhD/sim_results/llValMat_100.mat');
numMCSims = 100;
trainVecSize = 500:250:2000;

% average the monte-carlo simulations
llValMatIdxs = isfinite(llValMat);
llValMatFiniteIdxs = [];
% find MC sim number which doesn't have any NaN's or +/- INFs
for mcSimNumber=1:numMCSims
    if(sum(any(llValMatIdxs(:,:,mcSimNumber)==0))==0)
        llValMatFiniteIdxs = [llValMatFiniteIdxs mcSimNumber];
    else
        fprintf('Found bad at %d\n', mcSimNumber);
    end
end

llValsAvg = mean(llValMat(:,:,llValMatFiniteIdxs),3);
gaussRef = llValsAvg(1,:);
otherRef = llValsAvg(4,:);
set(gca,'fontsize',30)
hold on;
plot(trainVecSize, gaussRef./llValsAvg(1,:), 'b*-.', 'LineWidth',2);
plot(trainVecSize, gaussRef./llValsAvg(2,:), 'r*-.', 'LineWidth',2);
plot(trainVecSize, gaussRef./llValsAvg(3,:), 'k*-.', 'LineWidth',2);
plot(trainVecSize, otherRef./llValsAvg(4,:), 'b+-.', 'LineWidth',2);
plot(trainVecSize, otherRef./llValsAvg(5,:), 'r+-.', 'LineWidth',2);
plot(trainVecSize, otherRef./llValsAvg(6,:), 'k+-.', 'LineWidth',2);
grid on; xlabel('# Training Samples'); 
legend('CLG (Gaussian)', ...
       'MTE (Gaussian)', ...
       'HCBN (Gaussian)', ...
       'CLG (Other)', ...
       'MTE (Other)', ...
       'HCBN (Other)')
hold off;
save('/home/kiran/ownCloud/PhD/synthDataSim.mat')

toc