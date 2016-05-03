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

%% synthetic data simulation for HCBN only to catch and squish bugs
clear;
clc;

fid = fopen('/home/kiran/ownCloud/PhD/sim_results/synthDataSim.diary', 'w');

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
bntPath = '../bnt';
discreteNodeNames = {'A','B'};

hcbn_K = 25;        % WARNING - change this in actual HCBN code for now :(

% Generate the synthetic data set

discreteType = {};
nodeA = [0.4 0.3 0.2 0.1]; discreteType{1} = nodeA;
nodeB = [0.6 0.1 0.05 0.25]; discreteType{2} = nodeB;

% instantiate the CLG object
numMCSims = 5; numTest = 1000;
trainVecSize = 1000;
M = max(trainVecSize)+numTest; 
% (1,:) -> CLG-Gaussian, (2,:) -> MTE-Gaussian, (3,:) -> HCBN-Gaussian
% (4,:) -> CLG-Other     (5,:) -> MTE-Other     (6,:) -> HCBN-Other
llValMat = zeros(6,length(trainVecSize),numMCSims);   

dispstat('','init'); % One time only initialization
dispstat(sprintf('Begining the simulation...'),'keepthis','timestamp');
progressIdx = 1;
for mcSimNumber=1:numMCSims
    idx = 1;
    for numTrain=trainVecSize

        numTotalLoops = numMCSims*length(trainVecSize);
        progress = (progressIdx/numTotalLoops)*100;
        progressStr = sprintf('ABC Progress: Training Size=%d %0.02f%%',numTrain, progress);
        dispstat(progressStr,'timestamp');
        fprintf(fid, '%s\n', progressStr);

%         continuousType = 'Gaussian';
%         X = genSynthData2(discreteType, continuousType, M);
%         X_train_full = X(1:max(trainVecSize),:);
%         X_test = X(max(trainVecSize)+1:end,:);
%         X_train = X_train_full(1:numTrain,:);
%         hcbnObj = hcbn(bntPath, X_train, nodeNames, discreteNodeNames, dag);
%         hcbnObj.setSimNum(mcSimNumber); hcbnObj.setTypeName('Gaussian');
%         [hcbnLLVal_Gaussian,Rc_val_mat_Gaussian,Rc_num_mat_Gaussian,Rc_den_mat_Gaussian] = hcbnObj.dataLogLikelihood(X_test);
%         clgObj = clgbn(X_train, discreteNodes, dag); clgLLVal_Gaussian = clgObj.dataLogLikelihood(X_test);
%         mteObj = mtebn(X_train, discreteNodes, dag); mteLLVal_Gaussian = mteObj.dataLogLikelihood(X_test);
        
        continuousType = 'other';
        X = genSynthData2(discreteType, continuousType, M);
        X_train_full = X(1:max(trainVecSize),:);
        X_test = X(max(trainVecSize)+1:end,:);
        X_train = X_train_full(1:numTrain,:);
        hcbnObj = hcbn(bntPath, X_train, nodeNames, discreteNodeNames, dag); 
        hcbnObj.setSimNum(mcSimNumber); hcbnObj.setTypeName('Other');
        [hcbnLLVal_Other,Rc_val_mat_Other,Rc_num_mat_Other,Rc_den_mat_Other] = hcbnObj.dataLogLikelihood(X_test);
        clgObj = clgbn(X_train, discreteNodes, dag); clgLLVal_Other = clgObj.dataLogLikelihood(X_test);
        mteObj = mtebn(X_train, discreteNodes, dag); mteLLVal_Other = mteObj.dataLogLikelihood(X_test);
        
%         if(mean(Rc_val_mat_Gaussian(1,:))~=1 || mean(Rc_val_mat_Gaussian(2,:))~=1 || ...
%            mean(Rc_num_mat_Gaussian(1,:))~=1 || mean(Rc_num_mat_Gaussian(2,:))~=1 || ...
%            mean(Rc_den_mat_Gaussian(1,:))~=1 || mean(Rc_den_mat_Gaussian(2,:))~=1 || ...
%            mean(Rc_val_mat_Other(1,:))~=1 || mean(Rc_val_mat_Other(2,:))~=1 || ...
%            mean(Rc_num_mat_Other(1,:))~=1 || mean(Rc_num_mat_Other(2,:))~=1 || ...
%            mean(Rc_den_mat_Other(1,:))~=1 || mean(Rc_den_mat_Other(2,:))~=1 )
%             error('Node A and/or B error!');
%         end
        
%         subplot(1,2,1); plot(1:numTest, Rc_val_mat_Gaussian(3,:), 'b*', ...
%                              1:numTest, Rc_num_mat_Gaussian(3,:), 'r', ...
%                              1:numTest, Rc_den_mat_Gaussian(3,:), 'k'); grid on;        
%                          title(sprintf('Gaussian - Node C - %f',hcbnLLVal_Gaussian));
%         
%         subplot(1,2,2); plot(1:numTest, Rc_val_mat_Other(3,:), 'b*', ...
%                              1:numTest, Rc_num_mat_Other(3,:), 'r', ...
%                              1:numTest, Rc_den_mat_Other(3,:), 'k'); grid on;        
%                          title(sprintf('Other - Node C - %f',hcbnLLVal_Other));
%         pause(1);
        
%         llValMat(1, idx, mcSimNumber) = clgLLVal_Gaussian;
%         llValMat(2, idx, mcSimNumber) = mteLLVal_Gaussian;
%         llValMat(3, idx, mcSimNumber) = hcbnLLVal_Gaussian;
        llValMat(4, idx, mcSimNumber) = clgLLVal_Other;
        llValMat(5, idx, mcSimNumber) = mteLLVal_Other;
        llValMat(6, idx, mcSimNumber) = hcbnLLVal_Other;
                
        idx = idx + 1;
        progressIdx = progressIdx + 1;
    end

end
dispstat('Finished.','keepprev');

llValsAvg = mean(llValMat,3)

% fig1 = figure;

% subplot(2,1,1); 
% hold on;
% plot(trainVecSize, llValsAvg(1,:)); xlabel('Training Vector Size'); ylabel('Gaussian Data'); grid on
% plot(trainVecSize, llValsAvg(2,:)); 
% plot(trainVecSize, llValsAvg(3,:)); 
% legend('CLG','MTE', sprintf('HCBN - K=%d', hcbn_K));
% title('A->C B->C Structure')

% subplot(2,1,2); 
% hold on;
% plot(trainVecSize, llValsAvg(4,:)); xlabel('Training Vector Size'); ylabel('Non-Gaussian Likelihoods'); grid on
% plot(trainVecSize, llValsAvg(5,:)); 
% plot(trainVecSize, llValsAvg(6,:)); 
% legend('CLG','MTE', sprintf('HCBN - K=%d', hcbn_K));
% 
% set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.
% print(sprintf('/home/kiran/ownCloud/PhD/sim_results/simpleABC'),'-dpng')
% close(fig1);
fclose(fid);

%%

% synthetic data simulation for HCBN, MTE and CLG
clear;
clc;
fid = fopen('/home/kiran/ownCloud/PhD/sim_results/synthDataSim.diary', 'a');

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

hcbn_K = 25;        % WARNING - change this in actual HCBN code for now :(

% Generate the synthetic data set

discreteType = {};
nodeA = [0.4 0.3 0.2 0.1]; discreteType{1} = nodeA;
nodeB = [0.6 0.1 0.05 0.25]; discreteType{2} = nodeB;

% perform CLG/HCBN/MTE modeling, parametric to train/test size

% instantiate the CLG object
numMCSims = 100; numTest = 1000;
trainVecSize = 100:100:2000;
M = max(trainVecSize)+numTest; 
% (1,:) -> CLG-Gaussian, (2,:) -> MTE-Gaussian, (3,:) -> HCBN-Gaussian
% (4,:) -> CLG-Other     (5,:) -> MTE-Other     (6,:) -> HCBN-Other
llValMat = zeros(6,length(trainVecSize),numMCSims);   

dispstat('','init'); % One time only initialization
dispstat(sprintf('Begining the simulation...'),'keepthis','timestamp');
progressIdx = 1;
for mcSimNumber=1:numMCSims
    idx = 1;
    for numTrain=trainVecSize

        numTotalLoops = numMCSims*length(trainVecSize);
        progress = (progressIdx/numTotalLoops)*100;
        progressStr = sprintf('ABCDE Progress: Training Size=%d %0.02f%%',numTrain, progress);
        dispstat(progressStr,'timestamp');
        fprintf(fid, '%s\n', progressStr);

        continuousType = 'Gaussian';
        X = genSynthData(discreteType, continuousType, M);
        X_train_full = X(1:max(trainVecSize),:);
        X_test = X(max(trainVecSize)+1:end,:);
        X_train = X_train_full(1:numTrain,:);

        clgObj = clgbn(X_train, discreteNodes, dag); clgLLVal_Gaussian = clgObj.dataLogLikelihood(X_test);
        mteObj = mtebn(X_train, discreteNodes, dag); mteLLVal_Gaussian = mteObj.dataLogLikelihood(X_test);
        hcbnObj = hcbn(bntPath, X_train, nodeNames, discreteNodeNames, dag); hcbnLLVal_Gaussian = hcbnObj.dataLogLikelihood(X_test);

        continuousType = 'other';
        X = genSynthData(discreteType, continuousType, M);
        X_train_full = X(1:max(trainVecSize),:);
        X_test = X(max(trainVecSize)+1:end,:);
        X_train = X_train_full(1:numTrain,:);

        clgObj = clgbn(X_train, discreteNodes, dag); clgLLVal_Other = clgObj.dataLogLikelihood(X_test);
        mteObj = mtebn(X_train, discreteNodes, dag); mteLLVal_Other = mteObj.dataLogLikelihood(X_test);
        hcbnObj = hcbn(bntPath, X_train, nodeNames, discreteNodeNames, dag); hcbnLLVal_Other = hcbnObj.dataLogLikelihood(X_test);

        llValMat(1, idx, mcSimNumber) = clgLLVal_Gaussian;
        llValMat(2, idx, mcSimNumber) = mteLLVal_Gaussian;
        llValMat(3, idx, mcSimNumber) = hcbnLLVal_Gaussian;
        llValMat(4, idx, mcSimNumber) = clgLLVal_Other;
        llValMat(5, idx, mcSimNumber) = mteLLVal_Other;
        llValMat(6, idx, mcSimNumber) = hcbnLLVal_Other;
                
        idx = idx + 1;
        progressIdx = progressIdx + 1;
    end
    
end

dispstat('Finished.','keepprev');

llValsAvg = mean(llValMat,3);

fig1 = figure;

subplot(2,1,1); 
hold on;
plot(trainVecSize, llValsAvg(1,:)); xlabel('Training Vector Size'); ylabel('Gaussian Data'); grid on
plot(trainVecSize, llValsAvg(2,:)); 
plot(trainVecSize, llValsAvg(3,:)); 
legend('CLG','MTE', sprintf('HCBN - K=%d', hcbn_K));
title('A->C A->D B->D B->E Structure')

subplot(2,1,2); 
hold on;
plot(trainVecSize, llValsAvg(4,:)); xlabel('Training Vector Size'); ylabel('Non-Gaussian Likelihoods'); grid on
plot(trainVecSize, llValsAvg(5,:)); 
plot(trainVecSize, llValsAvg(6,:)); 
legend('CLG','MTE', sprintf('HCBN - K=%d', hcbn_K));

set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.
print(sprintf('/home/kiran/ownCloud/PhD/sim_results/simpleABCDE'),'-dpng')
close(fig1);

fclose(fid);

% %% 
% clear;
% clc;
% 
% load('/home/kiran/ownCloud/PhD/sim_results/llValMat_100.mat');
% numMCSims = 100;
% trainVecSize = 500:250:2000;
% 
% % average the monte-carlo simulations
% llValMatIdxs = isfinite(llValMat);
% llValMatFiniteIdxs = [];
% % find MC sim number which doesn't have any NaN's or +/- INFs
% for mcSimNumber=1:numMCSims
%     if(sum(any(llValMatIdxs(:,:,mcSimNumber)==0))==0)
%         llValMatFiniteIdxs = [llValMatFiniteIdxs mcSimNumber];
%     else
%         fprintf('Found bad at %d\n', mcSimNumber);
%     end
% end
% 
% llValsAvg = mean(llValMat(:,:,llValMatFiniteIdxs),3);
% gaussRef = llValsAvg(1,:);
% otherRef = llValsAvg(4,:);
% set(gca,'fontsize',30)
% hold on;
% plot(trainVecSize, gaussRef./llValsAvg(1,:), 'b*-.', 'LineWidth',2);
% plot(trainVecSize, gaussRef./llValsAvg(2,:), 'r*-.', 'LineWidth',2);
% plot(trainVecSize, gaussRef./llValsAvg(3,:), 'k*-.', 'LineWidth',2);
% plot(trainVecSize, otherRef./llValsAvg(4,:), 'b+-.', 'LineWidth',2);
% plot(trainVecSize, otherRef./llValsAvg(5,:), 'r+-.', 'LineWidth',2);
% plot(trainVecSize, otherRef./llValsAvg(6,:), 'k+-.', 'LineWidth',2);
% grid on; xlabel('# Training Samples'); 
% legend('CLG (Gaussian)', ...
%        'MTE (Gaussian)', ...
%        'HCBN (Gaussian)', ...
%        'CLG (Other)', ...
%        'MTE (Other)', ...
%        'HCBN (Other)')
% hold off;
% save('/home/kiran/ownCloud/PhD/synthDataSim.mat')
% 
% toc