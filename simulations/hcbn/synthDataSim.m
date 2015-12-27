% synthetic data simulation for HCBN, MTE and CLG
clear;
clc;

%% setup global parameters
M = 10000; D = 5;

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

%% Generate the synthetic data set

discreteType = {};
nodeA = [0.4 0.3 0.2 0.1]; discreteType{1} = nodeA;
nodeB = [0.6 0.1 0.05 0.25]; discreteType{2} = nodeB;
continuousType = 'Gaussian';
X = genSynthData(discreteType, continuousType, M);

X_train_full = X(1:9000,:);
X_test = X(9001:end,:);

%% perform CLG/HCBN/MTE modeling, parametric to train/test size

% instantiate the CLG object
trainVecSize = 500:500:5000;
llValVec = zeros(3,length(trainVecSize));   % (1,:) -> CLG, (2,:) -> HCBN, (3,:) -> MTE
idx = 1;
for numTrain=trainVecSize
    fprintf('Processing training size=%d\n', numTrain);
    X_train = X_train_full(1:numTrain,:);
    
    clgObj = clg(X_train, discreteNodes, dag);
    clgLLVal = clgObj.dataLogLikelihood(X_test);
    
    hcbnObj = hcbn(bntPath, X_train, nodeNames, discreteNodeNames, dag);
    hcbnLLval = hcbnObj.hcbnLogLikelihood(X_test);
    
    llValVec(1, idx) = clgLLVal;
    llValVec(2, idx) = hcbnLLval;
    idx = idx + 1;
end

plot(trainVecSize, llValVec(1,:), trainVecSize, llValVec(2,:)); 
grid on; xlabel('# Training Samples'); ylabel('Log-Likelihood')
legend('CLG', 'HCBN')