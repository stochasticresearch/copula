% synthetic data simulation for HCBN, MTE and CLG
clear;
clc;

discreteType = {};
nodeA = [0.4 0.3 0.2 0.1]; discreteType{1} = nodeA;
nodeB = [0.6 0.1 0.05 0.25]; discreteType{2} = nodeB;
continuousType = 'Gaussian';
M = 1000; D = 5;
X = genSynthData(discreteType, continuousType, M);

% instantiate the CLG object
%       A   B      
%      / \ / \
%     C   D   E
aa = 1; bb = 2; cc = 3; dd = 4; ee = 5;
dag = zeros(D,D);
dag(aa,[cc dd]) = 1;
dag(bb,[dd ee]) = 1;
discreteNodes = [aa bb];
clgObj = clg(X, discreteNodes, dag);