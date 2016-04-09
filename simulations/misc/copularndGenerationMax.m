% Function to see if copularnd ever generates any samples at the edges and
% what the max value it generates is, over some # of MC simulations

clear;
clc;

numMCSims = 1000;
alphaVec = 1:3:10;
RhoVecs_2D = cell(1,length(alphaVec)); 
RhoVecs_2D{1} = [1 -0.9; -0.9 1]; RhoVecs_2D{2} = [1 -0.65; -0.65 1];
RhoVecs_2D{3} = [1 0.35; 0.35 1]; RhoVecs_2D{4} = [1 0.1; 0.1 1];
copulaTypeVec = {'Frank', 'Gumbel', 'Clayton', 'Gaussian'}; 

M = 1000;

for copulaTypeIdx=1:length(copulaTypeVec)
    for alphaVecIdx=1:length(alphaVec)
        maxUVec = zeros(2,numMCSims);
        for mcSim=1:numMCSims
            copulaType = copulaTypeVec{copulaTypeIdx};
            if(strcmp(copulaType,'Gaussian'))
                alpha = RhoVecs_2D{alphaVecIdx};
                alphaValPrint = alpha(1,2);
            else
                alpha = alphaVec(alphaVecIdx);
                alphaValPrint = alpha;
            end
            U = copularnd(copulaType, alpha, M);
            maxUVec(1,mcSim) = max(U(:,1));
            maxUVec(2,mcSim) = max(U(:,2));
        end
        
        fprintf('copulaType=%s alpha=%f maxU1=%f maxU2=%f\n', ...
            copulaType, alphaValPrint, max(maxUVec(1,:)), max(maxUVec(2,:)));
    end
end