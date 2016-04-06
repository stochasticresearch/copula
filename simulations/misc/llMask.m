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

% A script to help us determine if there is a likelihood masking problem

%% Compare what the generative models's density of the multimodal distribution looks like
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

alpha = 4;

discreteType = {};
nodeA = [0.4 0.3 0.2 0.1]; discreteType{1} = nodeA;

% numTrainVec = 250:250:2000;
numTrainInt = 100;
numTrainVec = 100:numTrainInt:300;
numMC = 3;

trueDistConditionalLLMat = zeros(length(numTrainVec),numMC); clgConditionalLLMat = zeros(length(numTrainVec),numMC); mteConditionalLLMat = zeros(length(numTrainVec),numMC);
hcbnConditionalLLMat = zeros(length(numTrainVec),numMC);
trueDistLLMat = zeros(length(numTrainVec),numMC); clgLLMat = zeros(length(numTrainVec),numMC); mteLLMat = zeros(length(numTrainVec),numMC);
hcbnLLMat = zeros(length(numTrainVec), numMC);

LOG_MIN = 1e-5;
dispstat('','init'); % One time only initialization
dispstat(sprintf('Begining the simulation...'),'keepthis','timestamp');
for numTrain=numTrainVec
    numTest = numTrain;
    for mcSimNum=1:numMC
        
        progressStr = sprintf('Progress: numTrain=%d mcSimNum=%d',numTrain, mcSimNum);
%         dispstat(progressStr,'timestamp');

        M = max(numTrain)+numTest; 

        continuousType = 'other';
        X = genSynthData3(discreteType, continuousType, M);

        X_train_full = X(1:max(numTrain),:);
        xtestIdxStart = max(numTrain)+1;
        X_test = X(xtestIdxStart:end,:);
        X_train = X_train_full(1:numTrain,:);
        mteObj = mte(X_train, discreteNodes, dag); mteLLVal_Other = mteObj.dataLogLikelihood(X_test);
        clgObj = clg(X_train, discreteNodes, dag); clgLLVal_Other = clgObj.dataLogLikelihood(X_test);
        hcbnObj = hcbn(bntPath, X_train, nodeNames, discreteNodeNames, dag); hcbnTrueLL = clgObj.dataLogLikelihood(X_test);       
        
%         [uU1, uU2] = ndgrid(linspace(0,1,100));
%         actualCopula = reshape(copulapdf('Gumbel', [uU1(:) uU2(:)], 4) , 100, 100);
%         estimatedCopula = hcbnObj.copulaFamilies{2}.c;
%         subplot(1,2,1);
%         surf(uU1,uU2,estimatedCopula); grid on; title('Estimated Copula');
%         subplot(1,2,2);
%         surf(uU1,uU2,actualCopula); grid on; title('Actual Copula');        
%         pause;
        
        % create the "true" conditional distribution
        trueConditionalDistribution = cell(1,4);
        for combo=1:4
            X_continuous_subset = [];
            for jj=1:numTrain
                if(X_train(jj,1)==combo)
                    X_continuous_subset = [X_continuous_subset; X_train(jj,2)];
                end
            end
            % estimate the empirical PDF/CDF for this
            isdiscrete = 0;
            [f,xi] = emppdf(X_continuous_subset, isdiscrete);
            F = empcdf(X_continuous_subset, isdiscrete);
            trueConditionalDistribution{combo} = rvEmpiricalInfo(xi, f, F);
        end
        
        isdiscrete = 0;
        [f,xi] = emppdf(X_train(:,2), isdiscrete);
        F = empcdf(X_train(:,2), isdiscrete);
        continuousDistInfo = rvEmpiricalInfo(xi,f,F);
        
        isdiscrete = 1;
        [f,xi] = emppdf(X_train(:,1), isdiscrete);
        F = empcdf(X_train(:,1), isdiscrete);
        discreteDistInfo = rvEmpiricalInfo(xi,f, F);
        
        trueDistConditionalLL = 0; clgConditionalLL = 0; mteConditionalLL = 0;
        trueDistLL = 0; clgLL = 0; mteLL = 0; 
        hcbnConditionalLL = 0; hcbnLL = 0;
        for ii=1:numTest
            xi = X_test(ii,:);
            comboVal = xi(1);
            trueDistInfo = trueConditionalDistribution{comboVal};
            clgContinuousDistInfo = clgObj.bnParams{2}{comboVal};
            mteContinuousDistInfo = mteObj.bnParams{2}{comboVal}.mte_info;

            continuousVal = xi(2);

            trueDistConditionalLikelihood = trueDistInfo.queryDensity(continuousVal);
            clgConditionalLikelihood = normpdf(continuousVal,clgContinuousDistInfo.Mean,clgContinuousDistInfo.Covariance);
            mteConditionalLikelihood = mteContinuousDistInfo.queryDensity(continuousVal);
            % compute hcbn conditional log likelihood
            u = [discreteDistInfo.queryDistribution(comboVal) continuousDistInfo.queryDistribution(continuousVal)];
            u = fixU(u);
%             copVal = empcopula_val(estimatedCopula, u);
%             copVal = empcopula_val(actualCopula, u);
            copVal = copulapdf('Gumbel', u, alpha);
            if(copVal<LOG_MIN)
                copVal=LOG_MIN;
            end
            f_x2 = continuousDistInfo.queryDensity(continuousVal);
            if(f_x2<LOG_MIN)
                f_x2=LOG_MIN;
            end
            hcbnConditionalLikelihood = f_x2*copVal;

            trueDistLikelihood = trueDistConditionalLikelihood*nodeA(comboVal);
            clgLikelihood = clgConditionalLikelihood*clgObj.bnParams{1}.density(comboVal);
            mteLikelihood = mteConditionalLikelihood*mteObj.bnParams{1}.density(comboVal);
            hcbnLikelihood = hcbnConditionalLikelihood*hcbnObj.empInfo{1}.density(comboVal);

            % plot all 4 distributions on top of each other
            % compute the HCBN conditional density
            hcbnConditionalDensity = zeros(1,length(trueDistInfo.density));
            for zzz=1:length(hcbnConditionalDensity)
                uu = [discreteDistInfo.queryDistribution(comboVal) continuousDistInfo.queryDistribution(trueDistInfo.domain(zzz))];
                f_x2_val = continuousDistInfo.queryDensity(trueDistInfo.domain(zzz));
%                 hcbnConditionalDensity(zzz) = empcopula_val(estimatedCopula.c, uu)*f_x2_val;
                hcbnConditionalDensity(zzz) = copulapdf('Gumbel', uu,4)*f_x2_val;
            end
            
            fig1 = figure(1);
            plot(trueDistInfo.domain, trueDistInfo.density, ...
                 trueDistInfo.domain, normpdf(trueDistInfo.domain, clgContinuousDistInfo.Mean, clgContinuousDistInfo.Covariance), ...
                 mteContinuousDistInfo.domain, mteContinuousDistInfo.density, ...
                 trueDistInfo.domain, hcbnConditionalDensity);
            grid on; axis([min(trueDistInfo.domain) max(trueDistInfo.domain) 0 1]); 
            hold on;
            % plot the data points and their respective probabilities calculated    
            plot(continuousVal,trueDistLikelihood, '+'); 
            plot(continuousVal,trueDistConditionalLikelihood, '^');
            plot(continuousVal,clgLikelihood, '*'); 
            plot(continuousVal,clgConditionalLikelihood, 'v');
            plot(continuousVal,mteLikelihood, 'o'); 
            plot(continuousVal,mteConditionalLikelihood, 'h');
            plot(continuousVal,hcbnLikelihood, 'x');
            plot(continuousVal,hcbnConditionalLikelihood, 'p');
            
            legend('True','CLG','MTE', 'HCBN', 'TrueLL', 'TrueConditionaLL', 'CLG LL', 'CLG ConditionalLL', ...
                'MTE LL', 'MTE_Conditional LL', 'HCBN LL', 'HCBN Conditional LL');
            
            pause;
            clf(fig1);

            if(trueDistConditionalLikelihood<LOG_MIN)
                trueDistConditionalLikelihood=LOG_MIN;
            end
            if(clgConditionalLikelihood<LOG_MIN)
                clgConditionalLikelihood=LOG_MIN;
            end
            if(mteConditionalLikelihood<LOG_MIN)
                mteConditionalLikelihood=LOG_MIN;
            end
            if(hcbnConditionalLikelihood<LOG_MIN)
                hcbnConditionalLikelihood=LOG_MIN;
            end
            if(trueDistLikelihood<LOG_MIN)
                trueDistLikelihood=LOG_MIN;
            end
            if(clgLikelihood<LOG_MIN)
                clgLikelihood=LOG_MIN;
            end
            if(mteLikelihood<LOG_MIN)
                mteLikelihood=LOG_MIN;
            end
            if(hcbnLikelihood<LOG_MIN)
                hcbnLikelihood=LOG_MIN;
            end
        
            trueDistConditionalLL = trueDistConditionalLL + log(trueDistConditionalLikelihood);
            clgConditionalLL = clgConditionalLL + log(clgConditionalLikelihood);
            mteConditionalLL = mteConditionalLL + log(mteConditionalLikelihood);
            hcbnConditionalLL = hcbnConditionalLL + log(hcbnConditionalLikelihood);
            
%             fprintf('Conditional >> true=%0.02f clg=%0.02f mte=%0.02f hcbn=%0.02f\n', ...
%                 trueDistConditionalLL, clgConditionalLL, mteConditionalLL, hcbnConditionalLL);
            
            trueDistLL = trueDistLL + log(trueDistLikelihood);
            clgLL = clgLL + log(clgLikelihood);
            mteLL = mteLL + log(mteLikelihood);
            hcbnLL = hcbnLL + log(hcbnLikelihood);
        end
        
        trueDistConditionalLLMat(numTrain/numTrainInt, mcSimNum) = trueDistConditionalLL;
        clgConditionalLLMat(numTrain/numTrainInt, mcSimNum) = clgConditionalLL;
        mteConditionalLLMat(numTrain/numTrainInt, mcSimNum) = mteConditionalLL;
        hcbnConditionalLLMat(numTrain/numTrainInt, mcSimNum) = hcbnConditionalLL;
        
        trueDistLLMat(numTrain/numTrainInt, mcSimNum) = trueDistLL;
        clgLLMat(numTrain/numTrainInt, mcSimNum) = clgLL;
        mteLLMat(numTrain/numTrainInt, mcSimNum) = mteLL;
        hcbnLLMat(numTrain/numTrainInt, mcSimNum) = hcbnLL;
        fprintf('hcbnLLCalc=%f hcbnLLModule=%f\n', hcbnLL, hcbnTrueLL);
    end

end
dispstat('Finished.','keepprev');

trueDistConditionalLL_mean = mean(trueDistConditionalLLMat,2);
mteConditionalLL_mean = mean(mteConditionalLLMat,2);
clgConditionalLL_mean = mean(clgConditionalLLMat,2);
hcbnConditionalLL_mean = mean(hcbnConditionalLLMat,2);
trueDistLL_mean = mean(trueDistLLMat,2);
mteLL_mean = mean(mteLLMat,2);
clgLL_mean = mean(clgLLMat,2);
hcbnLL_mean = mean(hcbnLLMat, 2);

mteConditionalDiff = ((trueDistConditionalLL_mean-mteConditionalLL_mean)./abs(trueDistConditionalLL_mean)) * 100;
clgConditionalDiff = ((trueDistConditionalLL_mean-clgConditionalLL_mean)./abs(trueDistConditionalLL_mean)) * 100;
hcbnConditionalDiff = ((trueDistConditionalLL_mean-hcbnConditionalLL_mean)./abs(trueDistConditionalLL_mean)) * 100;

mteTotalDiff = ((trueDistLL_mean-mteLL_mean)./abs(trueDistLL_mean)) * 100;
clgTotalDiff = ((trueDistLL_mean-clgLL_mean)./abs(trueDistLL_mean)) * 100;
hcbnTotalDiff = ((trueDistLL_mean-hcbnLL_mean)./abs(trueDistLL_mean)) * 100;

plot(numTrainVec,mteConditionalDiff,...
    numTrainVec,clgConditionalDiff,...
    numTrainVec,hcbnConditionalDiff, ...
    numTrainVec,mteTotalDiff,...
    numTrainVec,clgTotalDiff,...
    numTrainVec,hcbnTotalDiff); grid on;
xlabel('Num Training Points'); ylabel('% Diff from True LL'); title(sprintf('%d MC Simulations', numMC));
legend('MTE Conditional LL', 'CLG Conditional LL', 'HCBN Conditional LL', 'MTE Full LL', 'CLG Full LL', 'HCBN Full LL');

% fprintf('Condtional(True LL)=%f Conditional(CLG LL)=%f Conditional(MTE LL)=%f\n', ...
%          trueDistConditionalLL_mean, clgConditionalLL_mean, mteConditionalLL_mean);
% % fprintf('True LL=%f CLG LL=%f MTE LL=%f CLG_LL_calc=%f MTE_LL_calc=%f\n', ...
% %          trueDistLL, clgLL, mteLL, clgLLVal_Other, mteLLVal_Other);
% fprintf('True LL=%f CLG LL=%f MTE LL=%f\n', ...
%          trueDistLL_mean, clgLL_mean, mteLL_mean);

% fprintf('Conditional(CLG Diff)=%0.02f Conditional(MTE Diff)=%0.02f CLG Diff=%0.02f MTE Diff=%0.02f\n', ...
%     clgConditionalDiff, mteConditionalDiff, clgTotalDiff, mteTotalDiff);