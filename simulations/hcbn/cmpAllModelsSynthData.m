function [llMat, llVarMat, llBiasMat, llMCCell] = cmpAllModelsSynthData( D, numMCSims, cfg, logFilename, plotOption, mVecInput )
%CMPALLMODELSSYNTHDATA - compares all the models for hybrid networks
% (CLG, HCBN, MTE, [CBN], [MULTINOMIAL]) with synthetic data generated with
% various different dependencies
%  Inputs:
%   D - the dimensionality of the data.  A value of D maps to a unique
%       graphical model from which the synthetic data was generated from.
%   numMCSims - the number of monte-carlo simulations to run each time to
%               generate average KL Div values
%   cfg - the configuration for the value of D, different configurations of
%         D test different configurations of the discrete variable
%
%  Outputs:
%   resultsMat - 
%     For D=2, resultsMat is a [4 x 4 x 4 x 8 x 7 x <CFG>] matrix of results, 
%     where each dimension is representative of the following:
%       1 >> the copula type.  (1)=Frank, (2)=Gumbel, (3)=Clayton, (4)=Gaussian
%       2 >> the type of distribution for the dependent RV that is
%            continuous, maps to:
%            (1)=Multinomial (2)=Uniform (3)=Gaussian (4)=ThickTailed
%       3 >> dependency amount
%            For archimedean copulas, 
%              (1)=1, (2)=4, (3)=7, (4)=10
%       4 >> the number of samples
%              (1)=250 (2)=500 (3)=750 (4)=1000
%              (5)=1250 (6)=1500 (7)=1750 (8)=2000
%       5 >> the model compared to generative model via KL Divergence
%              (1)-copula estimate w/ rank pseudo-obs, f2 estimate
%              (2)-copula estimate w/ rank pseudo-obs, f2 actual
%              (3)-copula actual, f2 estimate
%              (4)-hcbn (copula estimate w/ ecdf pseudo-obs, f2 estimate)
%              (5)-empirical kernel density estimate of data
%              (6)-mte
%              (7)-clg
%       6 >> the conditional value of X1 for which the KL Divergence was
%            estimated for X2.  The size of this depends on the
%            configuration input.  See below for different configurations
%            and the kind of data they generate
%
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

% define the parametrization variables
% parametrization variables
global mVec copulaTypeVec_2D alphaVec RhoVecs_2D RhoVecs_3D copulaTypeVec_3D 
global copulaTypeVec_4D RhoVecs_4D
global numModelsCompared numMC bntPath logFile K h continuousDistTypeVec
global plotFlag numTest
global HCBN_LL_MAT_IDX MTE_LL_MAT_IDX CLG_LL_MAT_IDX MULTINOMIAL_LL_MAT_IDX 
global CBN_LL_MAT_IDX CBN_GAUSSIAN_LL_MAT_IDX REF_LL_MAT_IDX 
global NUM_DISCRETE_INTERVALS
global HCBN_DEBUGALL_LL_MAT_IDX HCBN_DEBUGCOPULA_LL_MAT_IDX 
global HCBN_DEBUGEMPINFO_LL_MAT_IDX
global CDE_combinations C1C2C3_combinations dependency_combinations

if(nargin>4)
    plotFlag = plotOption;
else
    plotFlag = 0;
end
if(nargin>5)
    mVec = mVecInput;
else
    mVec = 250:250:1000;
end


K = 50; h = 0.05;      % beta kernel estimation parameters
NUM_DISCRETE_INTERVALS = 10;
bntPath = '../bnt'; addpath(genpath(bntPath));

copulaTypeVec_2D = {'Frank', 'Gaussian'};
copulaTypeVec_3D = {'Gaussian'};
copulaTypeVec_4D = {'Gaussian'};
alphaVec = [1 10];
RhoVecs_2D = cell(1,length(alphaVec)); 
RhoVecs_2D{1} = [1 0.1; 0.1 1]; RhoVecs_2D{2} = [1 -0.9; -0.9 1]; 
RhoVecs_2D{3} = [1 0.35; 0.35 1]; RhoVecs_2D{4} = [1 -0.65; -0.65 1];
RhoVecs_3D = cell(1,length(alphaVec));
RhoVecs_3D{1} = [1 0 .2; 0 1 -.8; .2 -.8 1];
RhoVecs_3D{2} = [1 0 .3; 0 1 -.6; .3 -.6 1];
RhoVecs_3D{3} = [1 0 .3; 0 1 -.1; .3 -.1 1];
RhoVecs_3D{4} = [1 0 -.3; 0 1 .1; -.3 .1 1];

RhoVecs_4D = cell(1,4);
RhoVecs_4D{1} = [1 0 0 .2; 0 1 0 -.8; 0 0 1 .3; 0.2 -0.8 0.3 1];
RhoVecs_4D{2} = [1 0 0 .5; 0 1 0 -.7; 0 0 1 .2; 0.5 -0.7 0.2 1];
RhoVecs_4D{3} = [1 0 0 .25; 0 1 0 -.4; 0 0 1 .6; 0.25 -0.4 0.6 1];

numTest = 1000; % the # of samples to generate to calculate likelihood
HCBN_LL_MAT_IDX = 1;
HCBN_DEBUGALL_LL_MAT_IDX = 2;
HCBN_DEBUGCOPULA_LL_MAT_IDX = 3;
HCBN_DEBUGEMPINFO_LL_MAT_IDX = 4;
MTE_LL_MAT_IDX = 5;
CLG_LL_MAT_IDX = 6;
MULTINOMIAL_LL_MAT_IDX = 7;
CBN_LL_MAT_IDX = 8;
CBN_GAUSSIAN_LL_MAT_IDX = 9;
REF_LL_MAT_IDX = 10;

continuousDistTypeVec = {'Gaussian', 'Uniform', 'Multimodal', 'ThickTailed'};
numModelsCompared = 10;
numMC = numMCSims;
logFile = logFilename;

% enumerate all the possible combinations for the 5-D test
CDE_combinations = cell(1,4);
CDE_combinations{1} = {'Gaussian', 'Gaussian', 'Gaussian'};
CDE_combinations{2} = {'Multimodal', 'Uniform', 'Multimodal'};
CDE_combinations{3} = {'Uniform', 'Uniform', 'ThickTailed'};
CDE_combinations{4} = {'Multimodal', 'Gaussian', 'Uniform'};

C1C2C3_combinations = cell(1,2);
C1C2C3_combinations{1} = {'Gaussian', 'Gaussian', 'Gaussian'};
C1C2C3_combinations{2} = {'Frank', 'Gaussian', 'Frank'};

dependency_combinations = cell(1,2);
dependency_combinations{1} = {'Strong', 'Strong', 'Strong'};
dependency_combinations{2} = {'Weak', 'Weak', 'Weak'};

if(exist(logFile,'file')==2)
    % delete the file
    delete(logFile);
end

% decide which function to run based on D and the configuration
switch D
    case 2
        switch cfg
            case 1
                [llMat, llVarMat, llBiasMat, llMCCell] = runD2CFG1();
            case 2
                [llMat, llVarMat, llBiasMat, llMCCell] = runD2CFG2();
            case 3
                [llMat, llVarMat, llBiasMat, llMCCell] = runD2CFG3();
            case 4
                [llMat, llVarMat, llBiasMat, llMCCell] = runD2CFG4();
            case 5      % BIAS/VARIANCE TEST CASE
                [llMat, llVarMat, llBiasMat, llMCCell] = runD2CFG5();
            case 6      % BIAS/VARIANCE TEST CASE
                [llMat, llVarMat, llBiasMat, llMCCell] = runD2CFG6();
            case 7      % Pathological BIAS/VARIANCE TEST CASE
                [llMat, llVarMat, llBiasMat, llMCCell] = runD2CFG7();
            otherwise
                error('CFG not valid for D=2!');
        end
    case 3
        switch cfg
            case 1
                [llMat, llVarMat, llBiasMat, llMCCell] = runD3CFG1();
            case 2
                [llMat, llVarMat, llBiasMat, llMCCell] = runD3CFG2();
            case 3
                [llMat, llVarMat, llBiasMat, llMCCell] = runD3CFG3();
            case 4
                [llMat, llVarMat, llBiasMat, llMCCell] = runD3CFG4();
            case 5      % BIAS/VARIANCE TEST CASE
                [llMat, llVarMat, llBiasMat, llMCCell] = runD3CFG5();
            case 6      % BIAS/VARIANCE TEST CASE
                [llMat, llVarMat, llBiasMat, llMCCell] = runD3CFG6();
            otherwise
                error('CFG not valid for D=3!');
        end
        case 4
            switch cfg
                case 1
                    [llMat, llVarMat, llBiasMat, llMCCell] = runD4CFG1();
                case 2
                    [llMat, llVarMat, llBiasMat, llMCCell] = runD4CFG2();
                case 3
                    [llMat, llVarMat, llBiasMat, llMCCell] = runD4CFG3();
                case 4
                    [llMat, llVarMat, llBiasMat, llMCCell] = runD4CFG4();
                case 5      % BIAS/VARIANCE TEST CASE
                    [llMat, llVarMat, llBiasMat, llMCCell] = runD4CFG5();
                case 6      % BIAS/VARIANCE TEST CASE
                    [llMat, llVarMat, llBiasMat, llMCCell] = runD4CFG6();
                otherwise
                    error('CFG not valid for D=4!');
            end
    case 5
        switch cfg
            case 1
                [llMat, llVarMat, llBiasMat, llMCCell] = runD5CFG1();
            case 2
                [llMat, llVarMat, llBiasMat, llMCCell] = runD5CFG2();
            case 3
                [llMat, llVarMat, llBiasMat, llMCCell] = runD5CFG3();
            case 4
                [llMat, llVarMat, llBiasMat, llMCCell] = runD5CFG4();
            case 5      % BIAS/VARIANCE TEST CASE
                [llMat, llVarMat, llBiasMat, llMCCell] = runD5CFG5();
            case 6      % BIAS/VARIANCE TEST CASE
                [llMat, llVarMat, llBiasMat, llMCCell] = runD5CFG6();
            otherwise
                error('CFG not valid for D=5!');
        end
    otherwise
        error('D != 2,3,5 not yet implemented!');
end

end

function [llMat, llVarMat, llBiasMat, llMCCell] = runD2CFG1()
% runD2CFG2 - the multinomial probabilities are skewed right
probs = [0.5 0.5];
[llMat, llVarMat, llBiasMat, llMCCell] = runD2(probs);
end

function [llMat, llVarMat, llBiasMat, llMCCell] = runD2CFG2()
% runD2CFG2 - the multinomial probabilities are skewed left
probs = [0.7 0.3];
[llMat, llVarMat, llBiasMat, llMCCell] = runD2(probs);
end

function [llMat, llVarMat, llBiasMat, llMCCell] = runD2CFG3()
% runD2CFG1 - the multinomial probabilities are evenly distributed
probs = [0.25 0.25 0.25 0.25];
[llMat, llVarMat, llBiasMat, llMCCell] = runD2(probs);
end

function [llMat, llVarMat, llBiasMat, llMCCell] = runD2CFG4()
% runD2CFG2 - the multinomial probabilities are skewed left
probs = [0.5 0.3 0.1 0.1];
[llMat, llVarMat, llBiasMat, llMCCell] = runD2(probs);
end

function [llMat, llVarMat, llBiasMat, llMCCell] = runD2CFG5()
% runD2CFG2 - the multinomial probabilities are uniformly distributed
numElem = 10;
probs = 1/numElem*ones(1,numElem);
[llMat, llVarMat, llBiasMat, llMCCell] = runD2(probs);
end

function [llMat, llVarMat, llBiasMat, llMCCell] = runD2CFG6()
% runD2CFG2 - the multinomial probabilities are skewed left
probs = [.25 .2 .15 .1 .05 .05 .05 .05 .05 .05];
[llMat, llVarMat, llBiasMat, llMCCell] = runD2(probs);
end

function [llMat, llVarMat, llBiasMat, llMCCell] = runD2CFG7()
% runD2CFG2 - the multinomial probabilities are uniformly distributed
numElem = 20;
probs = 1/numElem*ones(1,numElem);
[llMat, llVarMat, llBiasMat, llMCCell] = runD2(probs);
end

function [llMat, llVarMat, llBiasMat, llMCCell] = runD2(a_probs)
global mVec copulaTypeVec_2D alphaVec RhoVecs_2D continuousDistTypeVec 
global numModelsCompared numMC bntPath logFile K h
global plotFlag numTest
global HCBN_LL_MAT_IDX MTE_LL_MAT_IDX CLG_LL_MAT_IDX REF_LL_MAT_IDX
global MULTINOMIAL_LL_MAT_IDX NUM_DISCRETE_INTERVALS CBN_LL_MAT_IDX
global HCBN_DEBUGALL_LL_MAT_IDX HCBN_DEBUGCOPULA_LL_MAT_IDX 
global CBN_GAUSSIAN_LL_MAT_IDX HCBN_DEBUGEMPINFO_LL_MAT_IDX

a_dist = makedist('Multinomial','Probabilities',a_probs);

% setup the graphical model
aa = 1; bb = 2;
D = 2;
dag = zeros(D,D);
dag(aa,bb) = 1;
discreteNodes = [aa];
nodeNames = {'A', 'B'};
discreteNodeNames = {'A'};

llMCMat = zeros(numModelsCompared,numMC);
llMat = zeros(length(copulaTypeVec_2D),...
              length(continuousDistTypeVec),...
              length(alphaVec),...
              length(mVec),...
              numModelsCompared);
llVarMat = zeros(length(copulaTypeVec_2D),...
              length(continuousDistTypeVec),...
              length(alphaVec),...
              length(mVec),...
              numModelsCompared);
llBiasMat = zeros(length(copulaTypeVec_2D),...
              length(continuousDistTypeVec),...
              length(alphaVec),...
              length(mVec),...
              numModelsCompared);
llMCCell = cell(length(copulaTypeVec_2D),...
              length(continuousDistTypeVec),...
              length(alphaVec),...
              length(mVec),...
              1);
          
fid = fopen(logFile, 'a');
dispstat('','init'); % One time only initialization
dispstat(sprintf('Begining the simulation...\n'),'keepthis','timestamp');
numTotalLoops = length(copulaTypeVec_2D)*length(continuousDistTypeVec)*length(alphaVec)*length(mVec)*numMC;
progressIdx = 1;
for copulaTypeVecIdx=1:length(copulaTypeVec_2D)
    for continuousDistTypeVecIdx=1:length(continuousDistTypeVec)
        for alphaVecIdx=1:length(alphaVec)
            for mVecIdx=1:length(mVec)
                copulaType = copulaTypeVec_2D{copulaTypeVecIdx};
                continuousDistType = continuousDistTypeVec{continuousDistTypeVecIdx};
                if(strcmp(copulaType, 'Gaussian'))
                    alpha = RhoVecs_2D{alphaVecIdx};        % alpha is Rho here
                else
                    alpha = alphaVec(alphaVecIdx);
                end
                M = mVec(mVecIdx);
                
                %%%%%%%%%%% MAIN SIMULATION CODE %%%%%%%%%%
                progressAmt = progressIdx/numTotalLoops*100;
                if(strcmp(copulaType,'Gaussian'))
                    progressStr = sprintf('copulaType=%s x2DistType=%s rho=%f M=%d || Progress=%0.04f', ...
                                    copulaType, continuousDistType, alpha(1,2), M, progressAmt);
                else
                    progressStr = sprintf('copulaType=%s x2DistType=%s alpha=%d M=%d || Progress=%0.04f', ...
                                    copulaType, continuousDistType, alpha, M, progressAmt);
                end
                dispstat(progressStr,'keepthis','timestamp');
                fprintf(fid, progressStr);
                
                % make the actual copula density, and the partial
                % derivative of the copula function w.r.t. the
                % continuous variable (x_2)
                u = linspace(0,1,K);
                [U1,U2] = ndgrid(u);
                
                if(strcmpi(copulaType, 'Gaussian'))
                    cpdf = copulapdf('Gaussian', [U1(:) U2(:)], alpha);
                elseif(strcmpi(copulaType, 'Frank'))
                    cpdf = frankcopulapdf([U1(:) U2(:)], alpha);
                elseif(strcmpi(copulaType, 'Gumbel'))
                    cpdf = gumbelcopulapdf([U1(:) U2(:)], alpha);
                elseif(strcmpi(copulaType, 'Clayton'))
                    cpdf = claytoncopulapdf([U1(:) U2(:)], alpha);
                else
                    error('Unsupported copula type!');
                end
                c_actual = reshape(cpdf, K, K );
                C_actual_discrete_integrate = cumtrapz(u, c_actual, 1);
                
                % define the continuous distribution
                if(strcmp(continuousDistType, 'Multimodal'))
                    xContinuous = [normrnd(-2,0.3,1000,1); normrnd(2,0.8,1000,1)];
                elseif(strcmp(continuousDistType, 'Uniform'))
                    xContinuous = unifrnd(-2,2,2000,1);
                elseif(strcmp(continuousDistType, 'UnimodalSkewed'))
                    xContinuous = betarnd(2,5,2000,1);
                elseif(strcmp(continuousDistType, 'Gaussian'))
                    xContinuous = normrnd(2,0.5,2000,1);
                elseif(strcmp(continuousDistType, 'ThickTailed'))
                    xContinuous = trnd(3, 2000, 1);
                else
                    error('Unknown X2 Dist Type!');
                end
                xContinuous = xContinuous(randperm(2000),:);     % permute for evenness of samples

                isdiscrete = 0;
                [fContinous,xiContinuous] = emppdf(xContinuous,isdiscrete);
                FContinuous = empcdf(xContinuous,isdiscrete);
                continuousDistInfo = rvEmpiricalInfo(xiContinuous,fContinous,FContinuous,isdiscrete);
                
                copulaFamilies = cell(1,2);
                copulaFamilies{1} = [];
                tmp = cell(1,2); 
                tmp{1} = copulaType; tmp{2} = alpha;
                % NOTE: we don't need to do any permuting of the RHO matrix
                % because this is 2-D.  Rho(1,2) in the original matrix is
                % the correlation between parent and child.  Because
                % correlation is symmetric, the matrix does not need to be
                % permuted
                copulaFamilies{2} = tmp;
                empInfo = cell(1,2);
                empInfo{1} = a_dist; empInfo{2} = continuousDistInfo;
                
                for mcSimNum=1:numMC
                    dispstat(sprintf('MC Sim=%d', mcSimNum), 'timestamp');
                    U = copularnd(copulaType, alpha, M+numTest);
                    X_hybrid = zeros(M+numTest,2);

                    X_hybrid(:,1) = a_dist.icdf(U(:,1));
                    for ii=1:M+numTest
                        X_hybrid(ii,2) = continuousDistInfo.icdf(U(ii,2));
                    end
                    % generate train and test datasets
                    X_hybrid_test = X_hybrid(M+1:end,:);
                    X_hybrid = X_hybrid(1:M,:);
                    
                    isdiscrete = 1;
                    [fAest,xAest] = emppdf(X_hybrid(:,1),isdiscrete);
                    FAest = empcdf(X_hybrid(:,1),isdiscrete);
                    distAEst = rvEmpiricalInfo(xAest,fAest,FAest,isdiscrete);

                    isdiscrete = 0;
                    [fBest,xBest] = emppdf(X_hybrid(:,2),isdiscrete); FBest = empcdf(X_hybrid(:,2),isdiscrete);
                    distBEst = rvEmpiricalInfo(xBest,fBest,FBest,isdiscrete);

                    hcbnObj = hcbn(bntPath, X_hybrid, nodeNames, discreteNodeNames, K, h, dag); 
                    hcbnDebugAllObj = hcbn(bntPath, X_hybrid, nodeNames, discreteNodeNames, K, h, dag, 'copulaFamilyInput', copulaFamilies, 'empInfoInput', empInfo); 
                    hcbnDebugCopulaObj = hcbn(bntPath, X_hybrid, nodeNames, discreteNodeNames, K, h, dag, 'copulaFamilyInput', copulaFamilies); 
                    hcbnDebugEmpInfoObj = hcbn(bntPath, X_hybrid, nodeNames, discreteNodeNames, K, h, dag, 'empInfoInput', empInfo); 
                    mtebnObj = mtebn(X_hybrid, discreteNodes, dag);
                    clgbnObj = clgbn(X_hybrid, discreteNodes, dag);
                    multinomialbnObj = multinomialbn(X_hybrid, discreteNodes, dag, NUM_DISCRETE_INTERVALS);
                    cbnObj = cbn(bntPath, X_hybrid, nodeNames, dag);
                    cbnGaussianObj = cbn(bntPath, X_hybrid, nodeNames, dag, 1);

                    X_hybrid_continued = X_hybrid;
                    X_hybrid_continued(:,1) = continueRv(X_hybrid(:,1));
                    % generate pseudo-observations
                    U_hybrid_continued = pobs(X_hybrid_continued, 'ecdf', 100);
                    c_est = empcopulapdf(U_hybrid_continued, h, K, 'betak');
                    C_est_discrete_integrate = cumtrapz(u, c_est, 1);
                    
                    if(plotFlag)
                        % plot some conditional distributions from the known copula and known
                        % marginal distributions
                        x1_discrete_conditional_vec = 1:1:length(a_probs);
                        for x1_discrete_conditional=x1_discrete_conditional_vec
                            X_continuous_subset = [];
                            for jj=1:M
                                if(X_hybrid(jj,1)==x1_discrete_conditional)
                                    X_continuous_subset = [X_continuous_subset; X_hybrid(jj,2)];
                                end
                            end

                            % estimate the empirical PDF/CDF for this
                            isdiscrete = 0;
                            [f_kde,xi_kde] = emppdf(X_continuous_subset, isdiscrete);
                            F_kde = empcdf(X_continuous_subset, isdiscrete);
                            conditionalKDE = rvEmpiricalInfo(xi_kde, f_kde, F_kde,isdiscrete);

                            fx2_givenx1_generative = zeros(1,length(xiContinuous));
                            fx2_givenx1_copulaestf2est = zeros(1,length(xiContinuous));
                            fx2_givenx1_copulaestf2Actual = zeros(1,length(xiContinuous));
                            fx2_givenx1_copulaActualf2est = zeros(1,length(xiContinuous));
                            fx2_givenx1_copulahcbnf2hcbn = zeros(1,length(xiContinuous));
                            fx2_givenx1_copulahcbnf2Actual = zeros(1,length(xiContinuous));

                            fx2_givenx1_clg = zeros(1,length(xiContinuous));
                            fx2_givenx1_mte = zeros(1,length(xiContinuous));
                            fx2_givenx1_conditionalKDE = zeros(1,length(xiContinuous));

                            copulaActual = zeros(1,length(xiContinuous));
                            copulaEst = zeros(1,length(xiContinuous));
                            copulaHcbn = zeros(1,length(xiContinuous));
                            f2Actual = zeros(1,length(xiContinuous));
                            f2Est = zeros(1,length(xiContinuous));
                            f2Hcbn = zeros(1,length(xiContinuous));
                            u2Val = zeros(1,length(xiContinuous));

                            % compute the conditional density w/ the formula
                            % given by Eq.9 in Joint Regression 
                            % Analysis of Correlated Data Using Gaussian 
                            % Copulas -- Biometrics 2009
                            % Authors: Song, Li, Yuan
                            fX1 = a_dist.pdf(x1_discrete_conditional);
                            fX1_hat = distAEst.pdf(x1_discrete_conditional);
                            for ii=1:length(xiContinuous)
                                xiContinuous_val = xiContinuous(ii);
                                uuGenerative1 = [a_dist.cdf(x1_discrete_conditional) continuousDistInfo.cdf(xiContinuous_val)];
                                uuGenerative2 = [a_dist.cdf(x1_discrete_conditional-1) continuousDistInfo.cdf(xiContinuous_val)];

                                uuEst1 = [distAEst.cdf(x1_discrete_conditional) ...
                                          distBEst.cdf(xiContinuous_val)];
                                uuEst2 = [distAEst.cdf(x1_discrete_conditional-1) ...
                                          distBEst.cdf(xiContinuous_val)];
                                C_partial_X1X2 = empcopulaval(C_actual_discrete_integrate, uuGenerative1,1/K)-empcopulaval(C_actual_discrete_integrate, uuGenerative2,1/K);
                                C_partial_hat_X1X2 = empcopulaval(C_est_discrete_integrate, uuEst1,1/K) - empcopulaval(C_est_discrete_integrate, uuEst2,1/K);
                                C_partial_hcbn = empcopulaval(hcbnObj.copulaFamilies{2}.C_discrete_integrate, fliplr(uuEst1),1/K) - empcopulaval(hcbnObj.copulaFamilies{2}.C_discrete_integrate, fliplr(uuEst2),1/K);

                                fX2 = continuousDistInfo.pdf(xiContinuous_val);
                                fX2_hat = distBEst.pdf(xiContinuous_val);                                                                            % was estimated w/ [u2 u1]
                                f_x2_hcbn = hcbnObj.empInfo{2}.pdf(xiContinuous_val);

                                fx2_givenx1_generative(ii) = C_partial_X1X2*fX2/fX1;
                                fx2_givenx1_copulaestf2est(ii) = C_partial_hat_X1X2*fX2_hat/fX1_hat;
                                fx2_givenx1_copulaestf2Actual(ii) = C_partial_hat_X1X2*fX2/fX1;
                                fx2_givenx1_copulaActualf2est(ii) = C_partial_X1X2*fX2_hat/fX1_hat;
                                fx2_givenx1_copulahcbnf2hcbn(ii) = hcbnObj.computeMixedConditionalProbability_(...
                                    [x1_discrete_conditional xiContinuous_val], [bb aa], bb);
                                fx2_givenx1_copulahcbnf2Actual(ii) = C_partial_hcbn*fX2/fX1;

                                fx2_givenx1_clg(ii) = normpdf(xiContinuous_val, clgbnObj.bnParams{2}{x1_discrete_conditional}.Mean, clgbnObj.bnParams{2}{x1_discrete_conditional}.Covariance);
                                fx2_givenx1_mte(ii) = mtebnObj.bnParams{2}{x1_discrete_conditional}.mte_info.pdf(xiContinuous_val);
                                fx2_givenx1_conditionalKDE(ii) = conditionalKDE.pdf(xiContinuous_val);

                                copulaActual(ii) = C_partial_X1X2;
                                copulaEst(ii) = C_partial_hat_X1X2;
                                copulaHcbn(ii) = C_partial_hcbn;
                                f2Actual(ii) = fX2;
                                f2Est(ii) = fX2_hat;
                                f2Hcbn(ii) = f_x2_hcbn;
                                u2Val(ii) = uuGenerative1(2);
                            end

                            % plot the actual vs the copula version and compare the differences
                            fig1 = figure(1); subplot(2,2,1);
                            plot(xiContinuous, fx2_givenx1_generative, 'b*-', ...
                                 xiContinuous, fx2_givenx1_copulaestf2est, ...
                                 xiContinuous, fx2_givenx1_copulaestf2Actual, ...
                                 xiContinuous, fx2_givenx1_copulaActualf2est, ...
                                 xiContinuous, fx2_givenx1_copulahcbnf2hcbn, 'o-', ...
                                 xiContinuous, fx2_givenx1_copulahcbnf2Actual); 
                            grid on; title(sprintf('X_1=%d',x1_discrete_conditional));
                            h_legend1 = legend('$c*f(x_2)$', ...
                                '$\hat{c}_{ECDF}*\hat{f}(x_2)$', ...
                                '$\hat{c}_{ECDF}*f(x_2)$', ...
                                '$c*\hat{f}(x_2)$', ...
                                '$\hat{c}_{HCBN-RANK}*\hat{f}(x_2)_{HCBN}$', ...
                                '$\hat{c}_{HCBN-RANK}*f(x_2)$');
                            set(h_legend1,'FontSize',12);
                            set(h_legend1,'Interpreter','latex')

                            subplot(2,2,2);
                            plot(xiContinuous, fx2_givenx1_generative, 'b*-', ...
                                 xiContinuous, fx2_givenx1_conditionalKDE, ...
                                 xiContinuous, fx2_givenx1_mte, ...
                                 xiContinuous, fx2_givenx1_clg); 
                            grid on; title(sprintf('X_1=%d EstimationSize=%d',x1_discrete_conditional, length(X_continuous_subset)));
                            h_legend2 = legend('$c*f(x_2)$', 'KDE', 'MTE', 'CLG');
                            set(h_legend2,'FontSize',12);
                            set(h_legend2,'Interpreter','latex')

                            subplot(2,2,3);
                            plot(u2Val, copulaActual, ...
                                 u2Val, copulaEst, ...
                                 u2Val, copulaHcbn);
                            grid on;
                            if(strcmpi(copulaType,'Gaussian'))
                                alphaDisp = alpha(1,2);
                            else
                                alphaDisp = alpha;
                            end
                            title(sprintf('%s(%0.02f) Marginal Copula u_1=%0.02f', copulaType, alphaDisp, uuGenerative1(1)));
                            h_legend3 = legend('$c$', '$\hat{c}_{ECDF}$', '$\hat{c}_{HCBN-RANK}$');
                            xlabel('u_2');
                            set(h_legend3,'FontSize',12);
                            set(h_legend3,'Interpreter','latex')

                            subplot(2,2,4);
                            plot(xiContinuous, f2Actual, ...
                                 xiContinuous, f2Est, 'k',...
                                 xiContinuous, f2Hcbn, 'rp');
                            grid on;
                            h_legend4 = legend('$f(x_2)$', '$\hat{f}(x_2)$', '$\hat{f}(x_2)_{HCBN}$');
                            set(h_legend4,'FontSize',12);
                            set(h_legend4,'Interpreter','latex')
                            xlabel('x_2');
                            title(sprintf('%s M=%d\n', continuousDistType, M));

                            pause;
                            clf(fig1);
                        end
                    end
                    
                    % calculate the reference LL
                    refLL = 0;
                    totalProbVecRef = zeros(1,numTest);
                    for ii=1:numTest
                        xx = X_hybrid_test(ii,:);
                        uu_continuous = continuousDistInfo.cdf(xx(2));
                        uu1 = [a_dist.cdf(xx(1)) uu_continuous];
                        uu2 = [a_dist.cdf(xx(1)-1) uu_continuous];
                        totalProb = continuousDistInfo.pdf(xx(2)) * ...
                            (empcopulaval(C_actual_discrete_integrate, uu1,1/K) - empcopulaval(C_actual_discrete_integrate, uu2,1/K));
                        if(totalProb<1e-5)
                            totalProb = 1e-5;   % PUT BREAK POINT HERE IF YOU WANT TO DEBUG
                        end
                        refLL = refLL + log(totalProb);
                        totalProbVecRef(ii) = totalProb;
                    end
                    
                    % calculate LL values and assign to llDivMCMat
                    hcbnLL = hcbnObj.dataLogLikelihood(X_hybrid_test);
                    hcbnDebugAllLL = hcbnDebugAllObj.dataLogLikelihood(X_hybrid_test);
                    hcbnDebugCopulaLL = hcbnDebugCopulaObj.dataLogLikelihood(X_hybrid_test);
                    hcbnDebugEmpInfoLL = hcbnDebugEmpInfoObj.dataLogLikelihood(X_hybrid_test);
                    mteLL = mtebnObj.dataLogLikelihood(X_hybrid_test);
                    clgLL = clgbnObj.dataLogLikelihood(X_hybrid_test);
                    multinomialLL = multinomialbnObj.dataLogLikelihood(X_hybrid_test);
                    cbnLL = cbnObj.dataLogLikelihood(X_hybrid_test);
                    cbnGaussianLL = cbnGaussianObj.dataLogLikelihood(X_hybrid_test);
                    
                    % assign LL values to matrix
                    llMCMat(HCBN_LL_MAT_IDX,mcSimNum) = hcbnLL;
                    llMCMat(HCBN_DEBUGALL_LL_MAT_IDX,mcSimNum) = hcbnDebugAllLL;
                    llMCMat(HCBN_DEBUGCOPULA_LL_MAT_IDX,mcSimNum) = hcbnDebugCopulaLL;
                    llMCMat(HCBN_DEBUGEMPINFO_LL_MAT_IDX,mcSimNum) = hcbnDebugEmpInfoLL;
                    llMCMat(MTE_LL_MAT_IDX,mcSimNum) = mteLL;
                    llMCMat(CLG_LL_MAT_IDX,mcSimNum) = clgLL;
                    llMCMat(MULTINOMIAL_LL_MAT_IDX,mcSimNum) = multinomialLL;
                    llMCMat(REF_LL_MAT_IDX,mcSimNum) = refLL;
                    llMCMat(CBN_LL_MAT_IDX,mcSimNum) = cbnLL;
                    llMCMat(CBN_GAUSSIAN_LL_MAT_IDX,mcSimNum) = cbnGaussianLL;
                    
                    progressIdx = progressIdx + 1;
                end
                meanLLDivMCMat = mean(llMCMat,2);
                varLLDivMCMat  = mean(llMCMat.^2,2)-meanLLDivMCMat.^2;
                biasLLDivMCMat = mean( repmat(llMCMat(REF_LL_MAT_IDX,:), numModelsCompared, 1) - llMCMat, 2 );
                
                progressStr = sprintf('MEAN{ref}=%f VAR{ref}=%f BIAS{ref}=%f', ...
                    meanLLDivMCMat(REF_LL_MAT_IDX), varLLDivMCMat(REF_LL_MAT_IDX), biasLLDivMCMat(REF_LL_MAT_IDX));
                dispstat(progressStr,'timestamp','keepthis','timestamp');
                progressStr = sprintf('MEAN{hcbn}=%f VAR{hcbn}=%f BIAS{hcbn}=%f', ...
                    meanLLDivMCMat(HCBN_LL_MAT_IDX), varLLDivMCMat(HCBN_LL_MAT_IDX), biasLLDivMCMat(HCBN_LL_MAT_IDX));
                dispstat(progressStr,'timestamp','keepthis','timestamp');
                progressStr = sprintf('MEAN{hcbnDebugCopula}=%f VAR{hcbnDebugCopula}=%f BIAS{hcbnDebugCopula}=%f', ...
                    meanLLDivMCMat(HCBN_DEBUGCOPULA_LL_MAT_IDX), varLLDivMCMat(HCBN_DEBUGCOPULA_LL_MAT_IDX), biasLLDivMCMat(HCBN_DEBUGCOPULA_LL_MAT_IDX));
                dispstat(progressStr,'timestamp','keepthis','timestamp');
                progressStr = sprintf('MEAN{hcbnDebugEmpInfo}=%f VAR{hcbnDebugEmpInfo}=%f BIAS{hcbnDebugEmpInfo}=%f', ...
                    meanLLDivMCMat(HCBN_DEBUGEMPINFO_LL_MAT_IDX), varLLDivMCMat(HCBN_DEBUGEMPINFO_LL_MAT_IDX), biasLLDivMCMat(HCBN_DEBUGEMPINFO_LL_MAT_IDX));
                dispstat(progressStr,'timestamp','keepthis','timestamp');
                progressStr = sprintf('MEAN{cbn}=%f VAR{cbn}=%f BIAS{cbn}=%f', ...
                    meanLLDivMCMat(CBN_LL_MAT_IDX), varLLDivMCMat(CBN_LL_MAT_IDX), biasLLDivMCMat(CBN_LL_MAT_IDX));
                dispstat(progressStr,'timestamp','keepthis','timestamp');
                progressStr = sprintf('MEAN{cbnGaussian}=%f VAR{cbnGaussian}=%f BIAS{cbnGaussian}=%f', ...
                    meanLLDivMCMat(CBN_GAUSSIAN_LL_MAT_IDX), varLLDivMCMat(CBN_GAUSSIAN_LL_MAT_IDX), biasLLDivMCMat(CBN_GAUSSIAN_LL_MAT_IDX));
                dispstat(progressStr,'timestamp','keepthis','timestamp');
                progressStr = sprintf('MEAN{mte}=%f VAR{mte}=%f BIAS{mte}=%f', ...
                    meanLLDivMCMat(MTE_LL_MAT_IDX), varLLDivMCMat(MTE_LL_MAT_IDX), biasLLDivMCMat(MTE_LL_MAT_IDX));
                dispstat(progressStr,'timestamp','keepthis','timestamp');
                progressStr = sprintf('MEAN{clg}=%f VAR{clg}=%f BIAS{clg}=%f', ...
                    meanLLDivMCMat(CLG_LL_MAT_IDX), varLLDivMCMat(CLG_LL_MAT_IDX), biasLLDivMCMat(CLG_LL_MAT_IDX));
                dispstat(progressStr,'timestamp','keepthis','timestamp');
                progressStr = sprintf('MEAN{multinomial}=%f VAR{multinomial}=%f BIAS{multinomial}=%f', ...
                    meanLLDivMCMat(MULTINOMIAL_LL_MAT_IDX), varLLDivMCMat(MULTINOMIAL_LL_MAT_IDX), biasLLDivMCMat(MULTINOMIAL_LL_MAT_IDX));
                dispstat(progressStr,'timestamp','keepthis','timestamp');
                progressStr = sprintf('mean{hcbnDebug}==mean{ref}=%d\n', ...
                    abs(meanLLDivMCMat(HCBN_DEBUGALL_LL_MAT_IDX)-meanLLDivMCMat(REF_LL_MAT_IDX)) < .1 );
                dispstat(progressStr,'timestamp','keepthis','timestamp');
                
                llMat(copulaTypeVecIdx, ...
                        continuousDistTypeVecIdx,...
                        alphaVecIdx,...
                        mVecIdx,:) = meanLLDivMCMat(:);
                llVarMat(copulaTypeVecIdx, ...
                         continuousDistTypeVecIdx,...
                         alphaVecIdx,...
                         mVecIdx,:) = varLLDivMCMat(:);
                llBiasMat(copulaTypeVecIdx, ...
                          continuousDistTypeVecIdx,...
                          alphaVecIdx,...
                          mVecIdx,:) = biasLLDivMCMat(:);  
                llMCCell{copulaTypeVecIdx, ...
                          continuousDistTypeVecIdx,...
                          alphaVecIdx,...
                          mVecIdx,1} = llMCMat;
                %%%%%%%%%%% END OF MAIN SIMULATION CODE %%%%%%%%%%
            end
        end
    end
end

dispstat('Finished.','keepprev');
fclose(fid);

end

function [llMat, llVarMat, llBiasMat, llMCCell] = runD3CFG1()
% runD2CFG2 - the multinomial probabilities are skewed right
probsA = [0.5 0.5];
probsB = [0.5 0.5];
[llMat, llVarMat, llBiasMat, llMCCell] = runD3(probsA, probsB);
end

function [llMat, llVarMat, llBiasMat, llMCCell] = runD3CFG2()
% runD2CFG2 - the multinomial probabilities are skewed right and left
% oppositely
probsA = [0.3 0.7];
probsB = fliplr([0.3 0.7]);
[llMat, llVarMat, llBiasMat, llMCCell] = runD3(probsA, probsB);
end

function [llMat, llVarMat, llBiasMat, llMCCell] = runD3CFG3()
% runD2CFG1 - the multinomial probabilities are evenly distributed
probsA = [0.25 0.25 0.25 0.25];
probsB = [0.25 0.25 0.25 0.25];
[llMat, llVarMat, llBiasMat, llMCCell] = runD3(probsA, probsB);
end

function [llMat, llVarMat, llBiasMat, llMCCell] = runD3CFG4()
% runD2CFG2 - the multinomial probabilities are skewed left and right
% oppositely
probsA = [0.5 0.3 0.1 0.1];
probsB = fliplr([0.5 0.3 0.1 0.1]);
[llMat, llVarMat, llBiasMat, llMCCell] = runD3(probsA, probsB);
end

function [llMat, llVarMat, llBiasMat, llMCCell] = runD3CFG5()
% runD2CFG1 - the multinomial probabilities are evenly distributed
probsA = .1*ones(1,10);
probsB = .1*ones(1,10);
[llMat, llVarMat, llBiasMat, llMCCell] = runD3(probsA, probsB);
end

function [llMat, llVarMat, llBiasMat, llMCCell] = runD3CFG6()
% runD2CFG1 - both multinomial probabilities are skewed left & right
probsA = [.25 .2 .15 .1 .05 .05 .05 .05 .05 .05];
probsB = fliplr([.25 .2 .15 .1 .05 .05 .05 .05 .05 .05]);
[llMat, llVarMat, llBiasMat, llMCCell] = runD3(probsA, probsB);
end

function [llMat, llVarMat, llBiasMat, llMCCell] = runD3(a_probs, b_probs)
global mVec copulaTypeVec_3D alphaVec RhoVecs_3D continuousDistTypeVec 
global numModelsCompared numMC bntPath logFile K h plotFlag numTest
global HCBN_LL_MAT_IDX MTE_LL_MAT_IDX CLG_LL_MAT_IDX REF_LL_MAT_IDX
global MULTINOMIAL_LL_MAT_IDX NUM_DISCRETE_INTERVALS CBN_LL_MAT_IDX
global HCBN_DEBUGALL_LL_MAT_IDX HCBN_DEBUGCOPULA_LL_MAT_IDX HCBN_DEBUGEMPINFO_LL_MAT_IDX
global CBN_GAUSSIAN_LL_MAT_IDX

a_dist = makedist('Multinomial','Probabilities', a_probs);
b_dist = makedist('Multinomial','Probabilities', b_probs);

% setup the graphical model
aa = 1; bb = 2; cc = 3;
D = 3;
dag = zeros(D,D);
dag(aa,cc) = 1;
dag(bb,cc) = 1;
discreteNodes = [aa bb];
nodeNames = {'A', 'B', 'C'};
discreteNodeNames = {'A','B'};

llMCMat = zeros(numModelsCompared,numMC);
llMat = zeros(length(copulaTypeVec_3D),...
              length(continuousDistTypeVec),...
              length(alphaVec),...
              length(mVec),...
              numModelsCompared);
llVarMat = zeros(length(copulaTypeVec_3D),...
              length(continuousDistTypeVec),...
              length(alphaVec),...
              length(mVec),...
              numModelsCompared);
llBiasMat = zeros(length(copulaTypeVec_3D),...
              length(continuousDistTypeVec),...
              length(alphaVec),...
              length(mVec),...
              numModelsCompared);
llMCCell = cell(length(copulaTypeVec_3D),...
              length(continuousDistTypeVec),...
              length(alphaVec),...
              length(mVec),...
              1);
          
fid = fopen(logFile, 'a');
dispstat('','init'); % One time only initialization
dispstat(sprintf('Begining the simulation...'),'keepthis','timestamp');
numTotalLoops = length(copulaTypeVec_3D)*length(continuousDistTypeVec)*length(alphaVec)*length(mVec)*numMC;
progressIdx = 1;
for copulaTypeVecIdx=1:length(copulaTypeVec_3D)
    for continuousDistTypeVecIdx=1:length(continuousDistTypeVec)
        for alphaVecIdx=1:length(alphaVec)
            for mVecIdx=1:length(mVec)
                copulaType = copulaTypeVec_3D{copulaTypeVecIdx};
                continuousDistType = continuousDistTypeVec{continuousDistTypeVecIdx};
                if(strcmp(copulaType, 'Gaussian'))
                    Rho = RhoVecs_3D{alphaVecIdx};        % alpha is Rho here
                else
                    error('Unrecognized Copula Type!');
                end
                M = mVec(mVecIdx);
                
                progressAmt = progressIdx/numTotalLoops*100;
                if(strcmp(copulaType,'Gaussian'))
                    progressStr = sprintf('copulaType=%s x3DistType=%s rho=%f M=%d || Progress=%0.04f', ...
                                    copulaType, continuousDistType, Rho(1,3), M, progressAmt);
                else
                    error('Unrecognized Copula Type!');
                end
                dispstat(progressStr,'keepthis','timestamp');
                fprintf(fid, progressStr);
                
                if(strcmp(continuousDistType, 'Multimodal'))
                    xContinuous = [normrnd(-2,0.3,1000,1); normrnd(2,0.8,1000,1)];
                elseif(strcmp(continuousDistType, 'Uniform'))
                    xContinuous = unifrnd(-2,2,2000,1);
                elseif(strcmp(continuousDistType, 'UnimodalSkewed'))
                    xContinuous = betarnd(2,5,2000,1);
                elseif(strcmp(continuousDistType, 'Gaussian'))
                    xContinuous = normrnd(2,0.5,2000,1);
                elseif(strcmp(continuousDistType, 'ThickTailed'))
                    xContinuous = trnd(3, 2000, 1);
                else
                    error('Unknown X3 Dist Type!');
                end
                xContinuous = xContinuous(randperm(2000),:);     % permute for evenness of samples

                isdiscrete = 0;
                [fContinous,xiContinuous] = emppdf(xContinuous,isdiscrete);
                FContinuous = empcdf(xContinuous,isdiscrete);
                continuousDistInfo = rvEmpiricalInfo(xiContinuous,fContinous,FContinuous,isdiscrete);
                
                copulaFamilies = cell(1,3);
                copulaFamilies{1} = [];
                copulaFamilies{2} = [];
                tmp = cell(1,2); 
                tmp{1} = copulaType; 
                if(strcmpi(copulaType, 'Gaussian'))
                    % permute the correlation matrix such that it is in the
                    % format of Child,Parent1,Parent2
                    tmp{2} = circshift(circshift(Rho,1,1),1,2);
                else
                    error('Unrecognized Copula Type!!');
                end
                copulaFamilies{3} = tmp;
                empInfo = cell(1,3);
                empInfo{1} = a_dist; empInfo{2} = b_dist;
                empInfo{3} = continuousDistInfo;
                
                u = linspace(0,1,K); [U1_3,U2_3,U3_3] = ndgrid(u); [U1_2,U2_2] = ndgrid(u);
                if(strcmp(copulaType, 'Gaussian'))
                    c_actual_X1X2X3 = reshape(copulapdf('Gaussian', [U1_3(:) U2_3(:) U3_3(:)], Rho),K,K,K);
                    c_actual_X1X2 = reshape(copulapdf('Gaussian', [U1_2(:) U2_2(:)], Rho(1:2,1:2)),K,K);
                else
                    error('Copula Type not recognized!\n');
                end
                C_actual_X1X2X3_discrete_integrate = cumtrapz(u,cumtrapz(u, c_actual_X1X2X3, 1),2);
                C_actual_X1X2_discrete_integrate = cumtrapz(u,cumtrapz(u, c_actual_X1X2, 1),2);
                
                %%%%%%%%%%% MAIN SIMULATION CODE %%%%%%%%%%
                for mcSimNum=1:numMC
                    dispstat(sprintf('MC Sim=%d', mcSimNum), 'timestamp');
                    X_hybrid = zeros(M+numTest,D);
                    
                    % generate the copula random variates
                    if(strcmp(copulaType, 'Gaussian'))
                        u_X1X2X3 = copularnd('Gaussian', Rho, M+numTest);
                    else
                        error('Copula Type not recognized!\n');
                    end
                    
                    X_hybrid(:,1) = a_dist.icdf(u_X1X2X3(:,1));
                    X_hybrid(:,2) = b_dist.icdf(u_X1X2X3(:,2));
                    for ii=1:M+numTest
                        X_hybrid(ii,3) = continuousDistInfo.icdf(u_X1X2X3(ii,3));
                    end
                    % generate train and test datasets
                    X_hybrid_test = X_hybrid(M+1:end,:);
                    X_hybrid = X_hybrid(1:M,:);
                    
                    isdiscrete = 1;
                    [fAest,xAest] = emppdf(X_hybrid(:,1),isdiscrete);
                    FAest = empcdf(X_hybrid(:,1),isdiscrete);
                    distAEst = rvEmpiricalInfo(xAest,fAest,FAest,isdiscrete);
                    
                    [fBest,xBest] = emppdf(X_hybrid(:,2),isdiscrete);
                    FBest = empcdf(X_hybrid(:,2),isdiscrete);
                    distBEst = rvEmpiricalInfo(xBest,fBest,FBest,isdiscrete);
                    
                    isdiscrete = 0;
                    [fCest,xCest] = emppdf(X_hybrid(:,3),isdiscrete); 
                    FCest = empcdf(X_hybrid(:,3),isdiscrete);
                    distCEst = rvEmpiricalInfo(xCest,fCest,FCest,isdiscrete);
                    
                    hcbnObj = hcbn(bntPath, X_hybrid, nodeNames, discreteNodeNames, K, h, dag);  
                    hcbnDebugAllObj = hcbn(bntPath, X_hybrid, nodeNames, discreteNodeNames, K, h, dag, 'copulaFamilyInput', copulaFamilies, 'empInfoInput', empInfo);  
                    hcbnDebugCopulaObj = hcbn(bntPath, X_hybrid, nodeNames, discreteNodeNames, K, h, dag, 'copulaFamilyInput', copulaFamilies); 
                    hcbnDebugEmpInfoObj = hcbn(bntPath, X_hybrid, nodeNames, discreteNodeNames, K, h, dag, 'empInfoInput', empInfo); 
                    mtebnObj = mtebn(X_hybrid, discreteNodes, dag);
                    clgbnObj = clgbn(X_hybrid, discreteNodes, dag);
                    multinomialbnObj = multinomialbn(X_hybrid, discreteNodes, dag, NUM_DISCRETE_INTERVALS);
                    cbnObj = cbn(bntPath, X_hybrid, nodeNames, dag);
                    cbnGaussianObj = cbn(bntPath, X_hybrid, nodeNames, dag, 1);
                    
                    X_hybrid_continued = X_hybrid;
                    X_hybrid_continued(:,1) = continueRv(X_hybrid(:,1));
                    X_hybrid_continued(:,2) = continueRv(X_hybrid(:,2));
                    
                    U_hybrid_continued = pobs(X_hybrid_continued);
                    
                    % setup all the copula calculations for querying after
                    c_est_X1X2X3 = empcopulapdf(U_hybrid_continued, h, K, 'betak');
                    C_est_X1X2X3 = cumtrapz(u,cumtrapz(u,cumtrapz(u,c_est_X1X2X3,1),2),3);
                    C_est_X1X2X3_discrete_integrate = cumtrapz(u,cumtrapz(u,c_est_X1X2X3,1),2);
                    c_est_X1X2 = empcopulapdf(U_hybrid_continued(:,1:2), h, K, 'betak');
                    C_est_X1X2 = cumtrapz(u,cumtrapz(u,c_est_X1X2,1),2);
                    C_est_X1X2_discrete_integrate = C_est_X1X2;
                    %%%%%%%%%%%%%%%%%%%%%%%
                    
                    if(plotFlag)
                        x1_discrete_conditional_vec = 1:1:length(a_probs);
                        x2_discrete_conditional_vec = 1:1:length(b_probs);
                        for x1_discrete_conditional=x1_discrete_conditional_vec
                            for x2_discrete_conditional=x2_discrete_conditional_vec

                                % generate the continuous subset
                                X_continuous_subset = [];
                                for jj=1:M
                                    if(X_hybrid(jj,1)==x1_discrete_conditional && X_hybrid(jj,2)==x2_discrete_conditional)
                                        X_continuous_subset = [X_continuous_subset; X_hybrid(jj,3)];
                                    end
                                end

                                % estimate the empirical PDF/CDF for this
                                isdiscrete = 0;
                                [f_kde,xi_kde] = emppdf(X_continuous_subset, isdiscrete);
                                F_kde = empcdf(X_continuous_subset, isdiscrete);
                                conditionalKDE = rvEmpiricalInfo(xi_kde, f_kde, F_kde,isdiscrete);

                                fx3_givenx1x2_copula = zeros(1,length(xiContinuous));
                                fx3_givenx1x2_copulaestf3est = zeros(1,length(xiContinuous));
                                fx3_givenx1x2_copulaestf3Actual = zeros(1,length(xiContinuous));
                                fx3_givenx1x2_copulaActualf3est = zeros(1,length(xiContinuous));
                                fx3_givenx1x2_hcbn = zeros(1,length(xiContinuous));
                                fx3_givenx1x2_clg = zeros(1,length(xiContinuous));
                                fx3_givenx1x2_mte = zeros(1,length(xiContinuous));
                                fx3_givenx1x2_conditionalKDE = zeros(1,length(xiContinuous));

                                % compute the conditional density w/ the formula
                                % given by Eq.9 in Joint Regression 
                                % Analysis of Correlated Data Using Gaussian 
                                % Copulas -- Biometrics 2009
                                % Authors: Song, Li, Yuan
                                for ii=1:length(xiContinuous)
                                    xiContinuous_val = xiContinuous(ii);
                                    uuGenerativeX1X2X3_1 = [a_dist.cdf(x1_discrete_conditional) ...
                                                            b_dist.cdf(x2_discrete_conditional) ...
                                                            continuousDistInfo.cdf(xiContinuous_val)];
                                    uuGenerativeX1X2X3_2 = [a_dist.cdf(x1_discrete_conditional-1) ...
                                                            b_dist.cdf(x2_discrete_conditional) ...
                                                            continuousDistInfo.cdf(xiContinuous_val)];
                                    uuGenerativeX1X2X3_3 = [a_dist.cdf(x1_discrete_conditional) ...
                                                            b_dist.cdf(x2_discrete_conditional-1) ...
                                                            continuousDistInfo.cdf(xiContinuous_val)];
                                    uuGenerativeX1X2X3_4 = [a_dist.cdf(x1_discrete_conditional-1) ...
                                                            b_dist.cdf(x2_discrete_conditional-1) ...
                                                            continuousDistInfo.cdf(xiContinuous_val)];
                                    uuGenerativeX1X2_1 = uuGenerativeX1X2X3_1(1:2);
                                    uuGenerativeX1X2_2 = uuGenerativeX1X2X3_2(1:2);
                                    uuGenerativeX1X2_3 = uuGenerativeX1X2X3_3(1:2);
                                    uuGenerativeX1X2_4 = uuGenerativeX1X2X3_4(1:2);

                                    uuEstX1X2X3_1 = [distAEst.cdf(x1_discrete_conditional) ...
                                                     distBEst.cdf(x2_discrete_conditional) ...
                                                     distCEst.cdf(xiContinuous_val)];
                                    uuEstX1X2X3_2 = [distAEst.cdf(x1_discrete_conditional-1) ...
                                                     distBEst.cdf(x2_discrete_conditional) ...
                                                     distCEst.cdf(xiContinuous_val)];
                                    uuEstX1X2X3_3 = [distAEst.cdf(x1_discrete_conditional) ...
                                                     distBEst.cdf(x2_discrete_conditional-1) ...
                                                     distCEst.cdf(xiContinuous_val)];
                                    uuEstX1X2X3_4 = [distAEst.cdf(x1_discrete_conditional-1) ...
                                                     distBEst.cdf(x2_discrete_conditional-1) ...
                                                     distCEst.cdf(xiContinuous_val)];
                                    uuEstX1X2_1 = uuEstX1X2X3_1(1:2);
                                    uuEstX1X2_2 = uuEstX1X2X3_2(1:2);
                                    uuEstX1X2_3 = uuEstX1X2X3_3(1:2);
                                    uuEstX1X2_4 = uuEstX1X2X3_4(1:2);

    % % %                                 uuHcbnX1X2X3 = [hcbnObj.empInfo{1}.distribution(x1_discrete_conditional) ...
    % % %                                                 hcbnObj.empInfo{2}.distribution(x2_discrete_conditional) ...
    % % %                                                 hcbnObj.empInfo{3}.cdf(xiContinuous_val)];
    % % %                                 uuHcbnX1X2X3 = uuHcbnX1X2X3;
    % % %                                 uuHcbnX1X2 = uuHcbnX1X2X3(1:2);

                                    C_actual_partial_X1X2X3 = empcopulaval(C_actual_X1X2X3_discrete_integrate, uuGenerativeX1X2X3_1,1/K) - ...
                                                       empcopulaval(C_actual_X1X2X3_discrete_integrate, uuGenerativeX1X2X3_2,1/K) - ...
                                                       empcopulaval(C_actual_X1X2X3_discrete_integrate, uuGenerativeX1X2X3_3,1/K) + ...
                                                       empcopulaval(C_actual_X1X2X3_discrete_integrate, uuGenerativeX1X2X3_4,1/K);
                                    C_est_partial_hat_X1X2X3 = empcopulaval(C_est_X1X2X3_discrete_integrate, uuEstX1X2X3_1,1/K) - ...
                                                           empcopulaval(C_est_X1X2X3_discrete_integrate, uuEstX1X2X3_2,1/K) - ...
                                                           empcopulaval(C_est_X1X2X3_discrete_integrate, uuEstX1X2X3_3,1/K) + ...
                                                           empcopulaval(C_est_X1X2X3_discrete_integrate, uuEstX1X2X3_4,1/K);
                                    C_actual_partial_X1X2 = empcopulaval(C_actual_X1X2_discrete_integrate, uuGenerativeX1X2_1,1/K) - ...
                                                     empcopulaval(C_actual_X1X2_discrete_integrate, uuGenerativeX1X2_2,1/K) - ...
                                                     empcopulaval(C_actual_X1X2_discrete_integrate, uuGenerativeX1X2_3,1/K) + ...
                                                     empcopulaval(C_actual_X1X2_discrete_integrate, uuGenerativeX1X2_4,1/K);
                                    C_est_partial_hat_X1X2 = empcopulaval(C_est_X1X2_discrete_integrate, uuEstX1X2_1,1/K) - ...
                                                         empcopulaval(C_est_X1X2_discrete_integrate, uuEstX1X2_2,1/K) - ...
                                                         empcopulaval(C_est_X1X2_discrete_integrate, uuEstX1X2_3,1/K) + ...
                                                         empcopulaval(C_est_X1X2_discrete_integrate, uuEstX1X2_4,1/K);

    % % %                                 c_hcbn_X1X2X3 = empcopulaval(hcbnObj.copulaFamilies{3}.c, [uuHcbnX1X2X3(3) uuHcbnX1X2X3(1:2)]);   % reverse ordering of uuHcbn to be consistent                                                                                      % w/ how hcbn code estimates copula
    % % %                                 c_hcbn_X1X2 = empcopulaval(hcbnObj.copulaFamilies{3}.c_parents, uuHcbnX1X2);       % I don't think any ordering needs to be reversed here

                                    fX3 = continuousDistInfo.pdf(xiContinuous_val);
                                    fX3_hat = distCEst.pdf(xiContinuous_val);
    % % %                                 fX3_hcbn = hcbnObj.empInfo{3}.pdf(xiContinuous_val);

                                    % assign all copula based valued
                                    fx3_givenx1x2_copula(ii) = C_actual_partial_X1X2X3*fX3/C_actual_partial_X1X2;
                                    fx3_givenx1x2_copulaestf3est(ii) = C_est_partial_hat_X1X2X3*fX3_hat/C_est_partial_hat_X1X2;
                                    fx3_givenx1x2_copulaestf3Actual(ii) = C_est_partial_hat_X1X2X3*fX3/C_est_partial_hat_X1X2;
                                    fx3_givenx1x2_copulaActualf3est(ii) = C_actual_partial_X1X2X3*fX3_hat/C_actual_partial_X1X2;
                                    fx3_givenx1x2_hcbn(ii) = hcbnObj.computeMixedConditionalProbability_( ...
                                        [x1_discrete_conditional x2_discrete_conditional xiContinuous_val], [cc aa bb], cc);

                                    % assign KDE
                                    fx3_givenx1x2_conditionalKDE(ii) = conditionalKDE.pdf(xiContinuous_val);

                                    % assign CLG
                                    % iterate through the CLG nodeBnParams var
                                    % and find the combo of interest
                                    combo = [x1_discrete_conditional x2_discrete_conditional];
                                    nodeBnParams = clgbnObj.bnParams{3};
                                    for clgIdx=1:length(nodeBnParams)
                                        nodeBnParam = nodeBnParams{clgIdx};
                                        if(isequal(nodeBnParam.parentCombination,combo))
                                            % Get the Mean and Covariance
                                            % parameters and break out of loop
                                            Mean = nodeBnParam.Mean;
                                            Covariance = nodeBnParam.Covariance;
                                            break;
                                        end
                                    end
                                    fx3_givenx1x2_clg(ii) = normpdf(xiContinuous_val, Mean, Covariance);

                                    % assign MTE
                                    nodeBnParams = mtebnObj.bnParams{3};
                                    for mteIdx=1:length(nodeBnParams)
                                        nodeBnParam = nodeBnParams{mteIdx};
                                        if(isequal(nodeBnParam.parentCombination,combo))
                                            mteInfo = nodeBnParam.mte_info;
                                            break;
                                        end
                                    end
                                    fx3_givenx1x2_mte(ii) = mteInfo.pdf(xiContinuous_val);

                                end

                                % plot the actual vs the copula version and compare the differences
                                fig1 = figure(1);
                                subplot(1,2,1);
                                plot(xiContinuous, fx3_givenx1x2_copula, 'b*-', ...
                                     xiContinuous, fx3_givenx1x2_copulaestf3est, ...
                                     xiContinuous, fx3_givenx1x2_copulaestf3Actual, ...
                                     xiContinuous, fx3_givenx1x2_copulaActualf3est, ...
                                     xiContinuous, fx3_givenx1x2_hcbn, 'o-'); 
                                grid on; title(sprintf('X_1=%d X_2=%d',x1_discrete_conditional, x2_discrete_conditional));
                                h_legend = legend('$c*f(x_3)$', ...
                                    '$\hat{c}_{RANK}*\hat{f}(x_3)$', ...
                                    '$\hat{c}_{RANK}*f(x_3)$', ...
                                    '$c*\hat{f}(x_3)$', ...
                                    '$\hat{c}_{HCBN-ECDF}$' );
                                set(h_legend,'FontSize',10);
                                set(h_legend,'Interpreter','latex')

                                subplot(1,2,2);
                                plot(xiContinuous, fx3_givenx1x2_copula, 'b*-', ...
                                     xiContinuous, fx3_givenx1x2_conditionalKDE, ...
                                     xiContinuous, fx3_givenx1x2_mte, ...
                                     xiContinuous, fx3_givenx1x2_clg); 
                                grid on; title(sprintf('X_1=%d X_2=%d',x1_discrete_conditional, x2_discrete_conditional));
                                h_legend = legend('$c*f(x_3)$', 'KDE', 'MTE','CLG');
                                set(h_legend,'FontSize',10);
                                set(h_legend,'Interpreter','latex')

                                pause;
                                clf(fig1);
                            end
                        end
                    end
                    
                    % calculate the reference LL
                    refLL = 0;
                    totalProbRefVec = zeros(1,numTest);
                    for ii=1:numTest
                        xx = X_hybrid_test(ii,:);
                        uuGenerativeX1X2X3_1 = [a_dist.cdf(xx(1)) ...
                                                b_dist.cdf(xx(2)) ...
                                                continuousDistInfo.cdf(xx(3))];
                        uuGenerativeX1X2X3_2 = [a_dist.cdf(xx(1)-1) ...
                                                b_dist.cdf(xx(2)) ...
                                                continuousDistInfo.cdf(xx(3))];
                        uuGenerativeX1X2X3_3 = [a_dist.cdf(xx(1)) ...
                                                b_dist.cdf(xx(2)-1) ...
                                                continuousDistInfo.cdf(xx(3))];
                        uuGenerativeX1X2X3_4 = [a_dist.cdf(xx(1)-1) ...
                                                b_dist.cdf(xx(2)-1) ...
                                                continuousDistInfo.cdf(xx(3))];
                        uuGenerativeX1X2_1 = uuGenerativeX1X2X3_1(1:2);
                        uuGenerativeX1X2_2 = uuGenerativeX1X2X3_2(1:2);
                        uuGenerativeX1X2_3 = uuGenerativeX1X2X3_3(1:2);
                        uuGenerativeX1X2_4 = uuGenerativeX1X2X3_4(1:2);
                        
                        fX3 = continuousDistInfo.pdf(xx(3));
                        C_actual_partial_X1X2X3 = empcopulaval(C_actual_X1X2X3_discrete_integrate, uuGenerativeX1X2X3_1,1/K) - ...
                                                  empcopulaval(C_actual_X1X2X3_discrete_integrate, uuGenerativeX1X2X3_2,1/K) - ...
                                                  empcopulaval(C_actual_X1X2X3_discrete_integrate, uuGenerativeX1X2X3_3,1/K) + ...
                                                  empcopulaval(C_actual_X1X2X3_discrete_integrate, uuGenerativeX1X2X3_4,1/K);
                        f_X3X1X2 = C_actual_partial_X1X2X3*fX3;
                        
                        f_X1 = a_dist.pdf(xx(1));
                        f_X2 = b_dist.pdf(xx(2));
                        f_X1X2 = empcopulaval(C_actual_X1X2_discrete_integrate, uuGenerativeX1X2_1, 1/K) - ...
                                 empcopulaval(C_actual_X1X2_discrete_integrate, uuGenerativeX1X2_2, 1/K) - ...
                                 empcopulaval(C_actual_X1X2_discrete_integrate, uuGenerativeX1X2_3, 1/K) + ...
                                 empcopulaval(C_actual_X1X2_discrete_integrate, uuGenerativeX1X2_4, 1/K);
                        f_X3_given_X1X2 = f_X3X1X2/f_X1X2;
                        totalProb = f_X1*f_X2*f_X3_given_X1X2;
                        if(totalProb<1e-5)
                            totalProb = 1e-5;   % PUT BREAK POINT HERE IF YOU WANT TO DEBUG
                        end
                        refLL = refLL + log(totalProb);
                        totalProbRefVec(ii) = totalProb;
                    end
                    
                    % calculate LL values and assign to llDivMCMat
                    hcbnLL = hcbnObj.dataLogLikelihood(X_hybrid_test);
                    hcbnDebugAllLL = hcbnDebugAllObj.dataLogLikelihood(X_hybrid_test);
                    hcbnDebugCopulaLL = hcbnDebugCopulaObj.dataLogLikelihood(X_hybrid_test);
                    hcbnDebugEmpInfoLL = hcbnDebugEmpInfoObj.dataLogLikelihood(X_hybrid_test);
                    mteLL = mtebnObj.dataLogLikelihood(X_hybrid_test);
                    clgLL = clgbnObj.dataLogLikelihood(X_hybrid_test);
                    multinomialLL = multinomialbnObj.dataLogLikelihood(X_hybrid_test);
                    cbnLL = cbnObj.dataLogLikelihood(X_hybrid_test);
                    cbnGaussianLL = cbnGaussianObj.dataLogLikelihood(X_hybrid_test);
                    
                    if(any(isnan([hcbnLL mteLL clgLL multinomialLL cbnLL])))
                        1;      % interactive debugging
                    end
                    
                    % assign LL values to matrix
                    llMCMat(HCBN_LL_MAT_IDX,mcSimNum) = hcbnLL;
                    llMCMat(HCBN_DEBUGALL_LL_MAT_IDX,mcSimNum) = hcbnDebugAllLL;
                    llMCMat(HCBN_DEBUGCOPULA_LL_MAT_IDX,mcSimNum) = hcbnDebugCopulaLL;
                    llMCMat(HCBN_DEBUGEMPINFO_LL_MAT_IDX,mcSimNum) = hcbnDebugEmpInfoLL;
                    llMCMat(MTE_LL_MAT_IDX,mcSimNum) = mteLL;
                    llMCMat(CLG_LL_MAT_IDX,mcSimNum) = clgLL;
                    llMCMat(MULTINOMIAL_LL_MAT_IDX,mcSimNum) = multinomialLL;
                    llMCMat(CBN_LL_MAT_IDX,mcSimNum) = cbnLL;
                    llMCMat(CBN_GAUSSIAN_LL_MAT_IDX,mcSimNum) = cbnGaussianLL;
                    llMCMat(REF_LL_MAT_IDX,mcSimNum) = refLL;
                    
                    progressIdx = progressIdx + 1;
                end
                % average the results from the MC simulation and store
                meanLLDivMCMat = mean(llMCMat,2);
                varLLDivMCMat = mean(llMCMat.^2,2)-meanLLDivMCMat.^2;
                biasLLDivMCMat = mean( repmat(llMCMat(REF_LL_MAT_IDX,:), numModelsCompared, 1) - llMCMat, 2 );
                
                progressStr = sprintf('MEAN{ref}=%f VAR{ref}=%f BIAS{ref}=%f', ...
                    meanLLDivMCMat(REF_LL_MAT_IDX), varLLDivMCMat(REF_LL_MAT_IDX), biasLLDivMCMat(REF_LL_MAT_IDX));
                dispstat(progressStr,'timestamp','keepthis','timestamp');
                progressStr = sprintf('MEAN{hcbn}=%f VAR{hcbn}=%f BIAS{hcbn}=%f', ...
                    meanLLDivMCMat(HCBN_LL_MAT_IDX), varLLDivMCMat(HCBN_LL_MAT_IDX), biasLLDivMCMat(HCBN_LL_MAT_IDX));
                dispstat(progressStr,'timestamp','keepthis','timestamp');
                progressStr = sprintf('MEAN{hcbnDebugCopula}=%f VAR{hcbnDebugCopula}=%f BIAS{hcbnDebugCopula}=%f', ...
                    meanLLDivMCMat(HCBN_DEBUGCOPULA_LL_MAT_IDX), varLLDivMCMat(HCBN_DEBUGCOPULA_LL_MAT_IDX), biasLLDivMCMat(HCBN_DEBUGCOPULA_LL_MAT_IDX));
                dispstat(progressStr,'timestamp','keepthis','timestamp');
                progressStr = sprintf('MEAN{hcbnDebugEmpInfo}=%f VAR{hcbnDebugEmpInfo}=%f BIAS{hcbnDebugEmpInfo}=%f', ...
                    meanLLDivMCMat(HCBN_DEBUGEMPINFO_LL_MAT_IDX), varLLDivMCMat(HCBN_DEBUGEMPINFO_LL_MAT_IDX), biasLLDivMCMat(HCBN_DEBUGEMPINFO_LL_MAT_IDX));
                dispstat(progressStr,'timestamp','keepthis','timestamp');
                progressStr = sprintf('MEAN{cbn}=%f VAR{cbn}=%f BIAS{cbn}=%f', ...
                    meanLLDivMCMat(CBN_LL_MAT_IDX), varLLDivMCMat(CBN_LL_MAT_IDX), biasLLDivMCMat(CBN_LL_MAT_IDX));
                dispstat(progressStr,'timestamp','keepthis','timestamp');
                progressStr = sprintf('MEAN{cbnGaussian}=%f VAR{cbnGaussian}=%f BIAS{cbnGaussian}=%f', ...
                    meanLLDivMCMat(CBN_GAUSSIAN_LL_MAT_IDX), varLLDivMCMat(CBN_GAUSSIAN_LL_MAT_IDX), biasLLDivMCMat(CBN_GAUSSIAN_LL_MAT_IDX));
                dispstat(progressStr,'timestamp','keepthis','timestamp');
                progressStr = sprintf('MEAN{mte}=%f VAR{mte}=%f BIAS{mte}=%f', ...
                    meanLLDivMCMat(MTE_LL_MAT_IDX), varLLDivMCMat(MTE_LL_MAT_IDX), biasLLDivMCMat(MTE_LL_MAT_IDX));
                dispstat(progressStr,'timestamp','keepthis','timestamp');
                progressStr = sprintf('MEAN{clg}=%f VAR{clg}=%f BIAS{clg}=%f', ...
                    meanLLDivMCMat(CLG_LL_MAT_IDX), varLLDivMCMat(CLG_LL_MAT_IDX), biasLLDivMCMat(CLG_LL_MAT_IDX));
                dispstat(progressStr,'timestamp','keepthis','timestamp');
                progressStr = sprintf('MEAN{multinomial}=%f VAR{multinomial}=%f BIAS{multinomial}=%f', ...
                    meanLLDivMCMat(MULTINOMIAL_LL_MAT_IDX), varLLDivMCMat(MULTINOMIAL_LL_MAT_IDX), biasLLDivMCMat(MULTINOMIAL_LL_MAT_IDX));
                dispstat(progressStr,'timestamp','keepthis','timestamp');
                progressStr = sprintf('mean{hcbnDebug}==mean{ref}=%d\n', ...
                    abs(meanLLDivMCMat(HCBN_DEBUGALL_LL_MAT_IDX)-meanLLDivMCMat(REF_LL_MAT_IDX)) < .1 );
                dispstat(progressStr,'timestamp','keepthis','timestamp');
                
                llMat(copulaTypeVecIdx, ...
                        continuousDistTypeVecIdx,...
                        alphaVecIdx,...
                        mVecIdx,:) = meanLLDivMCMat(:);
                llVarMat(copulaTypeVecIdx, ...
                         continuousDistTypeVecIdx,...
                         alphaVecIdx,...
                         mVecIdx,:) = varLLDivMCMat(:);
                llBiasMat(copulaTypeVecIdx, ...
                          continuousDistTypeVecIdx,...
                          alphaVecIdx,...
                          mVecIdx,:) = biasLLDivMCMat(:);  
                llMCCell{copulaTypeVecIdx, ...
                          continuousDistTypeVecIdx,...
                          alphaVecIdx,...
                          mVecIdx,1} = llMCMat;
            end
        end
    end
end

end

function [llMat, llVarMat, llBiasMat, llMCCell] = runD4CFG1()
% runD2CFG1 - the multinomial probabilities are evenly distributed
probsA = [0.5 0.5];
probsB = [0.5 0.5];
probsC = [0.5 0.5];
[llMat, llVarMat, llBiasMat, llMCCell] = runD4(probsA, probsB, probsC);
end

function [llMat, llVarMat, llBiasMat, llMCCell] = runD4CFG2()
% runD2CFG1 - the multinomial probabilities are evenly distributed
probsA = [0.7 0.3];
probsB = [0.5 0.5];
probsC = [0.3 0.7];
[llMat, llVarMat, llBiasMat, llMCCell] = runD4(probsA, probsB, probsC);
end

function [llMat, llVarMat, llBiasMat, llMCCell] = runD4CFG3()
% runD2CFG1 - the multinomial probabilities are evenly distributed
probsA = [0.25 0.25 0.25 0.25];
probsB = [0.25 0.25 0.25 0.25];
probsC = [0.25 0.25 0.25 0.25];
[llMat, llVarMat, llBiasMat, llMCCell] = runD4(probsA, probsB, probsC);
end

function [llMat, llVarMat, llBiasMat, llMCCell] = runD4CFG4()
% runD2CFG2 - the multinomial probabilities are skewed left
probsA = [0.5 0.3 0.1 0.1];
probsB = [0.5 0.3 0.1 0.1];
probsC = fliplr([0.5 0.3 0.1 0.1]);
[llMat, llVarMat, llBiasMat, llMCCell] = runD4(probsA, probsB, probsC);
end

function [llMat, llVarMat, llBiasMat, llMCCell] = runD4CFG5()
% runD2CFG1 - the multinomial probabilities are evenly distributed
probsA = .1*ones(1,10);
probsB = .1*ones(1,10);
probsC = .1*ones(1,10);
[llMat, llVarMat, llBiasMat, llMCCell] = runD4(probsA, probsB, probsC);
end

function [llMat, llVarMat, llBiasMat, llMCCell] = runD4CFG6()
% runD2CFG1 - both multinomial probabilities are skewed left & right
probsA = [.25 .2 .15 .1 .05 .05 .05 .05 .05 .05];
probsB = fliplr([.25 .2 .15 .1 .05 .05 .05 .05 .05 .05]);
probsC = [.25 .2 .15 .1 .05 .05 .05 .05 .05 .05];
[llMat, llVarMat, llBiasMat, llMCCell] = runD4(probsA, probsB, probsC);
end

function [llMat, llVarMat, llBiasMat, llMCCell] = runD4(a_probs, b_probs, c_probs)
global mVec copulaTypeVec_4D alphaVec RhoVecs_4D continuousDistTypeVec 
global numModelsCompared numMC bntPath logFile K h numTest
global HCBN_LL_MAT_IDX MTE_LL_MAT_IDX CLG_LL_MAT_IDX REF_LL_MAT_IDX
global MULTINOMIAL_LL_MAT_IDX NUM_DISCRETE_INTERVALS CBN_LL_MAT_IDX
global HCBN_DEBUGALL_LL_MAT_IDX HCBN_DEBUGCOPULA_LL_MAT_IDX HCBN_DEBUGEMPINFO_LL_MAT_IDX
global CBN_GAUSSIAN_LL_MAT_IDX

a_dist = makedist('Multinomial','Probabilities', a_probs);
b_dist = makedist('Multinomial','Probabilities', b_probs);
c_dist = makedist('Multinomial','Probabilities', c_probs);

% setup the graphical model
aa = 1; bb = 2; cc = 3; dd = 4;
D = 4;
dag = zeros(D,D);
dag(aa,dd) = 1;
dag(bb,dd) = 1;
dag(cc,dd) = 1;
discreteNodes = [aa bb cc];
nodeNames = {'A', 'B', 'C', 'D'};
discreteNodeNames = {'A','B', 'C'};

llMCMat = zeros(numModelsCompared,numMC);
llMat = zeros(length(copulaTypeVec_4D),...
              length(continuousDistTypeVec),...
              length(alphaVec),...
              length(mVec),...
              numModelsCompared);
llVarMat = zeros(length(copulaTypeVec_4D),...
              length(continuousDistTypeVec),...
              length(alphaVec),...
              length(mVec),...
              numModelsCompared);
llBiasMat = zeros(length(copulaTypeVec_4D),...
              length(continuousDistTypeVec),...
              length(alphaVec),...
              length(mVec),...
              numModelsCompared);
llMCCell = cell(length(copulaTypeVec_3D),...
              length(continuousDistTypeVec),...
              length(alphaVec),...
              length(mVec),...
              1);
          
fid = fopen(logFile, 'a');
dispstat('','init'); % One time only initialization
dispstat(sprintf('Begining the simulation...'),'keepthis','timestamp');
numTotalLoops = length(copulaTypeVec_4D)*length(continuousDistTypeVec)*length(alphaVec)*length(mVec)*numMC;
progressIdx = 1;
for copulaTypeVecIdx=1:length(copulaTypeVec_4D)
    for continuousDistTypeVecIdx=1:length(continuousDistTypeVec)
        for alphaVecIdx=1:length(alphaVec)
            for mVecIdx=1:length(mVec)
                copulaType = copulaTypeVec_4D{copulaTypeVecIdx};
                continuousDistType = continuousDistTypeVec{continuousDistTypeVecIdx};
                if(strcmp(copulaType, 'Gaussian'))
                    Rho = RhoVecs_4D{alphaVecIdx};
                else
                    error('Unrecognized Copula type!');
                end
                M = mVec(mVecIdx);
                
                progressAmt = progressIdx/numTotalLoops*100;
                if(strcmp(copulaType,'Gaussian'))
                    progressStr = sprintf('copulaType=%s x4DistType=%s rho=%f M=%d || Progress=%0.04f', ...
                                    copulaType, continuousDistType, Rho(1,4), M, progressAmt);
                else
                    error('Unrecognized Copula Type!');
                end
                dispstat(progressStr,'keepthis','timestamp');
                fprintf(fid, progressStr);
                
                if(strcmp(continuousDistType, 'Multimodal'))
                    xContinuous = [normrnd(-2,0.3,1000,1); normrnd(2,0.8,1000,1)];
                elseif(strcmp(continuousDistType, 'Uniform'))
                    xContinuous = unifrnd(-2,2,2000,1);
                elseif(strcmp(continuousDistType, 'UnimodalSkewed'))
                    xContinuous = betarnd(2,5,2000,1);
                elseif(strcmp(continuousDistType, 'Gaussian'))
                    xContinuous = normrnd(2,0.5,2000,1);
                elseif(strcmp(continuousDistType, 'ThickTailed'))
                    xContinuous = trnd(3, 2000, 1);
                else
                    error('Unknown X4 Dist Type!');
                end
                xContinuous = xContinuous(randperm(2000),:);     % permute for evenness of samples

                isdiscrete = 0;
                [fContinous,xiContinuous] = emppdf(xContinuous,isdiscrete);
                FContinuous = empcdf(xContinuous,isdiscrete);
                continuousDistInfo = rvEmpiricalInfo(xiContinuous,fContinous,FContinuous,isdiscrete);
                
                copulaFamilies = cell(1,4);
                copulaFamilies{1} = [];
                copulaFamilies{2} = [];
                copulaFamilies{3} = [];
                tmp = cell(1,2); 
                tmp{1} = copulaType; 
                if(strcmpi(copulaType, 'Gaussian'))
                    % permute the correlation matrix such that it is in the
                    % format of Child,Parent1,Parent2,Parent3
                    tmp{2} = circshift(circshift(Rho,1,1),1,2);
                else
                    error('Unrecognized Copula Type!!');
                end
                copulaFamilies{4} = tmp;
                empInfo = cell(1,4);
                empInfo{1} = a_dist; empInfo{2} = b_dist; empInfo{3} = c_dist;
                empInfo{4} = continuousDistInfo;
                
                %%%%%%%%%%% MAIN SIMULATION CODE %%%%%%%%%%
                for mcSimNum=1:numMC
                    dispstat(sprintf('MC Sim=%d', mcSimNum), 'timestamp');
                    X_hybrid = zeros(M+numTest,D);
                    
                    % generate the copula random variates
                    if(strcmp(copulaType, 'Gaussian'))
                        u_X1X2X3X4 = copularnd('Gaussian', Rho, M+numTest);
                    else
                        error('Copula Type not recognized!\n');
                    end
                    
                    X_hybrid(:,1) = a_dist.icdf(u_X1X2X3X4(:,1));
                    X_hybrid(:,2) = b_dist.icdf(u_X1X2X3X4(:,2));
                    X_hybrid(:,3) = c_dist.icdf(u_X1X2X3X4(:,3));
                    for ii=1:M+numTest
                        X_hybrid(ii,4) = continuousDistInfo.icdf(u_X1X2X3X4(ii,3));
                    end
                    % generate train and test datasets
                    X_hybrid_test = X_hybrid(M+1:end,:);
                    X_hybrid = X_hybrid(1:M,:);
                    
                    hcbnObj = hcbn(bntPath, X_hybrid, nodeNames, discreteNodeNames, K, h, dag);  
                    hcbnDebugAllObj = hcbn(bntPath, X_hybrid, nodeNames, discreteNodeNames, K, h, dag, 'copulaFamilyInput', copulaFamilies, 'empInfoInput', empInfo);  
                    hcbnDebugCopulaObj = hcbn(bntPath, X_hybrid, nodeNames, discreteNodeNames, K, h, dag, 'copulaFamilyInput', copulaFamilies); 
                    hcbnDebugEmpInfoObj = hcbn(bntPath, X_hybrid, nodeNames, discreteNodeNames, K, h, dag, 'empInfoInput', empInfo); 
                    mtebnObj = mtebn(X_hybrid, discreteNodes, dag);
                    clgbnObj = clgbn(X_hybrid, discreteNodes, dag);
                    multinomialbnObj = multinomialbn(X_hybrid, discreteNodes, dag, NUM_DISCRETE_INTERVALS);
                    cbnObj = cbn(bntPath, X_hybrid, nodeNames, dag);
                    cbnGaussianObj = cbn(bntPath, X_hybrid, nodeNames, dag, 1);
                    
                    %%%%%%%%%%%%%%%%%%%%%%%
                    
                    % calculate LL values and assign to llDivMCMat
                    hcbnLL = hcbnObj.dataLogLikelihood(X_hybrid_test);
                    hcbnDebugAllLL = hcbnDebugAllObj.dataLogLikelihood(X_hybrid_test);
                    hcbnDebugCopulaLL = hcbnDebugCopulaObj.dataLogLikelihood(X_hybrid_test);
                    hcbnDebugEmpInfoLL = hcbnDebugEmpInfoObj.dataLogLikelihood(X_hybrid_test);
                    mteLL = mtebnObj.dataLogLikelihood(X_hybrid_test);
                    clgLL = clgbnObj.dataLogLikelihood(X_hybrid_test);
                    multinomialLL = multinomialbnObj.dataLogLikelihood(X_hybrid_test);
                    cbnLL = cbnObj.dataLogLikelihood(X_hybrid_test);
                    cbnGaussianLL = cbnGaussianObj.dataLogLikelihood(X_hybrid_test);
                    
                    if(any(isnan([hcbnLL mteLL clgLL multinomialLL cbnLL])))
                        1;      % interactive debugging
                    end
                    
                    % assign LL values to matrix
                    llMCMat(HCBN_LL_MAT_IDX,mcSimNum) = hcbnLL;
                    llMCMat(HCBN_DEBUGALL_LL_MAT_IDX,mcSimNum) = hcbnDebugAllLL;
                    llMCMat(HCBN_DEBUGCOPULA_LL_MAT_IDX,mcSimNum) = hcbnDebugCopulaLL;
                    llMCMat(HCBN_DEBUGEMPINFO_LL_MAT_IDX,mcSimNum) = hcbnDebugEmpInfoLL;
                    llMCMat(MTE_LL_MAT_IDX,mcSimNum) = mteLL;
                    llMCMat(CLG_LL_MAT_IDX,mcSimNum) = clgLL;
                    llMCMat(MULTINOMIAL_LL_MAT_IDX,mcSimNum) = multinomialLL;
                    llMCMat(CBN_LL_MAT_IDX,mcSimNum) = cbnLL;
                    llMCMat(CBN_GAUSSIAN_LL_MAT_IDX,mcSimNum) = cbnGaussianLL;
                    llMCMat(REF_LL_MAT_IDX,mcSimNum) = hcbnDebugAllLL;      % WARNING -- we are using HCBN_DEBUG as REFERENCE!
                                                                            % ENSURE SANITY OF RESULTS BY HAND
                    
                    progressIdx = progressIdx + 1;
                end
                % average the results from the MC simulation and store
                meanLLDivMCMat = mean(llMCMat,2);
                varLLDivMCMat = mean(llMCMat.^2,2)-meanLLDivMCMat.^2;
                biasLLDivMCMat = mean( repmat(llMCMat(REF_LL_MAT_IDX,:), numModelsCompared, 1) - llMCMat, 2 );
                
                progressStr = sprintf('MEAN{ref}=%f VAR{ref}=%f BIAS{ref}=%f', ...
                    meanLLDivMCMat(REF_LL_MAT_IDX), varLLDivMCMat(REF_LL_MAT_IDX), biasLLDivMCMat(REF_LL_MAT_IDX));
                dispstat(progressStr,'timestamp','keepthis','timestamp');
                progressStr = sprintf('MEAN{hcbn}=%f VAR{hcbn}=%f BIAS{hcbn}=%f', ...
                    meanLLDivMCMat(HCBN_LL_MAT_IDX), varLLDivMCMat(HCBN_LL_MAT_IDX), biasLLDivMCMat(HCBN_LL_MAT_IDX));
                dispstat(progressStr,'timestamp','keepthis','timestamp');
                progressStr = sprintf('MEAN{hcbnDebugCopula}=%f VAR{hcbnDebugCopula}=%f BIAS{hcbnDebugCopula}=%f', ...
                    meanLLDivMCMat(HCBN_DEBUGCOPULA_LL_MAT_IDX), varLLDivMCMat(HCBN_DEBUGCOPULA_LL_MAT_IDX), biasLLDivMCMat(HCBN_DEBUGCOPULA_LL_MAT_IDX));
                dispstat(progressStr,'timestamp','keepthis','timestamp');
                progressStr = sprintf('MEAN{hcbnDebugEmpInfo}=%f VAR{hcbnDebugEmpInfo}=%f BIAS{hcbnDebugEmpInfo}=%f', ...
                    meanLLDivMCMat(HCBN_DEBUGEMPINFO_LL_MAT_IDX), varLLDivMCMat(HCBN_DEBUGEMPINFO_LL_MAT_IDX), biasLLDivMCMat(HCBN_DEBUGEMPINFO_LL_MAT_IDX));
                dispstat(progressStr,'timestamp','keepthis','timestamp');
                progressStr = sprintf('MEAN{cbn}=%f VAR{cbn}=%f BIAS{cbn}=%f', ...
                    meanLLDivMCMat(CBN_LL_MAT_IDX), varLLDivMCMat(CBN_LL_MAT_IDX), biasLLDivMCMat(CBN_LL_MAT_IDX));
                dispstat(progressStr,'timestamp','keepthis','timestamp');
                progressStr = sprintf('MEAN{cbnGaussian}=%f VAR{cbnGaussian}=%f BIAS{cbnGaussian}=%f', ...
                    meanLLDivMCMat(CBN_GAUSSIAN_LL_MAT_IDX), varLLDivMCMat(CBN_GAUSSIAN_LL_MAT_IDX), biasLLDivMCMat(CBN_GAUSSIAN_LL_MAT_IDX));
                dispstat(progressStr,'timestamp','keepthis','timestamp');
                progressStr = sprintf('MEAN{mte}=%f VAR{mte}=%f BIAS{mte}=%f', ...
                    meanLLDivMCMat(MTE_LL_MAT_IDX), varLLDivMCMat(MTE_LL_MAT_IDX), biasLLDivMCMat(MTE_LL_MAT_IDX));
                dispstat(progressStr,'timestamp','keepthis','timestamp');
                progressStr = sprintf('MEAN{clg}=%f VAR{clg}=%f BIAS{clg}=%f', ...
                    meanLLDivMCMat(CLG_LL_MAT_IDX), varLLDivMCMat(CLG_LL_MAT_IDX), biasLLDivMCMat(CLG_LL_MAT_IDX));
                dispstat(progressStr,'timestamp','keepthis','timestamp');
                progressStr = sprintf('MEAN{multinomial}=%f VAR{multinomial}=%f BIAS{multinomial}=%f', ...
                    meanLLDivMCMat(MULTINOMIAL_LL_MAT_IDX), varLLDivMCMat(MULTINOMIAL_LL_MAT_IDX), biasLLDivMCMat(MULTINOMIAL_LL_MAT_IDX));
                dispstat(progressStr,'timestamp','keepthis','timestamp');
                progressStr = sprintf('mean{hcbnDebug}==mean{ref}=%d\n', ...
                    abs(meanLLDivMCMat(HCBN_DEBUGALL_LL_MAT_IDX)-meanLLDivMCMat(REF_LL_MAT_IDX)) < .1 );
                dispstat(progressStr,'timestamp','keepthis','timestamp');
                
                llMat(copulaTypeVecIdx, ...
                        continuousDistTypeVecIdx,...
                        alphaVecIdx,...
                        mVecIdx,:) = meanLLDivMCMat(:);
                llVarMat(copulaTypeVecIdx, ...
                         continuousDistTypeVecIdx,...
                         alphaVecIdx,...
                         mVecIdx,:) = varLLDivMCMat(:);
                llBiasMat(copulaTypeVecIdx, ...
                          continuousDistTypeVecIdx,...
                          alphaVecIdx,...
                          mVecIdx,:) = biasLLDivMCMat(:);  
                llMCCell{copulaTypeVecIdx, ...
                          continuousDistTypeVecIdx,...
                          alphaVecIdx,...
                          mVecIdx,1} = llMCMat;
            end
        end
    end
end

end

function [llMat, llVarMat, llBiasMat, llMCCell] = runD5CFG1()
% runD2CFG2 - the multinomial probabilities are skewed right
probsA = [0.5 0.5];
probsB = [0.5 0.5];
[llMat, llVarMat, llBiasMat, llMCCell] = runD5(probsA, probsB);
end

function [llMat, llVarMat, llBiasMat, llMCCell] = runD5CFG2()
% runD2CFG2 - the multinomial probabilities are skewed right and left
% oppositely
probsA = [0.3 0.7];
probsB = fliplr([0.3 0.7]);
[llMat, llVarMat, llBiasMat, llMCCell] = runD5(probsA, probsB);
end

function [llMat, llVarMat, llBiasMat, llMCCell] = runD5CFG3()
% runD2CFG1 - the multinomial probabilities are evenly distributed
probsA = [0.25 0.25 0.25 0.25];
probsB = [0.25 0.25 0.25 0.25];
[llMat, llVarMat, llBiasMat, llMCCell] = runD5(probsA, probsB);
end

function [llMat, llVarMat, llBiasMat, llMCCell] = runD5CFG4()
% runD2CFG2 - the multinomial probabilities are skewed left and right
% oppositely
probsA = [0.5 0.3 0.1 0.1];
probsB = fliplr([0.5 0.3 0.1 0.1]);
[llMat, llVarMat, llBiasMat, llMCCell] = runD5(probsA, probsB);
end

function [llMat, llVarMat, llBiasMat, llMCCell] = runD5CFG5()
% runD2CFG1 - the multinomial probabilities are evenly distributed
probsA = .1*ones(1,10);
probsB = .1*ones(1,10);
[llMat, llVarMat, llBiasMat, llMCCell] = runD5(probsA, probsB);
end

function [llMat, llVarMat, llBiasMat, llMCCell] = runD5CFG6()
% runD2CFG1 - both multinomial probabilities are skewed left & right
probsA = [.25 .2 .15 .1 .05 .05 .05 .05 .05 .05];
probsB = fliplr([.25 .2 .15 .1 .05 .05 .05 .05 .05 .05]);
[llMat, llVarMat, llBiasMat, llMCCell] = runD5(probsA, probsB);
end

function [llMat, llVarMat, llBiasMat, llMCCell] = runD5(a_probs, b_probs)
global mVec
global numModelsCompared numMC bntPath logFile K h numTest
global HCBN_LL_MAT_IDX MTE_LL_MAT_IDX CLG_LL_MAT_IDX REF_LL_MAT_IDX
global MULTINOMIAL_LL_MAT_IDX NUM_DISCRETE_INTERVALS CBN_LL_MAT_IDX
global CDE_combinations C1C2C3_combinations dependency_combinations
global HCBN_DEBUGALL_LL_MAT_IDX HCBN_DEBUGCOPULA_LL_MAT_IDX HCBN_DEBUGEMPINFO_LL_MAT_IDX
global CBN_GAUSSIAN_LL_MAT_IDX

a_dist = makedist('Multinomial','Probabilities', a_probs);
b_dist = makedist('Multinomial','Probabilities', b_probs);

% setup the graphical model
aa = 1; bb = 2; cc = 3; dd = 4; ee = 5;
D = 5;
dag = zeros(D,D);
dag(aa,cc) = 1; dag(aa,dd) = 1;
dag(bb,dd) = 1; dag(bb,ee) = 1;
discreteNodes = [aa bb];
nodeNames = {'A', 'B', 'C', 'D', 'E'};
discreteNodeNames = {'A','B'};

llMCMat = zeros(numModelsCompared,numMC);
llMat = zeros(length(CDE_combinations),...
              length(C1C2C3_combinations),...
              length(dependency_combinations),...
              length(mVec),...
              numModelsCompared);
llVarMat = zeros(length(CDE_combinations),...
              length(C1C2C3_combinations),...
              length(dependency_combinations),...
              length(mVec),...
              numModelsCompared);
llBiasMat = zeros(length(CDE_combinations),...
              length(C1C2C3_combinations),...
              length(dependency_combinations),...
              length(mVec),...
              numModelsCompared);
llMCCell = cell(length(CDE_combinations),...
              length(C1C2C3_combinations),...
              length(dependency_combinations),...
              length(mVec),...
              1);
          
fid = fopen(logFile, 'a');
dispstat('','init'); % One time only initialization
dispstat(sprintf('Begining the simulation...'),'keepthis','timestamp');
numTotalLoops = length(CDE_combinations)*length(C1C2C3_combinations)*length(dependency_combinations)*length(mVec)*numMC;
progressIdx = 1;

u = linspace(0,1,K);
[U1_2, U2_2] = ndgrid(u);
[U1_3, U2_3, U3_3] = ndgrid(u);
for cdeCombinationsVecIdx=1:length(CDE_combinations)
    for c1c2c3CombinationsVecIdx=1:length(C1C2C3_combinations)
        for dependencyCombinationsVecIdx=1:length(dependency_combinations)
            for mVecIdx=1:length(mVec)
                
                marginalDistributionCombinations = CDE_combinations{cdeCombinationsVecIdx};
                dependencyCombinations = dependency_combinations{dependencyCombinationsVecIdx};
                c1c2c3Types = C1C2C3_combinations{c1c2c3CombinationsVecIdx};
                
                % setup the empirical marginal distribution definitions
                empiricalDistSamples = 2000;
                continuousEmpiricalDists = cell(1,3);
                for ii=1:3
                    if(strcmpi(marginalDistributionCombinations{ii},'Gaussian'))
                        xx = normrnd(2,0.5,empiricalDistSamples,1);
                    elseif(strcmpi(marginalDistributionCombinations{ii},'Uniform'))
                        xx = unifrnd(-2,2,empiricalDistSamples,1);
                    elseif(strcmpi(marginalDistributionCombinations{ii},'Multimodal'))
                        xx = [normrnd(-2,0.3,empiricalDistSamples/2,1); normrnd(2,0.8,empiricalDistSamples/2,1)];
                    elseif(strcmpi(marginalDistributionCombinations{ii},'ThickTailed'))
                        xx = trnd(3, empiricalDistSamples, 1);
                    end
                    xx = xx(randperm(empiricalDistSamples),:);     % permute for evenness of samples
                    
                    isdiscrete = 0;
                    [fContinous,xiContinuous] = emppdf(xx,isdiscrete);
                    FContinuous = empcdf(xx,isdiscrete);
                    tmp = rvEmpiricalInfo(xiContinuous,fContinous,FContinuous,isdiscrete);
                    continuousEmpiricalDists{ii} = tmp;
                end
                
                % setup dependency amounts
                copulaDepParams = cell(1,3);
                for ii=1:3
                    if(strcmpi(dependencyCombinations{ii},'Strong'))
                        if(ii==2)
                            % deal w/ gaussian case separately b/c we have
                            % to create a correlation matrix
                            copulaDepParams{ii} = [1 0 -.8; 0 1 .59; -.8 .59 1];
                        else
                            if(strcmpi(c1c2c3Types{ii},'Gaussian'))
                                copulaDepParams{ii} = [1 -0.8; -0.8 1];
                            else
                                copulaDepParams{ii} = 10;       % alpha = 10
                            end
                        end
                    elseif(strcmpi(dependencyCombinations{ii}, 'Weak'))
                        if(ii==2)
                            % deal w/ gaussian case separately b/c we have
                            % to create a correlation matrix
                            copulaDepParams{ii} = [1 0 .2; 0 1 -.1; .2 -.1 1];
                        else
                            if(strcmpi(c1c2c3Types{ii},'Gaussian'))
                                copulaDepParams{ii} = [1 -0.1; -0.1 1];
                            else
                                copulaDepParams{ii} = 1;       % alpha = 1
                            end
                        end
                    end
                end
                dep_C1 = copulaDepParams{1};
                dep_C2 = copulaDepParams{2};
                dep_C3 = copulaDepParams{3};
                
                % compute actual copula's, which will be used for reference
                % likelihood calculations
                % compute c1 - Frank copula
                if(strcmpi(c1c2c3Types{1},'Frank'))    
                    c1 = reshape(frankcopulapdf([U1_2(:) U2_2(:)], dep_C1), K, K);
                elseif(strcmpi(c1c2c3Types{1},'Gaussian'))
                    c1 = reshape(copulapdf('Gaussian', [U1_2(:) U2_2(:)], dep_C1), K, K);
                else
                    error('Unsupported Copula Type!');
                end
                % compute c2 - always Gaussian
                c2 = reshape(copulapdf('Gaussian', [U1_3(:) U2_3(:) U3_3(:)], dep_C2), K, K, K);
                c2_parents = reshape(copulapdf('Gaussian', [U1_2(:) U2_2(:)], dep_C2(1:2,1:2)), K, K);
                % compute c3 - Frank copula
                if(strcmpi(c1c2c3Types{3},'Frank'))    
                    c3 = reshape(frankcopulapdf([U1_2(:) U2_2(:)], dep_C3), K, K);
                elseif(strcmpi(c1c2c3Types{3},'Gaussian'))
                    c3 = reshape(copulapdf('Gaussian', [U1_2(:) U2_2(:)], dep_C3), K, K);
                else
                    error('Unsupported Copula Type!');
                end
                % integrate discrete dimensions so we can calculate
                % discrete probabilities accurately w/ c-volume
                C1_partial_discrete_integrate = cumtrapz(u, c1, 1);
                C2_partial_discrete_integrate = cumtrapz(u, cumtrapz(u, c2, 1), 2);
                C2_parents_partial_discrete_integrate = cumtrapz(u, cumtrapz(u, c2_parents, 1), 2);
                C3_partial_discrete_integrate = cumtrapz(u, c3, 1);
                
                copulaFamilies = cell(1,5);
                copulaFamilies{1} = [];
                copulaFamilies{2} = [];
                tmp = cell(1,2);
                tmp{1} = c1c2c3Types{1}; tmp{2} = dep_C1;
                copulaFamilies{3} = tmp;
                tmp = cell(1,2);
                tmp{1} = c1c2c3Types{3}; tmp{2} = dep_C3;
                copulaFamilies{5} = tmp;
                tmp = cell(1,2); 
                tmp{1} = 'Gaussian'; 
                tmp{2} = circshift(circshift(dep_C2,1,1),1,2);
                copulaFamilies{4} = tmp;
                
                empInfo = cell(1,5);
                empInfo{1} = a_dist; empInfo{2} = b_dist;
                empInfo{3} = continuousEmpiricalDists{1};
                empInfo{4} = continuousEmpiricalDists{2};
                empInfo{5} = continuousEmpiricalDists{3};
                
                M = mVec(mVecIdx);
                
                progressAmt = progressIdx/numTotalLoops*100;
                progressStr = sprintf('f(C),f(D),f(E)=%s,%s,%s C1,C2,C3=%s,%s,%s dep=%s,%s,%s M=%d || Progress=%0.04f', ...
                                marginalDistributionCombinations{1}, marginalDistributionCombinations{2}, marginalDistributionCombinations{3}, ...
                                c1c2c3Types{1}, c1c2c3Types{2}, c1c2c3Types{3}, ...
                                dependencyCombinations{1}, dependencyCombinations{2}, dependencyCombinations{3}, ...
                                M, progressAmt);
                dispstat(progressStr,'keepthis','timestamp');
                fprintf(fid, progressStr);
                
                %%%%%%%%%%% MAIN SIMULATION CODE %%%%%%%%%%
                for mcSimNum=1:numMC
                    dispstat(sprintf('MC Sim=%d', mcSimNum), 'timestamp');
                    % Generate the data from the reference BN structure
                    U_C2 = copularnd('Gaussian', dep_C2, M+numTest);        % (:,1)=A
                                                                            % (:,2)=B
                                                                            % (:,3)=D
                    % generate U_C1
                    u1 = U_C2(:,1);
                    p = rand(M+numTest,1);
                    if(strcmpi(c1c2c3Types{1},'Frank'))    
                        U_C1_2 = -log((exp(-dep_C1.*u1).*(1-p)./p + exp(-dep_C1))./(1 + exp(-dep_C1.*u1).*(1-p)./p))./dep_C1;
                    elseif(strcmpi(c1c2c3Types{1},'Gaussian'))
                        x1 = norminv(u1, 0, 1);
                        U = chol(dep_C1,'upper');
                        x2 = [x1 normrnd(0, 1, length(u1), 1)]*U;
                        U_C1_2 = normcdf(x2(:,2));
                    else
                        error('Unsupported Copula Type!');
                    end
                    
                    % generate U_C3
                    u1 = U_C2(:,2);
                    p = rand(M+numTest,1);
                    if(strcmpi(c1c2c3Types{3},'Frank'))
                        U_C3_2 = -log((exp(-dep_C3.*u1).*(1-p)./p + exp(-dep_C3))./(1 + exp(-dep_C3.*u1).*(1-p)./p))./dep_C3;
                    elseif(strcmpi(c1c2c3Types{3},'Gaussian'))
                        x1 = norminv(u1, 0, 1);
                        U = chol(dep_C3,'upper');
                        x2 = [x1 normrnd(0, 1, length(u1), 1)]*U;
                        U_C3_2 = normcdf(x2(:,2));
                    else
                        error('Unsupported Copula Type!');
                    end
                    
                    % Combine all copula samples
                    U = [U_C2(:,1:2) U_C1_2 U_C2(:,3) U_C3_2];  % order the data as A,B,C,D,E
                                                                % U(:,1)=A
                                                                % U(:,2)=B
                                                                % U(:,3)=C
                                                                % U(:,4)=D
                                                                % U(:,5)=E
                                                                
                    % use probability integral transform to create the hybrid data
                    X_hybrid = zeros(size(U));
                    X_hybrid(:,1) = a_dist.icdf(U(:,1));
                    X_hybrid(:,2) = b_dist.icdf(U(:,2));
                    for ii=1:M+numTest
                        X_hybrid(ii,3) = continuousEmpiricalDists{1}.icdf(U(ii,3));
                        X_hybrid(ii,4) = continuousEmpiricalDists{2}.icdf(U(ii,4));
                        X_hybrid(ii,5) = continuousEmpiricalDists{3}.icdf(U(ii,5));
                    end
                    % generate train and test datasets
                    X_hybrid_test = X_hybrid(M+1:end,:);
                    X_hybrid = X_hybrid(1:M,:);
                    
                    % create models of the data
                    hcbnObj = hcbn(bntPath, X_hybrid, nodeNames, discreteNodeNames, K, h, dag);  
                    hcbnDebugAllObj = hcbn(bntPath, X_hybrid, nodeNames, discreteNodeNames, K, h, dag, 'copulaFamilyInput', copulaFamilies, 'empInfoInput', empInfo);  
                    hcbnDebugCopulaObj = hcbn(bntPath, X_hybrid, nodeNames, discreteNodeNames, K, h, dag, 'copulaFamilyInput', copulaFamilies); 
                    hcbnDebugEmpInfoObj = hcbn(bntPath, X_hybrid, nodeNames, discreteNodeNames, K, h, dag, 'empInfoInput', empInfo); 
                    mtebnObj = mtebn(X_hybrid, discreteNodes, dag);
                    clgbnObj = clgbn(X_hybrid, discreteNodes, dag);
                    multinomialbnObj = multinomialbn(X_hybrid, discreteNodes, dag, NUM_DISCRETE_INTERVALS);
                    cbnObj = cbn(bntPath, X_hybrid, nodeNames, dag);
                    cbnGaussianObj = cbn(bntPath, X_hybrid, nodeNames, dag, 1);
                    
                    % compute the reference likelihood
                    refLL = 0;
                    for ii=1:numTest
                        x_A = X_hybrid_test(ii,1); x_B = X_hybrid_test(ii,2);
                        x_C = X_hybrid_test(ii,3); x_D = X_hybrid_test(ii,4); 
                        x_E = X_hybrid_test(ii,5);
                        
                        uu_AB_D_1 = [a_dist.cdf(x_A) b_dist.cdf(x_B) ...
                                     continuousEmpiricalDists{2}.cdf(x_D)];
                        uu_AB_D_2 = [a_dist.cdf(x_A-1) b_dist.cdf(x_B) ...
                                     continuousEmpiricalDists{2}.cdf(x_D)];
                        uu_AB_D_3 = [a_dist.cdf(x_A) b_dist.cdf(x_B-1) ...
                                     continuousEmpiricalDists{2}.cdf(x_D)];
                        uu_AB_D_4 = [a_dist.cdf(x_A-1) b_dist.cdf(x_B-1) ...
                                     continuousEmpiricalDists{2}.cdf(x_D)];
                        uu_AB_1 = uu_AB_D_1(1:2);
                        uu_AB_2 = uu_AB_D_2(1:2);
                        uu_AB_3 = uu_AB_D_3(1:2);
                        uu_AB_4 = uu_AB_D_4(1:2);
                        
                        uu_A_C_1 = [a_dist.cdf(x_A) continuousEmpiricalDists{1}.cdf(x_C)];
                        uu_A_C_2 = [a_dist.cdf(x_A-1) continuousEmpiricalDists{1}.cdf(x_C)];
                        
                        uu_B_E_1 = [b_dist.cdf(x_B) continuousEmpiricalDists{3}.cdf(x_E)];
                        uu_B_E_2 = [b_dist.cdf(x_B-1) continuousEmpiricalDists{3}.cdf(x_E)];
                        
                        f_A = a_dist.pdf(x_A);
                        f_B = b_dist.pdf(x_B);
                        f_AB = empcopulaval(C2_parents_partial_discrete_integrate, uu_AB_1, 1/K) - ...
                               empcopulaval(C2_parents_partial_discrete_integrate, uu_AB_2, 1/K) - ...
                               empcopulaval(C2_parents_partial_discrete_integrate, uu_AB_3, 1/K) + ...
                               empcopulaval(C2_parents_partial_discrete_integrate, uu_AB_4, 1/K);
                        f_A_C = continuousEmpiricalDists{1}.pdf(x_C) * ...
                            (empcopulaval(C1_partial_discrete_integrate, uu_A_C_1, 1/K) - ...
                             empcopulaval(C1_partial_discrete_integrate, uu_A_C_2, 1/K));
                        f_C_given_A = f_A_C / f_A;
                        f_AB_D = (empcopulaval(C2_partial_discrete_integrate, uu_AB_D_1, 1/K) - ...
                                  empcopulaval(C2_partial_discrete_integrate, uu_AB_D_2, 1/K) - ...
                                  empcopulaval(C2_partial_discrete_integrate, uu_AB_D_3, 1/K) + ...
                                  empcopulaval(C2_partial_discrete_integrate, uu_AB_D_4, 1/K)) * ...
                                  continuousEmpiricalDists{2}.pdf(x_D);
                        f_D_given_AB = f_AB_D/f_AB;
                        f_B_E = continuousEmpiricalDists{3}.pdf(x_E) * ...
                            (empcopulaval(C3_partial_discrete_integrate, uu_B_E_1, 1/K) - ...
                             empcopulaval(C3_partial_discrete_integrate, uu_B_E_2, 1/K));
                        f_E_given_B = f_B_E/f_B;
                        
                        totalProb = f_A * f_B * f_C_given_A * f_D_given_AB * f_E_given_B;
                        if(totalProb<1e-5)
                            totalProb = 1e-5;   % PUT BREAK POINT HERE IF YOU WANT TO DEBUG
                        end
                        refLL = refLL + log(totalProb);
                    end
                    
                    % compute likelihood for each of the models
                    hcbnLL = hcbnObj.dataLogLikelihood(X_hybrid_test);
                    hcbnDebugAllLL = hcbnDebugAllObj.dataLogLikelihood(X_hybrid_test);
                    hcbnDebugCopulaLL = hcbnDebugCopulaObj.dataLogLikelihood(X_hybrid_test);
                    hcbnDebugEmpInfoLL = hcbnDebugEmpInfoObj.dataLogLikelihood(X_hybrid_test);
                    mteLL = mtebnObj.dataLogLikelihood(X_hybrid_test);
                    clgLL = clgbnObj.dataLogLikelihood(X_hybrid_test);
                    multinomialLL = multinomialbnObj.dataLogLikelihood(X_hybrid_test);
                    cbnLL = cbnObj.dataLogLikelihood(X_hybrid_test);
                    cbnGaussianLL = cbnGaussianObj.dataLogLikelihood(X_hybrid_test);
                    
                    % assign LL values to matrix
                    llMCMat(HCBN_LL_MAT_IDX,mcSimNum) = hcbnLL;
                    llMCMat(HCBN_DEBUGALL_LL_MAT_IDX,mcSimNum) = hcbnDebugAllLL;
                    llMCMat(HCBN_DEBUGCOPULA_LL_MAT_IDX,mcSimNum) = hcbnDebugCopulaLL;
                    llMCMat(HCBN_DEBUGEMPINFO_LL_MAT_IDX,mcSimNum) = hcbnDebugEmpInfoLL;
                    llMCMat(MTE_LL_MAT_IDX,mcSimNum) = mteLL;
                    llMCMat(CLG_LL_MAT_IDX,mcSimNum) = clgLL;
                    llMCMat(MULTINOMIAL_LL_MAT_IDX,mcSimNum) = multinomialLL;
                    llMCMat(CBN_LL_MAT_IDX,mcSimNum) = cbnLL;
                    llMCMat(CBN_GAUSSIAN_LL_MAT_IDX,mcSimNum) = cbnGaussianLL;
                    llMCMat(REF_LL_MAT_IDX,mcSimNum) = refLL;
                    
                    progressIdx = progressIdx + 1;
                end
                
                % average the results from the MC simulation and store
                meanLLDivMCMat = mean(llMCMat,2);
                varLLDivMCMat = mean(llMCMat.^2,2)-meanLLDivMCMat.^2;
                biasLLDivMCMat = mean( repmat(llMCMat(REF_LL_MAT_IDX,:), numModelsCompared, 1) - llMCMat, 2 );
                
                progressStr = sprintf('MEAN{ref}=%f VAR{ref}=%f BIAS{ref}=%f', ...
                    meanLLDivMCMat(REF_LL_MAT_IDX), varLLDivMCMat(REF_LL_MAT_IDX), biasLLDivMCMat(REF_LL_MAT_IDX));
                dispstat(progressStr,'timestamp','keepthis','timestamp');
                progressStr = sprintf('MEAN{hcbn}=%f VAR{hcbn}=%f BIAS{hcbn}=%f', ...
                    meanLLDivMCMat(HCBN_LL_MAT_IDX), varLLDivMCMat(HCBN_LL_MAT_IDX), biasLLDivMCMat(HCBN_LL_MAT_IDX));
                dispstat(progressStr,'timestamp','keepthis','timestamp');
                progressStr = sprintf('MEAN{hcbnDebugCopula}=%f VAR{hcbnDebugCopula}=%f BIAS{hcbnDebugCopula}=%f', ...
                    meanLLDivMCMat(HCBN_DEBUGCOPULA_LL_MAT_IDX), varLLDivMCMat(HCBN_DEBUGCOPULA_LL_MAT_IDX), biasLLDivMCMat(HCBN_DEBUGCOPULA_LL_MAT_IDX));
                dispstat(progressStr,'timestamp','keepthis','timestamp');
                progressStr = sprintf('MEAN{hcbnDebugEmpInfo}=%f VAR{hcbnDebugEmpInfo}=%f BIAS{hcbnDebugEmpInfo}=%f', ...
                    meanLLDivMCMat(HCBN_DEBUGEMPINFO_LL_MAT_IDX), varLLDivMCMat(HCBN_DEBUGEMPINFO_LL_MAT_IDX), biasLLDivMCMat(HCBN_DEBUGEMPINFO_LL_MAT_IDX));
                dispstat(progressStr,'timestamp','keepthis','timestamp');
                progressStr = sprintf('MEAN{cbn}=%f VAR{cbn}=%f BIAS{cbn}=%f', ...
                    meanLLDivMCMat(CBN_LL_MAT_IDX), varLLDivMCMat(CBN_LL_MAT_IDX), biasLLDivMCMat(CBN_LL_MAT_IDX));
                dispstat(progressStr,'timestamp','keepthis','timestamp');
                progressStr = sprintf('MEAN{cbnGaussian}=%f VAR{cbnGaussian}=%f BIAS{cbnGaussian}=%f', ...
                    meanLLDivMCMat(CBN_GAUSSIAN_LL_MAT_IDX), varLLDivMCMat(CBN_GAUSSIAN_LL_MAT_IDX), biasLLDivMCMat(CBN_GAUSSIAN_LL_MAT_IDX));
                dispstat(progressStr,'timestamp','keepthis','timestamp');
                progressStr = sprintf('MEAN{mte}=%f VAR{mte}=%f BIAS{mte}=%f', ...
                    meanLLDivMCMat(MTE_LL_MAT_IDX), varLLDivMCMat(MTE_LL_MAT_IDX), biasLLDivMCMat(MTE_LL_MAT_IDX));
                dispstat(progressStr,'timestamp','keepthis','timestamp');
                progressStr = sprintf('MEAN{clg}=%f VAR{clg}=%f BIAS{clg}=%f', ...
                    meanLLDivMCMat(CLG_LL_MAT_IDX), varLLDivMCMat(CLG_LL_MAT_IDX), biasLLDivMCMat(CLG_LL_MAT_IDX));
                dispstat(progressStr,'timestamp','keepthis','timestamp');
                progressStr = sprintf('MEAN{multinomial}=%f VAR{multinomial}=%f BIAS{multinomial}=%f', ...
                    meanLLDivMCMat(MULTINOMIAL_LL_MAT_IDX), varLLDivMCMat(MULTINOMIAL_LL_MAT_IDX), biasLLDivMCMat(MULTINOMIAL_LL_MAT_IDX));
                dispstat(progressStr,'timestamp','keepthis','timestamp');
                progressStr = sprintf('mean{hcbnDebug}==mean{ref}=%d\n', ...
                    abs(meanLLDivMCMat(HCBN_DEBUGALL_LL_MAT_IDX)-meanLLDivMCMat(REF_LL_MAT_IDX)) < .1 );
                dispstat(progressStr,'timestamp','keepthis','timestamp');
                
                llMat(cdeCombinationsVecIdx, ...
                     c1c2c3CombinationsVecIdx,...
                     dependencyCombinationsVecIdx,...
                     mVecIdx,:) = meanLLDivMCMat(:);
                llVarMat(cdeCombinationsVecIdx, ...
                     c1c2c3CombinationsVecIdx,...
                     dependencyCombinationsVecIdx,...
                     mVecIdx,:) = varLLDivMCMat(:);
                llBiasMat(cdeCombinationsVecIdx, ...
                     c1c2c3CombinationsVecIdx,...
                     dependencyCombinationsVecIdx,...
                     mVecIdx,:) = biasLLDivMCMat(:);  
                llMCCell{cdeCombinationsVecIdx, ...
                     c1c2c3CombinationsVecIdx,...
                     dependencyCombinationsVecIdx,...
                     mVecIdx,1} = llMCMat;
            end
        end
    end
end

end