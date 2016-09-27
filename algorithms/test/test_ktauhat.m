%**************************************************************************
%*                                                                        *
%* Copyright (C) 2016  Kiran Karra <kiran.karra@gmail.com>                *
%*                                                                        *
%* This program is free software: you can redistribute it and/or modify   *
%* it under the terms of the GNU General Public License as published by   *
%* the Free Software Foundation, either version 3 of the License, or      *
%* (at your option) any later version.                                    *
%*                                                                        *
%* This program is distributed in the hope that it will be useful,        *
%* but WITHOUT ANY WARRANTY; without even the implied warranty of         *
%* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          *
%* GNU General Public License for more details.                           *
%*                                                                        *
%* You should have received a copy of the GNU General Public License      *
%* along with this program.  If not, see <http://www.gnu.org/licenses/>.  *
%*                                                                        *
%**************************************************************************

%% Understand how different tau's work w/ discrete and hybrid functional dependencies
clear;
clc;
close all;

% rng(12345);

M = 500;

numDiscreteIntervals = 4;

% optimal parameters for MICe
mine_c = 15;
mine_alpha = 0.6;

% Optimal parameters for RDC
rdc_k = 20;
rdc_s = 1/6;

% Strictly monotonic
x = rand(M,1);
% y = x.^3;
y = x;
xx = discretizeRv(x,numDiscreteIntervals)';
yy = discretizeRv(y,numDiscreteIntervals)';
% yyy = xx.^3;        % both x and y are discrete here
yyy = xx;

% some required but not used configuration parameters
alpha_dontCare = 0.05;
wantplot_dontCare = 0;

% calculate the different kind's of tau for each of these
% continuous/discrete/hybrid scenarios
% we have 5 measures -- tau, tau_b, tau_cj, mic_e, rdc
dep_xC_yC = zeros(1,5);
dep_xC_yD = zeros(1,5);
dep_xD_yC = zeros(1,5);
dep_xD_yD = zeros(1,5);

tau = corr(x,y,'type','kendall');
taub = ktaub([x y], alpha_dontCare, wantplot_dontCare);
tau_hat = ktauhat(x, y);
minestats = mine(x',y',mine_alpha,mine_c,'mic_e');
rdcVal = rdc(x,y,rdc_k,rdc_s);
dep_xC_yC(1:5) = [tau taub tau_hat minestats.mic rdcVal];

tau = corr(x,yy,'type','kendall');
taub = ktaub([x yy], alpha_dontCare, wantplot_dontCare);
tau_hat = ktauhat(x, yy);
minestats = mine(x',yy',mine_alpha,mine_c,'mic_e');
rdcVal = rdc(x,yy,rdc_k,rdc_s);
dep_xC_yD(1:5) = [tau taub tau_hat minestats.mic rdcVal];

tau = corr(xx,y,'type','kendall');
taub = ktaub([xx y], alpha_dontCare, wantplot_dontCare);
tau_hat = ktauhat(xx, y);
minestats = mine(xx',y',mine_alpha,mine_c,'mic_e');
rdcVal = rdc(xx,y,rdc_k,rdc_s);
dep_xD_yC(1:5) = [tau taub tau_hat minestats.mic rdcVal];

tau = corr(xx,yyy,'type','kendall');
taub = ktaub([xx yyy], alpha_dontCare, wantplot_dontCare);
tau_hat = ktauhat(xx, yyy);
minestats = mine(xx',yyy',mine_alpha,mine_c,'mic_e');
rdcVal = rdc(xx,yyy,rdc_k,rdc_s);
dep_xD_yD(1:5) = [tau taub tau_hat minestats.mic rdcVal];

figure;
subplot(2,2,1);
scatter(x,y); grid on;
title(sprintf('$\\tau=%0.02f \\  \\tau_b=%0.02f \\ \\hat{\\tau}=%0.02f \\ MIC_e=%0.02f \\ RDC=%0.02f$', ...
    dep_xC_yC(1), dep_xC_yC(2), dep_xC_yC(3), dep_xC_yC(4), dep_xC_yC(5) ), 'Interpreter', 'Latex');

subplot(2,2,2);
scatter(x,yy); grid on;
title(sprintf('$\\tau=%0.02f \\ \\tau_b=%0.02f \\ \\hat{\\tau}=%0.02f \\ MIC_e=%0.02f \\  RDC=%0.02f$', ...
    dep_xC_yD(1), dep_xC_yD(2), dep_xC_yD(3), dep_xC_yD(4), dep_xC_yD(5) ), 'Interpreter', 'Latex');

subplot(2,2,3);
scatter(xx,y); grid on;
title(sprintf('$\\tau=%0.02f \\ \\tau_b=%0.02f \\ \\hat{\\tau}=%0.02f \\ MIC_e=%0.02f \\ RDC=%0.02f$', ...
    dep_xD_yC(1), dep_xD_yC(2), dep_xD_yC(3), dep_xD_yC(4), dep_xD_yC(5) ), 'Interpreter', 'Latex');

subplot(2,2,4);
scatter(xx,yyy); grid on;
title(sprintf('$\\tau=%0.02f \\ \\tau_b=%0.02f \\ \\hat{\\tau}=%0.02f \\ MIC_e=%0.02f \\ RDC=%0.02f$', ...
    dep_xD_yD(1), dep_xD_yD(2), dep_xD_yD(3), dep_xD_yD(4), dep_xD_yD(5) ), 'Interpreter', 'Latex');

% Strictly counter-monotonic
% Strictly monotonic
x = rand(M,1);
% y = -x.^3;
y = -x;
xx = discretizeRv(x,numDiscreteIntervals)';
yy = discretizeRv(y,numDiscreteIntervals)';
% yyy = -xx.^3;        % both x and y are discrete here
yyy = -xx;

% some required but not used configuration parameters
alpha_dontCare = 0.05;
wantplot_dontCare = 0;

% calculate the different kind's of tau for each of these
% continuous/discrete/hybrid scenarios
% we have 5 measures -- tau, tau_b, tau_cj, mic_e, rdc
dep_xC_yC = zeros(1,5);
dep_xC_yD = zeros(1,5);
dep_xD_yC = zeros(1,5);
dep_xD_yD = zeros(1,5);

tau = corr(x,y,'type','kendall');
taub = ktaub([x y], alpha_dontCare, wantplot_dontCare);
tau_hat = ktauhat(x, y);
minestats = mine(x',y',mine_alpha,mine_c,'mic_e');
rdcVal = rdc(x,y,rdc_k,rdc_s);
dep_xC_yC(1:5) = [tau taub tau_hat minestats.mic rdcVal];

tau = corr(x,yy,'type','kendall');
taub = ktaub([x yy], alpha_dontCare, wantplot_dontCare);
tau_hat = ktauhat(x, yy);
minestats = mine(x',yy',mine_alpha,mine_c,'mic_e');
rdcVal = rdc(x,yy,rdc_k,rdc_s);
dep_xC_yD(1:5) = [tau taub tau_hat minestats.mic rdcVal];

tau = corr(xx,y,'type','kendall');
taub = ktaub([xx y], alpha_dontCare, wantplot_dontCare);
tau_hat = ktauhat(xx, y);
minestats = mine(xx',y',mine_alpha,mine_c,'mic_e');
rdcVal = rdc(xx,y,rdc_k,rdc_s);
dep_xD_yC(1:5) = [tau taub tau_hat minestats.mic rdcVal];

tau = corr(xx,yyy,'type','kendall');
taub = ktaub([xx yyy], alpha_dontCare, wantplot_dontCare);
tau_hat = ktauhat(xx, yyy);
minestats = mine(xx',yyy',mine_alpha,mine_c,'mic_e');
rdcVal = rdc(xx,yyy,rdc_k,rdc_s);
dep_xD_yD(1:5) = [tau taub tau_hat minestats.mic rdcVal];

figure;
subplot(2,2,1);
scatter(x,y); grid on;
title(sprintf('$\\tau=%0.02f \\ \\tau_b=%0.02f \\ \\hat{\\tau}=%0.02f \\ MIC_e=%0.02f \\ RDC=%0.02f$', ...
    dep_xC_yC(1), dep_xC_yC(2), dep_xC_yC(3), dep_xC_yC(4), dep_xC_yC(5) ), 'Interpreter', 'Latex');

subplot(2,2,2);
scatter(x,yy); grid on;
title(sprintf('$\\tau=%0.02f \\ \\tau_b=%0.02f \\ \\hat{\\tau}=%0.02f \\ MIC_e=%0.02f \\ RDC=%0.02f$', ...
    dep_xC_yD(1), dep_xC_yD(2), dep_xC_yD(3), dep_xC_yD(4), dep_xC_yD(5) ), 'Interpreter', 'Latex');

subplot(2,2,3);
scatter(xx,y); grid on;
title(sprintf('$\\tau=%0.02f \\ \\tau_b=%0.02f \\ \\hat{\\tau}=%0.02f \\ MIC_e=%0.02f \\ RDC=%0.02f$', ...
    dep_xD_yC(1), dep_xD_yC(2), dep_xD_yC(3), dep_xD_yC(4), dep_xD_yC(5) ), 'Interpreter', 'Latex');

subplot(2,2,4);
scatter(xx,yyy); grid on;
title(sprintf('$\\tau=%0.02f \\ \\tau_b=%0.02f \\ \\hat{\\tau}=%0.02f \\ MIC_e=%0.02f \\ RDC=%0.02f$', ...
    dep_xD_yD(1), dep_xD_yD(2), dep_xD_yD(3), dep_xD_yD(4), dep_xD_yD(5) ), 'Interpreter', 'Latex');


%% Test tau-hat with some copula models for discrete data

%% Test tau-hat with some copula models for hybrid data
clear;
clc;

% rng(12345);

alpha = 5;
rhoVal = 0.5;
rho = [1 rhoVal; rhoVal 1];
M = 500;

% some required but not used configuration parameters
alpha_dontCare = 0.05;
wantplot_dontCare = 0;

% optimal parameters for MICe
mine_c = 15;
mine_alpha = 0.6;

% Optimal parameters for RDC
rdc_k = 20;
rdc_s = 1/6;

% the types of continuous distributions we will be generating for the
% continuous random variable in the bivariate hybrid model
copulaTypes = {'Gaussian', 'Frank', 'Gumbel', 'Clayton'};
continuousDistTypes = {'Gaussian', 'Uniform', 'ThickTailed'};
discreteDistTypes = {[0.25 0.25 0.25 0.25], ...
                       [0.5 0.3 0.1 0.1], ...
                       [0.1 0.1 0.3 0.5]};
% xyOrientationVec = [1 0];
xyOrientationVec = 1;

for copulaTypeIdx=1:length(copulaTypes)
    for continuousDistTypeIdx=1:length(continuousDistTypes)
        for discreteDistTypeIdx=1:length(discreteDistTypes)
            for xyOrientation=xyOrientationVec
                copulaType = copulaTypes{copulaTypeIdx};
                continuousDistType = continuousDistTypes{continuousDistTypeIdx};
                discreteDistType = discreteDistTypes{discreteDistTypeIdx};

                % make the continuous distribution type into a
                % probabilitydistribution object
                if(strcmpi(continuousDistType, 'Gaussian'))
                    pd1 = makedist('Normal');
                elseif(strcmpi(continuousDistType, 'Uniform'))
                    pd1 = makedist('Uniform');
                elseif(strcmpi(continuousDistType, 'ThickTailed'))
                    pd1 = makedist('tLocationScale','mu',0,'sigma',1,'nu',3);
                end

                % make the discrete distribution type into a
                % probabilitydistribution object
                pd2 = makedist('Multinomial','Probabilities',discreteDistType);

                % generate data from the copula
                if(strcmpi(copulaType, 'Gaussian'))
                    U = copularnd(copulaType, rho, M);
                    tauTrue = copulastat(copulaType, rhoVal);
                else
                    U = copularnd(copulaType, alpha, M);
                    tauTrue = copulastat(copulaType, alpha);
                end
                X = zeros(size(U));
                if(xyOrientation==1)
                    for ii=1:M
                        X(ii,1) = pd1.icdf(U(ii,1));
                        X(ii,2) = pd2.icdf(U(ii,2));
                    end
                else
                    for ii=1:M
                        X(ii,2) = pd1.icdf(U(ii,2));
                        X(ii,1) = pd2.icdf(U(ii,1));
                    end
                end
                
                x = X(:,1); y = X(:,2);
                
                tau = corr(x,y,'type','kendall');
%                 taub = ktaub([x y], alpha_dontCare, wantplot_dontCare);
                tau_hat = ktauhat(x, y);
                minestats = mine(x',y',mine_alpha,mine_c,'mic_e');
                rdcVal = rdc(x,y,rdc_k,rdc_s);
                
                subplot(2,1,1);
                scatter(U(:,1),U(:,2)); grid on;
                subplot(2,1,2);
                scatter(x, y); grid on;
                title(sprintf('$[%d,%d,%d] \\tau=%0.02f \\ \\hat{\\tau}=%0.02f \\ MIC_e=%0.02f \\ RDC=%0.02f \\ \\theta=%0.02f$', ...
                        copulaTypeIdx, continuousDistTypeIdx, discreteDistTypeIdx, ...
                        tau, tau_hat, minestats.mic, rdcVal, tauTrue), 'Interpreter', 'Latex');
                pause;
            end
        end
    end
end

%% Characterize the Bias and Variance of ktauhat for hybrid functional and copula dependencies
clear;
clc;

rng(123);

nsim = 300;
M = 500;
numDiscreteIntervals = 4;

correctionFactors = [1 2 3 4 5];
numTotalConfigurations = 7; % 7 entries, 1-5 are the different correction factors, 6 is matlab's built in tau, 7 is tau-b
% TODO: put monotonoic and comonotonic dependencies in here!
dependenciesVec = {'Gaussian', 'Frank', 'Gumbel', 'Clayton'};
subDependenciesVec = [linspace(0.01,0.99,10); ...
                       1:10; ...
                       1:10; ...
                       1:10];
xyOrientationVec = [1 2];       % use [1,2] so this can double as an index into a vector also
pdConfigurations = {'', '', ...
    {'Gaussian', 'Uniform', 'ThickTailed'}, ...
    {[0.25 0.25 0.25 0.25], [0.5 0.3 0.1 0.1], [0.1 0.1 0.3 0.5]} };

resultsXYCopula_bias = zeros( [size(subDependenciesVec) numTotalConfigurations ]);
resultsXYCopula_var = zeros( [size(subDependenciesVec) numTotalConfigurations ]);
resultsYXCopula_bias = zeros( [size(subDependenciesVec) numTotalConfigurations ]);
resultsYXCopula_var = zeros( [size(subDependenciesVec) numTotalConfigurations ]);

for dependencyIdx=1:length(dependenciesVec)
    dependency = dependenciesVec{dependencyIdx};
    subDependencies = subDependenciesVec(dependencyIdx,:);
    for subDependencyIdx=1:length(subDependencies)
        subDependency = subDependencies(subDependencyIdx);
        tauTrue = copulastat(dependency, subDependency);
        % TODO: add more kinds of distributions for pd1 and pd2...
        pd1 = makedist('Normal');
        pd2 = makedist('Multinomial','Probabilities',[0.25 0.25 0.25 0.25]);
        
        for xyOrientation=xyOrientationVec
            
            tau_hat_vec = zeros(numTotalConfigurations, nsim);
            for simnum=1:nsim
            
                % generate the data according to the configuration
                rho = subDependency;
                U = copularnd(dependency, subDependency, M);
                
                X = zeros(size(U));
                if(xyOrientation==1)
                    for ii=1:M
                        X(ii,1) = pd1.icdf(U(ii,1));
                        X(ii,2) = pd2.icdf(U(ii,2));
                    end
                else
                    for ii=1:M
                        X(ii,2) = pd1.icdf(U(ii,2));
                        X(ii,1) = pd2.icdf(U(ii,1));
                    end
                end
                x = X(:,1); y = X(:,2);
                
                % compute the tau
                for correctionFactor=correctionFactors
                    tau_hat = ktauhat(x,y,correctionFactor);
                    tau_hat_vec(correctionFactor, simnum) = tau_hat;
                end
                tau_hat_vec(6,simnum) = corr(x,y,'type','kendall');
                tau_hat_vec(7,simnum) =  ktaub([x y], 0.05, 0);

            end
            
            % compute standardized bias and variance
            tau_hat_bias = mean(tau_hat_vec,2)-tauTrue;
            tau_hat_var  = var(tau_hat_vec-tauTrue,0,2);    % w = 0
            
            % print the result
            fprintf('%s -- %0.02f -- %d\n', ...
                    dependency, subDependency, xyOrientation);
            fprintf('\t Bias >> [%0.02f %0.02f %0.02f %0.02f %0.02f]\n', ...
                tau_hat_bias(1), tau_hat_bias(2), tau_hat_bias(3), tau_hat_bias(4), tau_hat_bias(5));
            fprintf('\t Var  >> [%0.02f %0.02f %0.02f %0.02f %0.02f]\n', ...
                tau_hat_var(1), tau_hat_var(2), tau_hat_var(3), tau_hat_var(4), tau_hat_var(5));
            
            if(xyOrientation)
                resultsXYCopula_bias(dependencyIdx, subDependencyIdx, 1:numTotalConfigurations) = tau_hat_bias;
                resultsXYCopula_var(dependencyIdx, subDependencyIdx, 1:numTotalConfigurations) = tau_hat_var;
            else
                resultsYXCopula_bias(dependencyIdx, subDependencyIdx, 1:numTotalConfigurations) = tau_hat_bias;
                resultsYXCopula_var(dependencyIdx, subDependencyIdx, 1:numTotalConfigurations) = tau_hat_var;
            end
            
        end
    end
end

dependenciesVec_Functional = {'Monotonic', 'Comonotonic'};
subDependenciesVec_Functional = {'Linear', 'Quadratic', 'Exponential'};

resultsXYFunctional_bias = zeros( [size(subDependenciesVec_Functional) numTotalConfigurations ]);
resultsXYFunctional_var = zeros( [size(subDependenciesVec_Functional) numTotalConfigurations ]);
resultsYXFunctional_bias = zeros( [size(subDependenciesVec_Functional) numTotalConfigurations ]);
resultsYXFunctional_var = zeros( [size(subDependenciesVec_Functional) numTotalConfigurations ]);

for dependencyIdx=1:length(dependenciesVec_Functional)
    dependency = dependenciesVec_Functional{dependencyIdx};
    for subDependencyIdx=1:length(subDependenciesVec_Functional)
        subDependency = subDependenciesVec_Functional{subDependencyIdx};
        
        if(strcmpi(dependency, 'Monotonic'))
            tauTrue = 1;
        else
            tauTrue = -1;
        end
        
        for xyOrientation=xyOrientationVec
            tau_hat_vec = zeros(numTotalConfigurations, nsim);
            for simnum=1:nsim
                x = rand(M,1);
                if(strcmpi(subDependency, 'Linear'))
                    y = tauTrue*x;
                elseif(strcmpi(subDependency, 'Quadratic'))
                    y = tauTrue*x.^2;
                elseif(strcmpi(subDependency, 'Exponential'))
                    y = tauTrue*exp(x);
                end

                if(xyOrientation==1)
                    x = discretizeRv(x, numDiscreteIntervals)';
                else
                    y = discretizeRv(y, numDiscreteIntervals)';
                end

                % compute the metrics
                for correctionFactor=correctionFactors
                    tau_hat = ktauhat(x,y,correctionFactor);
                    tau_hat_vec(correctionFactor, simnum) = tau_hat;
                end
                tau_hat_vec(6,simnum) = corr(x,y,'type','kendall');
                tau_hat_vec(7,simnum) =  ktaub([x y], 0.05, 0);
            end
            
            tau_hat_bias = mean(tau_hat_vec,2)-tauTrue;
            tau_hat_var  = var(tau_hat_vec-tauTrue,0,2);    % w = 0
            
            if(xyOrientation==1)
                resultsXYFunctional_bias(dependencyIdx, subDependencyIdx, 1:numTotalConfigurations) = tau_hat_bias;
                resultsXYFunctional_var(dependencyIdx, subDependencyIdx, 1:numTotalConfigurations) = tau_hat_var;
            else
                resultsYXFunctional_bias(dependencyIdx, subDependencyIdx, 1:numTotalConfigurations) = tau_hat_bias;
                resultsYXFunctional_var(dependencyIdx, subDependencyIdx, 1:numTotalConfigurations) = tau_hat_var;
            end
            
            fprintf('%s -- %s -- %d\n', dependency, subDependency, xyOrientation);
        end
    end
end

% save the results
if(ispc)
    save('C:\\Users\\Kiran\\ownCloud\\PhD\\sim_results\\independence\\ktauhat_biasData.mat');
else
    save('/home/kiran/ownCloud/PhD/sim_results/independence/ktauhat_biasData');
end

% plot the data
metricsToPlot = [7 6 4];        % 7 = tau-b, 6 = tau, 4 = tau-h
% manually looking, CF4 seems to work best, plot against tau and tau-b
legendCell = {'\tau_b', '\tau_{CJ}', '\tau_h' };     % put tau-h last, it is the best performing but the bounds
                                                % are harder to see if we put it first
cmap = winter(3);
transparencyFactor = 0.15;

% Gaussian Copula Data
subplot(3,2,1);
depTypeIdx = 1;
x = subDependenciesVec(depTypeIdx,:);
y = squeeze(resultsXYCopula_bias(depTypeIdx,:,metricsToPlot))';
varTmp = squeeze(resultsXYCopula_var(depTypeIdx,:,metricsToPlot))';
e = zeros(size(subDependenciesVec,2), 2, length(metricsToPlot));
for ii=1:length(metricsToPlot)
    e(:,1,ii) = y(ii,:)-varTmp(ii,:)/2; e(:,2,ii) = y(ii,:)+varTmp(ii,:)/2;
end
boundedline(x,y,e,'cmap', cmap, 'transparency', transparencyFactor);
grid on; xlabel('$\rho$', 'Interpreter', 'Latex'); ylabel('$E[\hat{\tau}-\tau]$', 'Interpreter', 'Latex'); 
title('Gaussian Copula');

% Frank Copula Data
subplot(3,2,2);
depTypeIdx = 2;
x = subDependenciesVec(depTypeIdx,:);
y = squeeze(resultsXYCopula_bias(depTypeIdx,:,metricsToPlot))';
varTmp = squeeze(resultsXYCopula_var(depTypeIdx,:,metricsToPlot))';
e = zeros(size(subDependenciesVec,2), 2, 7);
for ii=1:length(metricsToPlot)
    e(:,1,ii) = y(ii,:)-varTmp(ii,:)/2; e(:,2,ii) = y(ii,:)+varTmp(ii,:)/2;
end
boundedline(x,y,e,'cmap', cmap, 'transparency', transparencyFactor);
grid on; xlabel('$\alpha$', 'Interpreter', 'Latex'); ylabel('$E[\hat{\tau}-\tau]$', 'Interpreter', 'Latex'); 
title('Frank Copula');

% Gumbel Copula Data
subplot(3,2,3);
depTypeIdx = 3;
x = subDependenciesVec(depTypeIdx,:);
y = squeeze(resultsXYCopula_bias(depTypeIdx,:,metricsToPlot))';
varTmp = squeeze(resultsXYCopula_var(depTypeIdx,:,metricsToPlot))';
e = zeros(size(subDependenciesVec,2), 2, 7);
for ii=1:length(metricsToPlot)
    e(:,1,ii) = y(ii,:)-varTmp(ii,:)/2; e(:,2,ii) = y(ii,:)+varTmp(ii,:)/2;
end
boundedline(x,y,e,'cmap', cmap, 'transparency', transparencyFactor);
grid on; xlabel('$\alpha$', 'Interpreter', 'Latex'); ylabel('$E[\hat{\tau}-\tau]$', 'Interpreter', 'Latex'); 
title('Gumbel Copula');

% Clayton Copula Data
subplot(3,2,4);
depTypeIdx = 4;
x = subDependenciesVec(depTypeIdx,:);
y = squeeze(resultsXYCopula_bias(depTypeIdx,:,metricsToPlot))';
varTmp = squeeze(resultsXYCopula_var(depTypeIdx,:,metricsToPlot))';
e = zeros(size(subDependenciesVec,2), 2, 7);
for ii=1:length(metricsToPlot)
    e(:,1,ii) = y(ii,:)-varTmp(ii,:)/2; e(:,2,ii) = y(ii,:)+varTmp(ii,:)/2;
end
boundedline(x,y,e,'cmap', cmap, 'transparency', transparencyFactor);
grid on; xlabel('$\alpha$', 'Interpreter', 'Latex'); ylabel('$E[\hat{\tau}-\tau]$', 'Interpreter', 'Latex'); 
title('Clayton Copula');

subplot(3,2,5);
biasValsMonotonicXY = abs(squeeze(resultsXYFunctional_bias(1,:,metricsToPlot)));
varValsMonotonicXY =  squeeze(resultsXYFunctional_var(1,:,metricsToPlot));

% biasValsXY = abs([biasValsMonotonicXY; biasValsComonotonicXY]);
% varValsXY  = [varValsMonotonicXY; varValsComonotonicXY];

% if we have any 0 bias, make it a very miniscule value so that the user
% can see it...
minVal = 0.05e-1;
biasValsMonotonicXY(biasValsMonotonicXY==0) = minVal;
stddevXY = sqrt(varValsMonotonicXY);
width = 1;
% groupnames = {'Linear', 'Quadratic', 'Exponential', '-1*Linear', '-1*Quadratic', '-1*Exponential'};
groupnames = {'Linear', 'Quadratic', 'Exponential'};
bw_title = 'Monotonic Dependency';
bw_xlabel = [];
bw_ylabel = '|E[\tau-\tau_t]|';
bw_colormap = winter;
gridstatus = 'xy';
bw_legend = legendCell;
error_sides = 2;        % change to 1 if you want 1-sided error-bars, but doesn't seem to work properly
legend_type = 'plot';

barweb(biasValsMonotonicXY,stddevXY,width,groupnames,bw_title, bw_xlabel, bw_ylabel, ...
       bw_colormap, gridstatus, bw_legend, error_sides, legend_type);

subplot(3,2,6);
biasValsComonotonicXY = abs(squeeze(resultsXYFunctional_bias(2,:,metricsToPlot)));
varValsComonotonicXY = squeeze(resultsXYFunctional_var(2,:,metricsToPlot));
biasValsComonotonicXY(biasValsComonotonicXY==0) = minVal;
stddevXY = sqrt(varValsComonotonicXY);

bw_title = 'Comonotonic Dependency';
barweb(biasValsComonotonicXY,stddevXY,width,groupnames,bw_title, bw_xlabel, bw_ylabel, ...
       bw_colormap, gridstatus, bw_legend, error_sides, legend_type);