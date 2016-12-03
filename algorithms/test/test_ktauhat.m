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
elseif(ismac)
    save('/Users/kiran/ownCloud/PhD/sim_results/independence/ktauhat_biasData');
else
    save('/home/kiran/ownCloud/PhD/sim_results/independence/ktauhat_biasData');
end

f = figure(1);

% plot the data
metricsToPlot = [7 6 4];        % 7 = tau-b, 6 = tau, 4 = tau-h
% manually looking, CF4 seems to work best, plot against tau and tau-b
legendCell = {'\tau_b', '\tau_{CJ}', '\tau_h' };     % put tau-h last, it is the best performing but the bounds
                                                % are harder to see if we put it first
cmap = winter(3);
transparencyFactor = 0.15;

% Gaussian Copula Data
h1 = subplot(3,2,1);
depTypeIdx = 1;
x = subDependenciesVec(depTypeIdx,:);
y = squeeze(resultsXYCopula_bias(depTypeIdx,:,metricsToPlot))';
varTmp = squeeze(resultsXYCopula_var(depTypeIdx,:,metricsToPlot))';
e = zeros(size(subDependenciesVec,2), 2, length(metricsToPlot));
for ii=1:length(metricsToPlot)
    e(:,1,ii) = y(ii,:)-varTmp(ii,:)/2; e(:,2,ii) = y(ii,:)+varTmp(ii,:)/2;
end
boundedline(x,y,e, 'cmap', cmap, 'transparency', transparencyFactor);
grid on; xlabel('{\boldmath$\rho$}', 'Interpreter', 'Latex', 'FontSize', 20); 
ylabel('{\boldmath$E[\hat{\tau}-\tau]$}', 'Interpreter', 'Latex', 'FontSize', 20); 
title('Gaussian Copula', 'FontSize', 24);
a = get(h1,'XTickLabel');
set(h1,'XTickLabel',a,'FontSize',20)
a = get(h1,'YTickLabel');
set(h1,'YTickLabel',a,'FontSize',20)

% Frank Copula Data
h2 = subplot(3,2,2);
depTypeIdx = 2;
x = subDependenciesVec(depTypeIdx,:);
y = squeeze(resultsXYCopula_bias(depTypeIdx,:,metricsToPlot))';
varTmp = squeeze(resultsXYCopula_var(depTypeIdx,:,metricsToPlot))';
e = zeros(size(subDependenciesVec,2), 2, 7);
for ii=1:length(metricsToPlot)
    e(:,1,ii) = y(ii,:)-varTmp(ii,:)/2; e(:,2,ii) = y(ii,:)+varTmp(ii,:)/2;
end
boundedline(x,y,e,'cmap', cmap, 'transparency', transparencyFactor);
grid on; xlabel('{\boldmath$\alpha$}', 'Interpreter', 'Latex', 'FontSize', 20); 
ylabel('{\boldmath$E[\hat{\tau}-\tau]$}', 'Interpreter', 'Latex', 'FontSize', 20); 
title('Frank Copula', 'FontSize', 24); 
a = get(h2,'XTickLabel');
set(h2,'XTickLabel',a,'FontSize',20)
a = get(h2,'YTickLabel');
set(h2,'YTickLabel',a,'FontSize',20)

% Gumbel Copula Data
h3 = subplot(3,2,3);
depTypeIdx = 3;
x = subDependenciesVec(depTypeIdx,:);
y = squeeze(resultsXYCopula_bias(depTypeIdx,:,metricsToPlot))';
varTmp = squeeze(resultsXYCopula_var(depTypeIdx,:,metricsToPlot))';
e = zeros(size(subDependenciesVec,2), 2, 7);
for ii=1:length(metricsToPlot)
    e(:,1,ii) = y(ii,:)-varTmp(ii,:)/2; e(:,2,ii) = y(ii,:)+varTmp(ii,:)/2;
end
boundedline(x,y,e,'cmap', cmap, 'transparency', transparencyFactor);
grid on; xlabel('{\boldmath$\alpha$}', 'Interpreter', 'Latex', 'FontSize', 20); 
ylabel('{\boldmath$E[\hat{\tau}-\tau]$}', 'Interpreter', 'Latex', 'FontSize', 20); 
title('Gumbel Copula', 'FontSize', 24); 
a = get(h3,'XTickLabel');
set(h3,'XTickLabel',a,'FontSize',20)
a = get(h3,'YTickLabel');
set(h3,'YTickLabel',a,'FontSize',20)

% Clayton Copula Data
h4 = subplot(3,2,4);
depTypeIdx = 4;
x = subDependenciesVec(depTypeIdx,:);
y = squeeze(resultsXYCopula_bias(depTypeIdx,:,metricsToPlot))';
varTmp = squeeze(resultsXYCopula_var(depTypeIdx,:,metricsToPlot))';
e = zeros(size(subDependenciesVec,2), 2, 7);
for ii=1:length(metricsToPlot)
    e(:,1,ii) = y(ii,:)-varTmp(ii,:)/2; e(:,2,ii) = y(ii,:)+varTmp(ii,:)/2;
end
boundedline(x,y,e,'cmap', cmap, 'transparency', transparencyFactor);
grid on; xlabel('{\boldmath$\alpha$}', 'Interpreter', 'Latex', 'FontSize', 20); 
ylabel('{\boldmath$E[\hat{\tau}-\tau]$}', 'Interpreter', 'Latex', 'FontSize', 20); 
title('Clayton Copula', 'FontSize', 24); 
a = get(h4,'XTickLabel');
set(h4,'XTickLabel',a,'FontSize',20)
a = get(h4,'YTickLabel');
set(h4,'YTickLabel',a,'FontSize',20)

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
bw_ylabel = '{\boldmath$|E[\hat{\tau}-\tau]|$}';
bw_colormap = winter(3);
gridstatus = 'xy';
bw_legend = legendCell;
error_sides = 2;        % change to 1 if you want 1-sided error-bars, but doesn't seem to work properly
legend_type = 'plot';
legendTextSize = 20;
labelTextSize = 20;
groupTextSize = 20;

barweb(biasValsMonotonicXY,stddevXY,width,groupnames,bw_title, bw_xlabel, bw_ylabel, ...
       bw_colormap, gridstatus, bw_legend, error_sides, legend_type, ...
       legendTextSize, labelTextSize, groupTextSize);

subplot(3,2,6);
biasValsComonotonicXY = abs(squeeze(resultsXYFunctional_bias(2,:,metricsToPlot)));
varValsComonotonicXY = squeeze(resultsXYFunctional_var(2,:,metricsToPlot));
biasValsComonotonicXY(biasValsComonotonicXY==0) = minVal;
stddevXY = sqrt(varValsComonotonicXY);

bw_title = 'Comonotonic Dependency';
barweb(biasValsComonotonicXY,stddevXY,width,groupnames,bw_title, bw_xlabel, bw_ylabel, ...
       bw_colormap, gridstatus, bw_legend, error_sides, legend_type, ...
       legendTextSize, labelTextSize, groupTextSize);
   
%subplotsqueeze(f,1.1);

%% Characterize the null distribution for ktauhat for continuous, discrete, and hybrid

clear;
clc;

rng(12);

nsim = 1000;
M_vec = 100:100:1000;

xMin = 0; xMax = 1;
yMin = 0; yMax = 1;

FIT_PLOTS = 0;

ktauhatNullDistributionResultsContinuous = zeros(nsim, length(M_vec));
ktauHatNullDistributionResultsDiscrete = zeros(nsim, length(M_vec));
ktauHatNullDistributionResultsHybrid1 = zeros(nsim, length(M_vec));
ktauHatNullDistributionResultsHybrid2 = zeros(nsim, length(M_vec));

numDiscreteIntervals = 4;

for ii=1:nsim
    parfor jj=1:length(M_vec)
        M = M_vec(jj);
        
        % continuous independent x & y
        x = rand(M,1)*(xMax-xMin)+xMin;
        y = rand(M,1)*(yMax-yMin)+yMin;
        
        % discrete independent x & y
        x_discrete = discretizeRv(x,numDiscreteIntervals)';
        y_discrete = discretizeRv(y,numDiscreteIntervals)';
        
        ktauhatNullDistributionResultsContinuous(ii,jj) = ktauhat(x,y);
        ktauhatNullDistributionResultsHybrid1(ii,jj) = ktauhat(x_discrete,y);
        ktauhatNullDistributionResultsHybrid2(ii,jj) = ktauhat(x,y_discrete);
        ktauhatNullDistributionResultsDiscrete(ii,jj) = ktauhat(x_discrete,y_discrete);
    end
end

% plot distribution of ktauhat under the null distribution 
legendCell = cell(1,length(M_vec));
subplot(2,2,1);
for ii=1:length(M_vec)
    [f,xi] = ksdensity(ktauhatNullDistributionResultsContinuous(:,ii));
    plot(xi,f); hold on;
    legendCell{ii} = sprintf('M=%d',M_vec(ii));
end
grid on;
legend(legendCell);
title('~ ktauhat, X-C,Y-C');

subplot(2,2,2);
for ii=1:length(M_vec)
    [f,xi] = ksdensity(ktauhatNullDistributionResultsHybrid1(:,ii));
    plot(xi,f); hold on;
    legendCell{ii} = sprintf('M=%d',M_vec(ii));
end
grid on;
legend(legendCell);
title('~ ktauhat, X-D,Y-C');

subplot(2,2,3);
for ii=1:length(M_vec)
    [f,xi] = ksdensity(ktauhatNullDistributionResultsHybrid2(:,ii));
    plot(xi,f); hold on;
    legendCell{ii} = sprintf('M=%d',M_vec(ii));
end
grid on;
legend(legendCell);
title('~ ktauhat, X-C,Y-D');

subplot(2,2,4);
for ii=1:length(M_vec)
    [f,xi] = ksdensity(ktauhatNullDistributionResultsDiscrete(:,ii));
    plot(xi,f); hold on;
    legendCell{ii} = sprintf('M=%d',M_vec(ii));
end
grid on;
legend(legendCell);
title('~ ktauhat, X-D,Y-D');

D_continuous_cell = cell(1,length(M_vec));  PD_continuous_cell = cell(1,length(M_vec));
D_hybrid1_cell = cell(1,length(M_vec));  PD_hybrid1_cell = cell(1,length(M_vec));
D_hybrid2_cell = cell(1,length(M_vec));  PD_hybrid2_cell = cell(1,length(M_vec));
D_discrete_cell = cell(1,length(M_vec));  PD_discrete_cell = cell(1,length(M_vec));
idx = 1;
for ii=1:length(M_vec)
    [D, PD] = allfitdist(ktauhatNullDistributionResultsContinuous(:,ii), 'PDF');
    D_continuous_cell{idx} = D;  PD_continuous_cell{idx} = PD;
    
    [D, PD] = allfitdist(ktauhatNullDistributionResultsHybrid1(:,ii), 'PDF');
    D_hybrid1_cell{idx} = D;  PD_hybrid1_cell{idx} = PD;
    
    [D, PD] = allfitdist(ktauhatNullDistributionResultsHybrid2(:,ii), 'PDF');
    D_hybrid2_cell{idx} = D;  PD_hybrid2_cell{idx} = PD;
    
    [D, PD] = allfitdist(ktauhatNullDistributionResultsDiscrete(:,ii), 'PDF');
    D_discrete_cell{idx} = D;  PD_discrete_cell{idx} = PD;
    
    idx = idx + 1;
end
if(~FIT_PLOTS)
    close all;      % close the generated plots
end

% for each PD type, compute the total BIC score for all sample sizes, and
% choose the best one in that fashion
distributions = {'Beta', 'Birnbaum-Saunders', 'Exponential', ...
                 'Extreme value', 'Gamma', 'Generalized extreme value', ...
                 'Generalized Pareto', 'Inverse Gaussian', 'Logistic', ...
                 'Log-logistic', 'Lognormal', 'Nakagami', 'Normal', ...
                 'Rayleigh', 'Rician', 't location-scale', 'Weibull'};
             
distScoresContinuous = zeros(4,length(distributions));
distScoresHybrid1 = zeros(4,length(distributions));
distScoresHybrid2 = zeros(4,length(distributions));
distScoresDiscrete = zeros(4,length(distributions));
for ii=1:length(distributions)
    dist = distributions{ii};
    % find this distribution in the fit and store the BIC, AIC, AICc scores
    % for all M
    NLogL_continuous = 0; BIC_continuous = 0; AIC_continuous = 0; AICc_continuous = 0;
    NLogL_hybrid1 = 0; BIC_hybrid1 = 0; AIC_hybrid1 = 0; AICc_hybrid1 = 0;
    NLogL_hybrid2 = 0; BIC_hybrid2 = 0; AIC_hybrid2 = 0; AICc_hybrid2 = 0;
    NLogL_discrete = 0; BIC_discrete = 0; AIC_discrete = 0; AICc_discrete = 0;
    for jj=1:length(M_vec)
        D = D_continuous_cell{jj};
        PD = PD_continuous_cell{jj};
        % find the distribution
        for kk=1:length(PD)
            if(strcmpi(PD{kk}.DistributionName, dist))
                break;
            end
        end
        NLogL_continuous = NLogL_continuous + D(kk).NLogL;
        BIC_continuous = BIC_continuous + D(kk).BIC;
        AIC_continuous = AIC_continuous + D(kk).AIC;
        AICc_continuous = AICc_continuous + D(kk).AICc;
        
        D = D_hybrid1_cell{jj};
        PD = PD_hybrid1_cell{jj};
        % find the distribution
        for kk=1:length(PD)
            if(strcmpi(PD{kk}.DistributionName, dist))
                break;
            end
        end
        NLogL_hybrid1 = NLogL_hybrid1 + D(kk).NLogL;
        BIC_hybrid1 = BIC_hybrid1 + D(kk).BIC;
        AIC_hybrid1 = AIC_hybrid1 + D(kk).AIC;
        AICc_hybrid1 = AICc_hybrid1 + D(kk).AICc;
        
        D = D_hybrid2_cell{jj};
        PD = PD_hybrid2_cell{jj};
        % find the distribution
        for kk=1:length(PD)
            if(strcmpi(PD{kk}.DistributionName, dist))
                break;
            end
        end
        NLogL_hybrid2 = NLogL_hybrid2 + D(kk).NLogL;
        BIC_hybrid2 = BIC_hybrid2 + D(kk).BIC;
        AIC_hybrid2 = AIC_hybrid2 + D(kk).AIC;
        AICc_hybrid2 = AICc_hybrid2 + D(kk).AICc;
        
        D = D_discrete_cell{jj};
        PD = PD_discrete_cell{jj};
        % find the distribution
        for kk=1:length(PD)
            if(strcmpi(PD{kk}.DistributionName, dist))
                break;
            end
        end
        NLogL_discrete = NLogL_discrete + D(kk).NLogL;
        BIC_discrete = BIC_discrete + D(kk).BIC;
        AIC_discrete = AIC_discrete + D(kk).AIC;
        AICc_discrete = AICc_discrete + D(kk).AICc;
    end
    
    distScoresContinuous(1,ii) = NLogL_continuous;
    distScoresContinuous(2,ii) = BIC_continuous;
    distScoresContinuous(3,ii) = AIC_continuous;
    distScoresContinuous(4,ii) = AICc_continuous;
    
    distScoresHybrid1(1,ii) = NLogL_hybrid1;
    distScoresHybrid1(2,ii) = BIC_hybrid1;
    distScoresHybrid1(3,ii) = AIC_hybrid1;
    distScoresHybrid1(4,ii) = AICc_hybrid1;
    
    distScoresHybrid2(1,ii) = NLogL_hybrid2;
    distScoresHybrid2(2,ii) = BIC_hybrid2;
    distScoresHybrid2(3,ii) = AIC_hybrid2;
    distScoresHybrid2(4,ii) = AICc_hybrid2;
    
    distScoresDiscrete(1,ii) = NLogL_discrete;
    distScoresDiscrete(2,ii) = BIC_discrete;
    distScoresDiscrete(3,ii) = AIC_discrete;
    distScoresDiscrete(4,ii) = AICc_discrete;
end

% save the data
if(ispc)
    save('C:\\Users\\Kiran\\ownCloud\\PhD\\sim_results\\independence\\ktauhatNullDistribution.mat');
elseif(ismac)
    save('/Users/Kiran/ownCloud/PhD/sim_results/independence/ktauhatNullDistribution.mat');
else
    save('/home/kiran/ownCloud/PhD/sim_results/independence/ktauhatNullDistribution.mat');
end

fprintf('*************** X & Y CONTINUOUS ****************\');
% Sort by NLogL
[~,I] = sort(distScoresContinuous(1,:), 'ascend');
fprintf('NLogL\n');
distributions{I(1)}

% Sort by BIC
[~,I] = sort(distScoresContinuous(2,:), 'ascend');
fprintf('BIC\n');
distributions{I(1)}

% Sort by AIC
[~,I] = sort(distScoresContinuous(3,:), 'ascend');
fprintf('AIC\n');
distributions{I(1)}

% Sort by AICc
[~,I] = sort(distScoresContinuous(4,:), 'ascend');
fprintf('AICc\n');
distributions{I(1)}
fprintf('************************************************\');

fprintf('*************** X & Y HYBRID 1 ****************\');
% Sort by NLogL
[~,I] = sort(distScoresHybrid1(1,:), 'ascend');
fprintf('NLogL\n');
distributions{I(1)}

% Sort by BIC
[~,I] = sort(distScoresHybrid1(2,:), 'ascend');
fprintf('BIC\n');
distributions{I(1)}

% Sort by AIC
[~,I] = sort(distScoresHybrid1(3,:), 'ascend');
fprintf('AIC\n');
distributions{I(1)}

% Sort by AICc
[~,I] = sort(distScoresHybrid1(4,:), 'ascend');
fprintf('AICc\n');
distributions{I(1)}
fprintf('************************************************\');

fprintf('*************** X & Y HYBRID 2 ****************\');
% Sort by NLogL
[~,I] = sort(distScoresHybrid2(1,:), 'ascend');
fprintf('NLogL\n');
distributions{I(1)}

% Sort by BIC
[~,I] = sort(distScoresHybrid2(2,:), 'ascend');
fprintf('BIC\n');
distributions{I(1)}

% Sort by AIC
[~,I] = sort(distScoresHybrid2(3,:), 'ascend');
fprintf('AIC\n');
distributions{I(1)}

% Sort by AICc
[~,I] = sort(distScoresHybrid2(4,:), 'ascend');
fprintf('AICc\n');
distributions{I(1)}
fprintf('************************************************\');

fprintf('*************** X & Y DISCRETE ****************\');
% Sort by NLogL
[~,I] = sort(distScoresDiscrete(1,:), 'ascend');
fprintf('NLogL\n');
distributions{I(1)}

% Sort by BIC
[~,I] = sort(distScoresDiscrete(2,:), 'ascend');
fprintf('BIC\n');
distributions{I(1)}

% Sort by AIC
[~,I] = sort(distScoresDiscrete(3,:), 'ascend');
fprintf('AIC\n');
distributions{I(1)}

% Sort by AICc
[~,I] = sort(distScoresDiscrete(4,:), 'ascend');
fprintf('AICc\n');
distributions{I(1)}
fprintf('************************************************\');