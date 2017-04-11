%**************************************************************************
%* 
%* Copyright (C) 2017  Kiran Karra <kiran.karra@gmail.com>
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
%* along with this program.  If not, see <http://www.gnu.org/licenses/>
%* 
%**************************************************************************

% A script to test the model selection algorithm

%% 2-D tests
clear;
clc;

r = 0.8;
Rho = [1 r; r 1];
alpha = 5;
M = 1000;

U_Gaussian = copularnd('Gaussian', Rho, M);
% U_Frank = frankcopularnd(M, 2, alpha);
U_Gumbel = gumbelcopularnd(M, 2, alpha);
U_Clayton = claytoncopularnd(M, 2, alpha);

% fprintf('----------------------------------------------------\n');
% [modelType, modelParams] = copmodelsel_helm(U_Gaussian)
% fprintf('----------------------------------------------------\n');
% [modelType, modelParams] = copmodelsel_helm(U_Frank)
% fprintf('----------------------------------------------------\n');
% [modelType, modelParams] = copmodelsel_helm(U_Gumbel)
% fprintf('----------------------------------------------------\n');
% [modelType, modelParams] = copmodelsel_helm(U_Clayton)
% fprintf('----------------------------------------------------\n');

fprintf('----------------------------------------------------\n');
modelType = copmodelsel_td(U_Gaussian)
fprintf('----------------------------------------------------\n');
modelType = copmodelsel_td(U_Gumbel)
fprintf('----------------------------------------------------\n');
modelType = copmodelsel_td(U_Clayton)
fprintf('----------------------------------------------------\n');

%% Classification Accuracy test for copmodelsel_td
clear;
clc;

M = 500;
r = 0.1:0.1:.9;
alpha = 1:1:9;
numMCSims = 100;

GaussAccuracy = zeros(1,length(r));
GumbelAccuracy = zeros(1,length(alpha));
ClaytonAccuracy = zeros(1,length(alpha));
for ii=1:length(r)
    rho = r(ii);
    fprintf('Processing Gaussian @ r=%0.02f\n', rho);
    classCorrect = 0;
    for jj=1:numMCSims
        U = copularnd('Gaussian',rho,M);
%         model = copmodelsel_td(U);
        model = copmodelsel_helm(U);  
        if(strcmpi(model,'Gaussian'))
            classCorrect = classCorrect + 1;
        end
    end
    GaussAccuracy(ii) = classCorrect/numMCSims;
end

for ii=1:length(alpha)
    theta = alpha(ii);
    fprintf('Processing Gumbel @ alpha=%0.02f\n', rho);
    classCorrect = 0;
    for jj=1:numMCSims
        U = copularnd('Gumbel',theta,M);
%         model = copmodelsel_td(U);
        model = copmodelsel_helm(U);
        if(strcmpi(model,'Gumbel'))
            classCorrect = classCorrect + 1;
        end
    end
    GumbelAccuracy(ii) = classCorrect/numMCSims;
end

for ii=1:length(alpha)
    theta = alpha(ii);
    fprintf('Processing Gumbel @ Clayton=%0.02f\n', rho);
    classCorrect = 0;
    for jj=1:numMCSims
        U = copularnd('Clayton',theta,M);
%         model = copmodelsel_td(U);
        model = copmodelsel_helm(U);
        if(strcmpi(model,'Clayton'))
            classCorrect = classCorrect + 1;
        end
    end
    ClaytonAccuracy(ii) = classCorrect/numMCSims;
end

plot(r,GaussAccuracy,alpha/10,GumbelAccuracy,alpha/10,ClaytonAccuracy);
xlabel('\rho / \alpha');
ylabel('Accuracy');
legend('Gaussian', 'Gumbel', 'Clayton');
grid on;

%% Generate Confusion Matrix for Copula Model Selection w/ Gridding Approach
clear;
clc;
dbstop if error;

cellfind = @(string)(@(cell_contents)(strcmp(string,cell_contents)));

dispstat('','init'); % One time only initialization
dispstat(sprintf('Begining the simulation...\n'),'keepthis','timestamp');

copulaFamilies = {'Gaussian','Frank','Gumbel','Clayton'};
numMCSims = 1000;
M = 1000;
K = 8;      % binning for identifying model
srhoVec = .25:.05:.95;
results = cell(1,length(srhoVec));
for ii=1:length(srhoVec)
    srho = srhoVec(ii);
    confusionMat = zeros(length(copulaFamilies));
    
    dispstat(sprintf('Computing for srho=%0.02f',srho),'keepthis', 'timestamp');
    
    for jj=1:numMCSims
        dispstat(sprintf('%0.02f',jj/numMCSims*100));
        % choose a copula family randomely
        copFamIdx = randi(length(copulaFamilies));
        cop_family = copulaFamilies{copFamIdx};
        cop_param = copulaparam(cop_family, srho, 'type','spearman');
        U = copularnd(cop_family, cop_param, M);
        selectedModel = copmodelsel_helm(U,K);       
        foundIdx = find(cellfun(cellfind(selectedModel),copulaFamilies));
        confusionMat(copFamIdx,foundIdx) = confusionMat(copFamIdx,foundIdx) + 1;
    end
    % normalize the confusion matrix
    normVec = sum(confusionMat,2);
    confusionMat = confusionMat./repmat(normVec,1,length(confusionMat));
    % store into overall results cell
    results{ii} = confusionMat;
end

if ispc
    saveDir = 'C:\\Users\\Kiran\\ownCloud\\PhD\\sim_results\\modelSelect';
elseif isunix
    % do something else
    saveDir = '/home/kiran/ownCloud/PhD/sim_results/modelSelect';
else
    % do a third thing
    saveDir = '/Users/Kiran/ownCloud/PhD/sim_results/modelSelect';
end
fileName = 'helm_confusion_matrix.mat';

save(fullfile(saveDir,fileName));

% do some analysis
