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

% generate plots for paper from the simulations.  Unfortunately, you will
% need to cross-reference cmpAllModelsSynthData.m in order to understand
% why the matrices are indexed the way they are :(

%% common setup

clear;
clc;

if(ispc)
    saveDir = 'C:\Users\Kiran\ownCloud\PhD\sim_results';
    logFile = 'C:\Users\Kiran\Desktop\out.log';
elseif(ismac)
    saveDir = '/Users/kiran/ownCloud/PhD/sim_results';
    logFile = '/tmp/out.log';
elseif(isunix)
    saveDir = '/home/kiran/ownCloud/PhD/sim_results';
    logFile = '/tmp/out.log';
else
end

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

%% Pathological Variance Case D=2
load(fullfile(saveDir, 'd2_cfg7_25mc_K50.mat'));

absBias = abs(llBiasMat);

copulaTypeVec_2D = {'Frank', 'Gaussian'};
alphaVec = [1 10];

% we have too much data to present effectively :(, we will suffice by
% averaging over different classes and presenting those results :(
absBias = squeeze(mean(absBias,2));     % idx2 = continuous dist type, 
                                        % so we are averaging over all 
                                        % different types of distributions 
                                        % for X2
                                        % the squeeze takes out the M
                                        % dimension, which is singular for
                                        % the pathological test case we
                                        % present
stdDevBias = squeeze(mean(sqrt(llVarMat),2));

x = 1:length(copulaTypeVec_2D)*length(alphaVec);
y_bias = [];
y_stddev = [];
for copulaTypeVecIdx=1:length(copulaTypeVec_2D)
    for alphaVecIdx=1:length(alphaVec)
            y_tmp_bias = [absBias(copulaTypeVecIdx,alphaVecIdx,HCBN_LL_MAT_IDX), ...
                          absBias(copulaTypeVecIdx,alphaVecIdx,MTE_LL_MAT_IDX), ...
                          absBias(copulaTypeVecIdx,alphaVecIdx,CLG_LL_MAT_IDX), ...
                          absBias(copulaTypeVecIdx,alphaVecIdx,MULTINOMIAL_LL_MAT_IDX), ...
                          absBias(copulaTypeVecIdx,alphaVecIdx,CBN_LL_MAT_IDX) ];
            y_tmp_stddev =  [stdDevBias(copulaTypeVecIdx,alphaVecIdx,HCBN_LL_MAT_IDX), ...
                          stdDevBias(copulaTypeVecIdx,alphaVecIdx,MTE_LL_MAT_IDX), ...
                          stdDevBias(copulaTypeVecIdx,alphaVecIdx,CLG_LL_MAT_IDX), ...
                          stdDevBias(copulaTypeVecIdx,alphaVecIdx,MULTINOMIAL_LL_MAT_IDX), ...
                          stdDevBias(copulaTypeVecIdx,alphaVecIdx,CBN_LL_MAT_IDX) ];
                      
            y_bias = [y_bias; y_tmp_bias];
            y_stddev  = [y_stddev;  y_tmp_stddev];
    end
end

% generate the plot
width = 1;
groupnames = {'NonLinear Weak Dependency', 'NonLinear Strong Dependency', ...
              'Linear Weak Dependency', 'Linear Strong Dependency'};
bw_title = '2-D Pathological Case';
bw_xlabel = [];
bw_ylabel = 'Bias';
bw_colormap = [];
gridstatus = 'xy';
bw_legend = {'HCBN', 'MTE', 'CLG', 'Multinomial', 'CBN'};
error_sides = 2;        % change to 1 if you want 1-sided error-bars, but doesn't seem to work properly
legend_type = 'plot';
barweb(y_bias,y_stddev,width,groupnames,bw_title, bw_xlabel, bw_ylabel, ...
       bw_colormap, gridstatus, bw_legend, error_sides, legend_type);
   
%% Make plots for D=5, CFG1
load(fullfile(saveDir, 'd5_cfg1_25.mat'));
numCDECombinations = 4;
numC1C2C3Combinations = 2;
numDependencyCombinations = 2;
numMVec = 4;

absBias = abs(llBiasMat);
absBias = squeeze(mean(absBias,1));     % average over all C/D/E combinations
stddevBias = squeeze(mean(sqrt(llVarMat),1));

y_bias = cell(1,numC1C2C3Combinations*numDependencyCombinations);
y_stddev = cell(1,numC1C2C3Combinations*numDependencyCombinations);

y_tmp_idx = 1;
for ii=1:numC1C2C3Combinations
    for jj=1:numDependencyCombinations
        y_tmp_bias = zeros(numMVec,5);       % we will compare 4 models
        y_tmp_stddev  = zeros(numMVec,5);       % we will compare 4 models
        for kk=1:numMVec
            y_tmp_bias(kk,:) = [absBias(ii,jj,kk,HCBN_LL_MAT_IDX), ...
                                       absBias(ii,jj,kk,MTE_LL_MAT_IDX), ...
                                       absBias(ii,jj,kk,CLG_LL_MAT_IDX), ...
                                       absBias(ii,jj,kk,MULTINOMIAL_LL_MAT_IDX), ...
                                       absBias(ii,jj,kk,CBN_LL_MAT_IDX)];
            y_tmp_stddev(kk,:)  = [stddevBias(ii,jj,kk,HCBN_LL_MAT_IDX), ...
                                   stddevBias(ii,jj,kk,MTE_LL_MAT_IDX), ...
                                   stddevBias(ii,jj,kk,CLG_LL_MAT_IDX), ...
                                   stddevBias(ii,jj,kk,MULTINOMIAL_LL_MAT_IDX), ...
                                   stddevBias(ii,jj,kk,CBN_LL_MAT_IDX)];
                                   
        end
        y_bias{y_tmp_idx} = y_tmp_bias;
        y_stddev{y_tmp_idx}  = y_tmp_stddev;
        y_tmp_idx = y_tmp_idx + 1;
    end
end

% define plot options
width = 1;
groupnames = {'M=250', 'M=500', 'M=750', 'M=1000'};
bw_title = [];
bw_xlabel = [];
bw_ylabel = 'Bias';
bw_colormap = parula;
gridstatus = 'xy';
% bw_legend = {'HCBN', 'MTE', 'CLG', 'Multinomial', 'CBN'};
bw_legend = [];     % we make a separate legend and add it to picture for clarity and conciseness
error_sides = 2;        % change to 1 if you want 1-sided error-bars, but doesn't seem to work properly
legend_type = 'plot';
dependencyCombos = {'Strong Linear Dependency', 'Weak Linear Dependency', ...
                    'Strong Non-Linear Dependency', 'Weak Non-Linear Dependency'};
legendTextSize = 12;
labelTextSize = 32;
groupTextSize = 32;
for kk=1:numC1C2C3Combinations*numDependencyCombinations
    subplot(2,2,kk);
    bw_title = sprintf('%s',dependencyCombos{kk});
    barweb(y_bias{kk},y_stddev{kk},width,groupnames,bw_title, bw_xlabel, bw_ylabel, ...
       bw_colormap, gridstatus, bw_legend, error_sides, legend_type, ...
       legendTextSize, labelTextSize, groupTextSize);
end
mtit('CFG=1')

%% Make Plots for D=5, CFG3
load(fullfile(saveDir, 'd5_cfg3_25mc_K25.mat'));
numCDECombinations = 4;
numC1C2C3Combinations = 2;
numDependencyCombinations = 2;
numMVec = 4;

absBias = abs(llBiasMat);
absBias = squeeze(mean(absBias,1));     % average over all C/D/E combinations
stddevBias = squeeze(mean(sqrt(llVarMat),1));

y_bias = cell(1,numC1C2C3Combinations*numDependencyCombinations);
y_stddev = cell(1,numC1C2C3Combinations*numDependencyCombinations);

y_tmp_idx = 1;
for ii=1:numC1C2C3Combinations
    for jj=1:numDependencyCombinations
        y_tmp_bias = zeros(numMVec,5);       % we will compare 4 models
        y_tmp_stddev  = zeros(numMVec,5);       % we will compare 4 models
        for kk=1:numMVec
            y_tmp_bias(kk,:) = [absBias(ii,jj,kk,HCBN_LL_MAT_IDX), ...
                                       absBias(ii,jj,kk,MTE_LL_MAT_IDX), ...
                                       absBias(ii,jj,kk,CLG_LL_MAT_IDX), ...
                                       absBias(ii,jj,kk,MULTINOMIAL_LL_MAT_IDX), ...
                                       absBias(ii,jj,kk,CBN_LL_MAT_IDX)];
            y_tmp_stddev(kk,:)  = [stddevBias(ii,jj,kk,HCBN_LL_MAT_IDX), ...
                                   stddevBias(ii,jj,kk,MTE_LL_MAT_IDX), ...
                                   stddevBias(ii,jj,kk,CLG_LL_MAT_IDX), ...
                                   stddevBias(ii,jj,kk,MULTINOMIAL_LL_MAT_IDX), ...
                                   stddevBias(ii,jj,kk,CBN_LL_MAT_IDX)];
                                   
        end
        y_bias{y_tmp_idx} = y_tmp_bias;
        y_stddev{y_tmp_idx}  = y_tmp_stddev;
        y_tmp_idx = y_tmp_idx + 1;
    end
end

% define plot options
width = 1;
groupnames = {'M=250', 'M=500', 'M=750', 'M=1000'};
bw_title = [];
bw_xlabel = [];
bw_ylabel = 'Bias';
bw_colormap = parula;
gridstatus = 'xy';
% bw_legend = {'HCBN', 'MTE', 'CLG', 'Multinomial', 'CBN'};
bw_legend = [];     % we make a separate legend and add it to picture for clarity and conciseness
error_sides = 2;        % change to 1 if you want 1-sided error-bars, but doesn't seem to work properly
legend_type = 'plot';
dependencyCombos = {'Strong Linear Dependency', 'Weak Linear Dependency', ...
                    'Strong Non-Linear Dependency', 'Weak Non-Linear Dependency'};
legendTextSize = 12;
labelTextSize = 32;
groupTextSize = 32;
for kk=1:numC1C2C3Combinations*numDependencyCombinations
    subplot(2,2,kk);
    bw_title = sprintf('%s',dependencyCombos{kk});
    barweb(y_bias{kk},y_stddev{kk},width,groupnames,bw_title, bw_xlabel, bw_ylabel, ...
       bw_colormap, gridstatus, bw_legend, error_sides, legend_type, ...
       legendTextSize, labelTextSize, groupTextSize);
end
mtit('CFG=3')

%% Make Plots for D=5, CFG5

%% Plots for D=5, CFG 1/3/5 combined

% load cfg1 and store
load(fullfile(saveDir, 'd5_cfg1_25mc_K25.mat'));
cfg1_llBiasMat = llBiasMat;
cfg1_llVarMat = llVarMat;

% we have to reset the save dir between these runs because saveDir is a
% variable that is loaded, possibly from a different OS :0
if(ispc)
    saveDir = 'C:\Users\Kiran\ownCloud\PhD\sim_results';
    logFile = 'C:\Users\Kiran\Desktop\out.log';
elseif(ismac)
    saveDir = '/Users/kiran/ownCloud/PhD/sim_results';
    logFile = '/tmp/out.log';
elseif(isunix)
    saveDir = '/home/kiran/ownCloud/PhD/sim_results';
    logFile = '/tmp/out.log';
else
end

% load cfg3 and store
load(fullfile(saveDir, 'd5_cfg3_25mc_K25.mat'));
cfg3_llBiasMat = llBiasMat;
cfg3_llVarMat = llVarMat;

if(ispc)
    saveDir = 'C:\Users\Kiran\ownCloud\PhD\sim_results';
    logFile = 'C:\Users\Kiran\Desktop\out.log';
elseif(ismac)
    saveDir = '/Users/kiran/ownCloud/PhD/sim_results';
    logFile = '/tmp/out.log';
elseif(isunix)
    saveDir = '/home/kiran/ownCloud/PhD/sim_results';
    logFile = '/tmp/out.log';
else
end

% average cfg1,cfg3,cfg5
llBiasMat = (cfg1_llBiasMat + cfg3_llBiasMat)/2;        % TODO; add CFG5
llVarMat  = (cfg1_llVarMat  + cfg3_llVarMat)/2;         % TODO: add CFG5

% generate the plot for the averaged CFG1/3/5 scenarios
numCDECombinations = 4;
numC1C2C3Combinations = 2;
numDependencyCombinations = 2;
numMVec = 4;

absBias = abs(llBiasMat);
absBias = squeeze(mean(absBias,1));     % average over all C/D/E combinations
stddevBias = squeeze(mean(sqrt(llVarMat),1));

y_bias = cell(1,numC1C2C3Combinations*numDependencyCombinations);
y_stddev = cell(1,numC1C2C3Combinations*numDependencyCombinations);

y_tmp_idx = 1;
for ii=1:numC1C2C3Combinations
    for jj=1:numDependencyCombinations
        y_tmp_bias = zeros(numMVec,4);       % we will compare 4 models
        y_tmp_stddev  = zeros(numMVec,4);       % we will compare 4 models
        for kk=1:numMVec
            y_tmp_bias(kk,:) = [absBias(ii,jj,kk,HCBN_LL_MAT_IDX), ...
                                       absBias(ii,jj,kk,MTE_LL_MAT_IDX), ...
                                       absBias(ii,jj,kk,CLG_LL_MAT_IDX), ...
                                       absBias(ii,jj,kk,MULTINOMIAL_LL_MAT_IDX)];
            y_tmp_stddev(kk,:)  = [stddevBias(ii,jj,kk,HCBN_LL_MAT_IDX), ...
                                   stddevBias(ii,jj,kk,MTE_LL_MAT_IDX), ...
                                   stddevBias(ii,jj,kk,CLG_LL_MAT_IDX), ...
                                   stddevBias(ii,jj,kk,MULTINOMIAL_LL_MAT_IDX)];
                                   
        end
        y_bias{y_tmp_idx} = y_tmp_bias;
        y_stddev{y_tmp_idx}  = y_tmp_stddev;
        y_tmp_idx = y_tmp_idx + 1;
    end
end

% define plot options
width = 1;
groupnames = {'M=250', 'M=500', 'M=750', 'M=1000'};
bw_title = [];
bw_xlabel = [];
bw_ylabel = 'Bias';
bw_colormap = parula;
gridstatus = 'xy';
bw_legend = {'HCBN', 'MTE', 'CLG', 'Multinomial'};
% bw_legend = [];     % we make a separate legend and add it to picture for clarity and conciseness
error_sides = 2;        % change to 1 if you want 1-sided error-bars, but doesn't seem to work properly
legend_type = 'plot';
dependencyCombos = {'Strong Linear Dependency', 'Weak Linear Dependency', ...
                    'Strong Non-Linear Dependency', 'Weak Non-Linear Dependency'};
legendTextSize = 32;
labelTextSize = 32;
groupTextSize = 32;
for kk=1:numC1C2C3Combinations*numDependencyCombinations
    subplot(2,2,kk);
    bw_title = sprintf('%s',dependencyCombos{kk});
    if(kk>1)
        bw_legend = [];
    end
    barweb(y_bias{kk},y_stddev{kk},width,groupnames,bw_title, bw_xlabel, bw_ylabel, ...
       bw_colormap, gridstatus, bw_legend, error_sides, legend_type, ...
       legendTextSize, labelTextSize, groupTextSize);
end