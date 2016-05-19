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
elseif(isunix)
    saveDir = '/home/kiran/ownCloud/PhD/sim_results';
    logFile = '/tmp/out.log';
elseif(ismac)
    saveDir = '/Users/kiran/ownCloud/PhD/sim_results';
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
load(fullfile(saveDir, 'd2_cfg7_25.mat'));

absBias = abs(llBiasMat);

copulaTypeVec_2D = {'Frank', 'Gaussian'};
alphaVec = [1 10];
continuousDistTypeVec = {'Gaussian', 'Uniform', 'Multimodal', 'ThickTailed'};
mVec = 250;

% we have too much data to present effectively :(, we will suffice by
% averaging over different classes and presenting those results :(



x = 1:length(copulaTypeVec_2D)*length(alphaVec)*length(continuousDistTypeVec)*length(mVec);
y = [];
for copulaTypeVecIdx=1:length(copulaTypeVec_2D)
    for continuousDistTypeVecIdx=1:length(continuousDistTypeVec)
        for alphaVecIdx=1:length(alphaVec)
            for mVecIdx=1:length(mVec)
                copulaType = copulaTypeVec_2D{copulaTypeVecIdx};
                continuousDistType = continuousDistTypeVec{continuousDistTypeVecIdx};
                M = mVec(mVecIdx);
                
                y_tmp = [absBias(copulaTypeVecIdx,continuousDistTypeVecIdx,alphaVecIdx,mVecIdx,HCBN_LL_MAT_IDX), ...
                         absBias(copulaTypeVecIdx,continuousDistTypeVecIdx,alphaVecIdx,mVecIdx,MTE_LL_MAT_IDX), ...
                         absBias(copulaTypeVecIdx,continuousDistTypeVecIdx,alphaVecIdx,mVecIdx,CLG_LL_MAT_IDX), ...
                         absBias(copulaTypeVecIdx,continuousDistTypeVecIdx,alphaVecIdx,mVecIdx,MULTINOMIAL_LL_MAT_IDX), ...
                         absBias(copulaTypeVecIdx,continuousDistTypeVecIdx,alphaVecIdx,mVecIdx,CBN_LL_MAT_IDX) ];
                y = [y; y_tmp];
            end
        end
    end
end
