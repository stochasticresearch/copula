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

% script which runs all the simulations and produces the plots!

%% Common Test Configuration;

clear;
clc;
numMCSims = 25;
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

%% D=2, CFG=1
D = 2;
cfg = 1;

diary off;
delete(fullfile(saveDir, 'd2_cfg1_25mc_K25.log')); 
diary(fullfile(saveDir, 'd2_cfg1_25mc_K25.log')); 
rng(12345); 
[llMat, llVarMat, llBiasMat, llMCCell] = cmpAllModelsSynthData(D,numMCSims,cfg,logFile);
% save the results off also
save(fullfile(saveDir, 'd2_cfg1_25mc_K25.mat'));

%% D=2, CFG=3
D = 2;
cfg = 3;

diary off;
delete(fullfile(saveDir, 'd2_cfg3_25mc_K25.log')); 
diary(fullfile(saveDir, 'd2_cfg3_25mc_K25.log')); 
rng(12345); 
[llMat, llVarMat, llBiasMat, llMCCell] = cmpAllModelsSynthData(D,numMCSims,cfg,logFile);
% save the results off also
save(fullfile(saveDir, 'd2_cfg3_25mc_K25.mat'));

%% D=2, CFG=7
D = 2;
cfg = 7;

diary off;
delete(fullfile(saveDir, 'd2_cfg7_25mc_K50.log')); 
diary(fullfile(saveDir, 'd2_cfg7_25mc_K50.log')); 
rng(12345); 

% WARNING - use K=50 for cfg=7!

[llMat, llVarMat, llBiasMat, llMCCell] = cmpAllModelsSynthData(D,numMCSims,cfg,logFile, 0, 250);
% save the results off also
save(fullfile(saveDir, 'd2_cfg7_25mc_K50.mat'));

%% D=3, CFG=1
D = 3;
cfg = 1;

diary off;
delete(fullfile(saveDir, 'd3_cfg1_25mc_K25.log')); 
diary(fullfile(saveDir, 'd3_cfg1_25mc_K25.log')); 
rng(12345); 
[llMat, llVarMat, llBiasMat, llMCCell] = cmpAllModelsSynthData(D,numMCSims,cfg,logFile);
% save the results off also
save(fullfile(saveDir, 'd3_cfg1_25mc_K25.mat'));

%% D=4, CFG=1
D = 4;
cfg = 1;

diary off;
delete(fullfile(saveDir, 'd4_cfg1_25mc_K25.log')); 
diary(fullfile(saveDir, 'd4_cfg1_25mc_K25.log')); 
rng(12345); 
[llMat, llVarMat, llBiasMat, llMCCell] = cmpAllModelsSynthData(D,numMCSims,cfg,logFile);
% save the results off also
save(fullfile(saveDir, 'd4_cfg1_25mc_K25.mat'));

%% D=5, CFG=1
D = 5;
cfg = 1;

diary off;
delete(fullfile(saveDir, 'd5_cfg1_25mc_K50.log')); 
diary(fullfile(saveDir, 'd5_cfg1_25mc_K50.log')); 
rng(12345); 
[llMat, llVarMat, llBiasMat, llMCCell] = cmpAllModelsSynthData(D,numMCSims,cfg,logFile);
% save the results off also
save(fullfile(saveDir, 'd5_cfg1_25mc_K50.mat'));

%% D=5, CFG=3
D = 5;
cfg = 3;

diary off;
delete(fullfile(saveDir, 'd5_cfg3_25mc_K50.log')); 
diary(fullfile(saveDir, 'd5_cfg3_25mc_K50.log')); 
rng(12345); 
[llMat, llVarMat, llBiasMat, llMCCell] = cmpAllModelsSynthData(D,numMCSims,cfg,logFile);
% save the results off also
save(fullfile(saveDir, 'd5_cfg3_25mc_K50.mat'));

%% D=5, CFG=5
D = 5;
cfg = 5;

diary off;
delete(fullfile(saveDir, 'd5_cfg5_25mc_K25.log')); 
diary(fullfile(saveDir, 'd5_cfg5_25mc_K25.log')); 
rng(12345); 
[llMat, llVarMat, llBiasMat, llMCCell] = cmpAllModelsSynthData(D,numMCSims,cfg,logFile);
% save the results off also
save(fullfile(saveDir, 'd5_cfg5_25mc_K25.mat'));