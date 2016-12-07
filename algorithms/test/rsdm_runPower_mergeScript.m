%% Merge the CoS/cCorr/TICe data
clear;
clc;
dbstop if error;

loRangeFilename = 'rsdmPower_CoS_cCorr_ticE_M_25_825.mat';
hiRangeFilename = 'rsdmPower_CoS_cCorr_ticE_M_850_1500.mat';
outputFilename  = 'rsdmPower_CoS_cCorr_ticE_M_25_1500.mat';

if(ispc)
    folder = 'C:\\Users\\Kiran\\ownCloud\\PhD\sim_results\\independence';
elseif(ismac)
    folder = '/Users/Kiran/ownCloud/PhD/sim_results/independence';
else % assume unix
    folder = '/home/kiran/ownCloud/PhD/sim_results/independence';
end

f1 = fullfile(folder, loRangeFilename);
f2 = fullfile(folder, hiRangeFilename);

load(f2);
ccorrPower_toMerge = ccorrPower;
cosPower_toMerge = cosPower;
ticePower_toMerge = ticePower;
load(f1);
mergeStartIdx = length(25:25:825);
numToMerge = length(850:25:1500);

% merge the powers
for typ=1:numDepTests
    for l=num_noise_test_min:num_noise_test_max
        for mm=1:numToMerge
            % for the submerge
            ccorrPower(typ,l,mm+mergeStartIdx) = ccorrPower_toMerge(typ,l,mm+mergeStartIdx);
            cosPower(typ,l,mm+mergeStartIdx) = cosPower_toMerge(typ,l,mm+mergeStartIdx);
            ticePower(typ,l,mm+mergeStartIdx) = ticePower_toMerge(typ,l,mm+mergeStartIdx);
        end
    end
end

% save as rsdmPower.mat
finalOutputFile = fullfile(folder, outputFilename);
clearvars loRangeFilename hiRangeFilename outputFilename
clearvars folder loFilename hiFilename
clearvars ccorrPower_toMerge cosPower_toMerge ticePower_toMerge
clearvars mergeStartIdx numToMerge 
save(finalOutputFile);
%% RSDM Power Merge Script
clear;
clc;
dbstop if error;

% loRangeFilename = 'rsdmPower_M_25_400.mat';
% hiRangeFilename = 'rsdmPower_M_425_750.mat';
% outputFilename  = 'rsdmPower_M_25_750.mat';
loRangeFilename = 'rsdmPower_M_25_750.mat';
hiRangeFilename = 'rsdmPower_M_775_1500.mat';
outputFilename  = 'rsdmPower.mat';

if(ispc)
    folder = 'C:\\Users\\Kiran\\ownCloud\\PhD\sim_results\\independence';
elseif(ismac)
    folder = '/Users/Kiran/ownCloud/PhD/sim_results/independence';
else % assume unix
    folder = '/home/kiran/ownCloud/PhD/sim_results/independence';
end

f1 = fullfile(folder, loRangeFilename);
f2 = fullfile(folder, hiRangeFilename);

load(f2);
corrPower_toMerge = corrPower;
dcorrPower_toMerge = dcorrPower;
micePower_toMerge = micePower;
rdcPower_toMerge = rdcPower;
rsdmPower_toMerge = rsdmPower;
load(f1);
% mergeStartIdx = length(25:25:400);
% numToMerge = length(425:25:750);
% M_vec = 25:25:750;
mergeStartIdx = length(25:25:750);
numToMerge = length(775:25:1500);
M_vec = 25:25:1500;


% Merge the M_vec
% no need to merge this, b/c we controlled the runs by teh idxs directly

% merge the powers
for typ=1:numDepTests
    for l=num_noise_test_min:num_noise_test_max
        for mm=1:numToMerge
            % for the submerge
            corrPower(typ,l,mm+mergeStartIdx) = corrPower_toMerge(typ,l,mm+mergeStartIdx);
            dcorrPower(typ,l,mm+mergeStartIdx) = dcorrPower_toMerge(typ,l,mm+mergeStartIdx);
            micePower(typ,l,mm+mergeStartIdx) = micePower_toMerge(typ,l,mm+mergeStartIdx);
            rdcPower(typ,l,mm+mergeStartIdx) = rdcPower_toMerge(typ,l,mm+mergeStartIdx);
            rsdmPower(typ,l,mm+mergeStartIdx) = rsdmPower_toMerge(typ,l,mm+mergeStartIdx);
        end
    end
end

% save as rsdmPower.mat
finalOutputFile = fullfile(folder, outputFilename);
clearvars loRangeFilename hiRangeFilename outputFilename
clearvars folder loFilename hiFilename
clearvars corrPower_toMerge dcorrPower_toMerge micePower_toMerge rdcPower_toMerge rsdmPower_toMerge
clearvars mergeStartIdx numToMerge 
save(finalOutputFile);

%% Merge the powers of rsdm/dCorr/MICe/rdc/corr with CoS/cCorr/TICe

clear;
clc;
dbstop if error;

f1name = 'power_rsdm_dCorr_rdc_corr_mice.mat';
f2name = 'power_CoS_cCorr_ticE.mat';
outputFilename = 'power_all.mat';

if(ispc)
    folder = 'C:\\Users\\Kiran\\ownCloud\\PhD\sim_results\\independence';
elseif(ismac)
    folder = '/Users/Kiran/ownCloud/PhD/sim_results/independence';
else % assume unix
    folder = '/home/kiran/ownCloud/PhD/sim_results/independence';
end

f1 = fullfile(folder, f1name);
f2 = fullfile(folder, f2name);

load(f2);
load(f1);

% save as rsdmPower.mat
finalOutputFile = fullfile(folder, outputFilename);
clearvars f1name f2name outputFilename
clearvars folder f1 f2
save(finalOutputFile);

%% Merge cosdv power into power_all

clear;
clc;
dbstop if error;

f1name = 'CoSPower_M_25_1500.mat';
f2name = 'power_all.mat';
outputFilename = 'power_all.mat';

if(ispc)
    folder = 'C:\\Users\\Kiran\\ownCloud\\PhD\sim_results\\independence';
elseif(ismac)
    folder = '/Users/Kiran/ownCloud/PhD/sim_results/independence';
else % assume unix
    folder = '/home/kiran/ownCloud/PhD/sim_results/independence';
end

f1 = fullfile(folder, f1name);
f2 = fullfile(folder, f2name);

load(f1);
cosdvPower = cosPower;
load(f2);
cosfpower = cosPower;       % save off old simulations for cosf
cosPower = cosdvPower;      % overwrite old cosf w/ cosdv simulations

finalOutputFile = fullfile(folder, outputFilename);
clearvars f1name f2name outputFilename
clearvars folder f1 f2
save(finalOutputFile);

%% Merge the cosdv power into power_M_500 (i.e. extract the relevant vector)

clear;
clc;
dbstop if error;

f1name = 'power_all.mat';
f2name = 'power_M_500.mat';
outputFilename = 'power_M_500.mat';

if(ispc)
    folder = 'C:\\Users\\Kiran\\ownCloud\\PhD\sim_results\\independence';
elseif(ismac)
    folder = '/Users/Kiran/ownCloud/PhD/sim_results/independence';
else % assume unix
    folder = '/home/kiran/ownCloud/PhD/sim_results/independence';
end

f1 = fullfile(folder, f1name);
f2 = fullfile(folder, f2name);

load(f1);
clearvars -except cosPower f1name f2name outputFilename folder f1 f2
cosdvPower = cosPower;
load(f2);
cosfpower = cosPower;   % save off old simulations for cosf
cosPower = squeeze(cosdvPower(:,:,20));  % overwrite old cosf w/ cosdv simulations, grab M=500 data
clearvars cosdvPower    % no need to store redundant information, makes for confusion later :D

finalOutputFile = fullfile(folder, outputFilename);
clearvars f1name f2name outputFilename
clearvars folder f1 f2
save(finalOutputFile);