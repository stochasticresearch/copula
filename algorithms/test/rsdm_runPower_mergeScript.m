% RSDM Power Merge Script

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

loFilename = fullfile(folder, loRangeFilename);
hiFilename = fullfile(folder, hiRangeFilename);

load(hiFilename);
corrPower_toMerge = corrPower;
dcorrPower_toMerge = dcorrPower;
micePower_toMerge = micePower;
rdcPower_toMerge = rdcPower;
rsdmPower_toMerge = rsdmPower;
load(loFilename);
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