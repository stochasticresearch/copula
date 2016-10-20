% RSDM Power Merge Script

clear;
clc;
dbstop if error;

if(ispc)
    load('C:\\Users\\Kiran\\ownCloud\\PhD\\sim_results\\independence\\rsdmPower_M_775_1500.mat');
    M_vec_toMerge = M_vec;
    corrPower_toMerge = corrPower;
    dcorrPower_toMerge = dcorrPower;
    micePower_toMerge = micePower;
    rdcPower_toMerge = rdcPower;
    rsdmPower_toMerge = rsdmPower;
    load('C:\\Users\\Kiran\\ownCloud\\PhD\\sim_results\\independence\\rsdmPower_M_25_750.mat');
elseif(ismac)
    load('/Users/Kiran/ownCloud/PhD/sim_results/independence/rsdmPower_M_775_1500.mat');
    M_vec_toMerge = M_vec;
    corrPower_toMerge = corrPower;
    dcorrPower_toMerge = dcorrPower;
    micePower_toMerge = micePower;
    rdcPower_toMerge = rdcPower;
    rsdmPower_toMerge = rsdmPower;
    load('/Users/Kiran/ownCloud/PhD/sim_results/independence/rsdmPower_M_25_750.mat');
else
    load('/home/kiran/ownCloud/PhD/sim_results/independence/rsdmPower_M_775_1500.mat');
    M_vec_toMerge = M_vec;
    corrPower_toMerge = corrPower;
    dcorrPower_toMerge = dcorrPower;
    micePower_toMerge = micePower;
    rdcPower_toMerge = rdcPower;
    rsdmPower_toMerge = rsdmPower;
    load('/home/kiran/ownCloud/PhD/sim_results/independence/rsdmPower_M_25_750.mat');
end

% Merge the M_vec
M_vec = [M_vec M_vec_toMerge];

% merge the powers
for typ=1:numDepTests
    for l=num_noise_test_min:num_noise_test_max
        for mm=1:30
            corrPower(typ,l,mm+30) = corrPower_toMerge(typ,l,mm);
            dcorrPower(typ,l,mm+30) = dcorrPower_toMerge(typ,l,mm);
            micePower(typ,l,mm+30) = micePower_toMerge(typ,l,mm);
            rdcPower(typ,l,mm+30) = rdcPower_toMerge(typ,l,mm);
            rsdmPower(typ,l,mm+30) = rsdmPower_toMerge(typ,l,mm);
        end
    end
end

% save as rsdmPower.mat
if(ispc)
    save('C:\\Users\\Kiran\\ownCloud\\PhD\\sim_results\\independence\\rsdmPower.mat');
elseif(ismac)
    save('/Users/Kiran/ownCloud/PhD/sim_results/independence/rsdmPower.mat');
else
    save('/home/kiran/ownCloud/PhD/sim_results/independence/rsdmPower.mat');
end