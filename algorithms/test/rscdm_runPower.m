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

%% Conditional Dependence Parametric test
clear;
clc;

rng(21);

M = 500;
nsim = 500;

gammaVec = 0:0.1:1;

numDepTypes = 6;
rsdmResultsMat = zeros(length(gammaVec), numDepTypes, nsim);        % for comparision of acceptance rates
rscdmResultsMat = zeros(length(gammaVec), numDepTypes, nsim);
cmaSurrResultsMat = zeros(length(gammaVec), numDepTypes, nsim);     % for comparision of aceptance rate
cmaResultsMat = zeros(length(gammaVec), numDepTypes, nsim);

hdResultsMat = zeros(length(gammaVec), numDepTypes, nsim);
hsicResultsMat = zeros(length(gammaVec), numDepTypes, nsim);
pdcorrResultsMat = zeros(length(gammaVec), numDepTypes, nsim);
pdcorrPValMat = zeros(length(gammaVec), numDepTypes, nsim);
pcorrResultsMat = zeros(length(gammaVec), numDepTypes, nsim);
pcorrPValMat = zeros(length(gammaVec), numDepTypes, nsim);

rsdmResultsVec = zeros(1,nsim);
rscdmResultsVec = zeros(1,nsim);
cmaSurrResultsVec = zeros(1,nsim);
cmaResultsVec = zeros(1,nsim);
hdResultsVec = zeros(1,nsim);
hsicResultsVec = zeros(1,nsim);
pdcorrResultsVec = zeros(1,nsim);
pdcorrPValVec = zeros(1,nsim);
pcorrResultsVec = zeros(1,nsim);
pcorrPValVec = zeros(1,nsim);

dispstat('','init'); % One time only initialization
dispstat(sprintf('Begining the simulation...\n'),'keepthis','timestamp');
for gammaIdx=1:length(gammaVec)
    gamma = gammaVec(gammaIdx);
    dispstat(sprintf('Computing for gamma=%0.02f',gamma),'keepthis', 'timestamp');
    for jj=1:numDepTypes
        for ii=1:nsim
            dispstat(sprintf('Simulating -- %0.02f %%', ii/nsim*100),'timestamp');
            Y = rand(M,1);
            Z = rand(M,1);
            eps = randn(M,1);
            
            switch(jj)
                case 1
                    X = gamma*(Y+Z) + (1-gamma)*eps;
                case 2
                    X = gamma*((Y-0.5).^2 + (Z-0.5).^2) + (1-gamma)*eps;
                case 3
                    X = gamma*(sin(4*pi*Y)+cos(4*pi*Z)) + (1-gamma)*eps;
                case 4
                    X = gamma*(nthroot(Y,4)+nthroot(Z,4)) + (1-gamma)*eps;
                case 5
                    X = gamma*(Y + (Z-0.5).^2) + (1-gamma)*eps;
                case 6
                    X = gamma*((Y-0.5).^2 + cos(4*pi*Z)) + (1-gamma)*eps;
            end
            
            rsdmVal = rsdm(Y,Z);        % for calculating the CI/CD test statistic
            rscdmVal = rscdm(Y,Z,X);
            
            data = struct();
            data.X = Y; data.Y = Z; data.Z = X;
            cmaVal = cassor(data);
            cmaSurrVal = gensurr(data);
            
            hdVal = hd(Y,Z,X);
            hsicVal = hsncic(Y,Z,X);
            [pdcorrVal, pdcorr_pval] = pdcorr_R(Y, Z, X);
            [pcorrVal, pcorr_pval] = partialcorr(Y, Z, X);
            
            rsdmResultsVec(ii) = rsdmVal;
            rscdmResultsVec(ii) = rscdmVal;
            
            cmaSurrResultsVec(ii) = cmaSurrVal;
            cmaResultsVec(ii) = cmaVal;
            
            hdResultsVec(ii) = hdVal;
            hsicResultsVec(ii) = hsicVal;
            pdcorrResultsVec(ii) = pdcorrVal;
            pdcorrPValVec(ii) = pdcorr_pval;
            
            pcorrResultsVec(ii) = pcorrVal;
            pcorrPValVec(ii) = pcorr_pval;
        end
        
        rsdmResultsMat(gammaIdx, jj, :) = rsdmResultsVec;
        rscdmResultsMat(gammaIdx, jj, :) = rscdmResultsVec;
        
        cmaSurrResultsMat(gammaIdx, jj, :) = cmaSurrResultsVec;
        cmaResultsMat(gammaIdx, jj, :) = cmaResultsVec;
        
        hdResultsMat(gammaIdx, jj, :) = hdResultsVec;
        hsicResultsMat(gammaIdx, jj, :) = hsicResultsVec;
        
        pdcorrResultsMat(gammaIdx, jj, :) = pdcorrResultsVec;
        pdcorrPValMat(gammaIdx, jj, :) = pdcorrPValVec;
        
        pcorrResultsMat(gammaIdx, jj, :) = pcorrResultsVec;
        pcorrPValMat(gammaIdx, jj, :) = pcorrPValVec;
    end
end

% TODO: calculate acceptace rates for RSDM/RSCDM based on calculating the
% test statistic ...

% To understand the acceptance for CMA, look at demo.m in the cmi folder to
% see how to use the "surrogate" versus the actual test statistic to see if
% there is an accept or reject.

% save the results
if(ispc)
    save('C:\\Users\\Kiran\\ownCloud\\PhD\\sim_results\\independence\\rscdm_CD.mat');
elseif(ismac)
    save('/Users/Kiran/ownCloud/PhD/sim_results/independence/rscdm_CD.mat');
elseif(isunix)
    save('/home/kiran/ownCloud/PhD/sim_results/independence/rscdm_CD.mat');
end

% plot the dependence metric results versus gamma for each dep type
figure;

h1 = subplot(3,2,1);
depIdx = 1;
meanRscdmRes = mean(squeeze(rscdmResultsMat(:,depIdx,:)),2);
meanPdcorrRes = mean(squeeze(pdcorrResultsMat(:,depIdx,:)),2);
meanPcorrRes = mean(squeeze(pcorrResultsMat(:,depIdx,:)),2);
meanCmaRes = mean(squeeze(cmaResultsMat(:,depIdx,:)),2);
plot(gammaVec, meanRscdmRes, 'o-.', ...
     gammaVec, meanPdcorrRes, 'x-.', ...
     gammaVec, meanPcorrRes, 'd-.', ...
     gammaVec, meanCmaRes, '+-.');
xlabel('\gamma', 'FontSize', 20); ylabel('\mu[DEP({X,Y}|Z)]', 'FontSize', 20); grid on;
h1.FontSize = 20; 

h2 = subplot(3,2,2);
depIdx = 2;
meanRscdmRes = mean(squeeze(rscdmResultsMat(:,depIdx,:)),2);
meanPdcorrRes = mean(squeeze(pdcorrResultsMat(:,depIdx,:)),2);
meanPcorrRes = mean(squeeze(pcorrResultsMat(:,depIdx,:)),2);
meanCmaRes = mean(squeeze(cmaResultsMat(:,depIdx,:)),2);
plot(gammaVec, meanRscdmRes, 'o-.', ...
     gammaVec, meanPdcorrRes, 'x-.', ...
     gammaVec, meanPcorrRes, 'd-.', ...
     gammaVec, meanCmaRes, '+-.');
xlabel('\gamma', 'FontSize', 20); ylabel('\mu[DEP({X,Y}|Z)]', 'FontSize', 20); grid on;
h2.FontSize = 20;

h3 = subplot(3,2,3);
depIdx = 3;
meanRscdmRes = mean(squeeze(rscdmResultsMat(:,depIdx,:)),2);
meanPdcorrRes = mean(squeeze(pdcorrResultsMat(:,depIdx,:)),2);
meanPcorrRes = mean(squeeze(pcorrResultsMat(:,depIdx,:)),2);
meanCmaRes = mean(squeeze(cmaResultsMat(:,depIdx,:)),2);
plot(gammaVec, meanRscdmRes, 'o-.', ...
     gammaVec, meanPdcorrRes, 'x-.', ...
     gammaVec, meanPcorrRes, 'd-.', ...
     gammaVec, meanCmaRes, '+-.');
xlabel('\gamma', 'FontSize', 20); ylabel('\mu[DEP({X,Y}|Z)]', 'FontSize', 20); grid on;
h3.FontSize = 20;

h4 = subplot(3,2,4);
depIdx = 4;
meanRscdmRes = mean(squeeze(rscdmResultsMat(:,depIdx,:)),2);
meanPdcorrRes = mean(squeeze(pdcorrResultsMat(:,depIdx,:)),2);
meanPcorrRes = mean(squeeze(pcorrResultsMat(:,depIdx,:)),2);
meanCmaRes = mean(squeeze(cmaResultsMat(:,depIdx,:)),2);
plot(gammaVec, meanRscdmRes, 'o-.', ...
     gammaVec, meanPdcorrRes, 'x-.', ...
     gammaVec, meanPcorrRes, 'd-.', ...
     gammaVec, meanCmaRes, '+-.');
xlabel('\gamma', 'FontSize', 20); ylabel('\mu[DEP({X,Y}|Z)]', 'FontSize', 20); grid on;
h4.FontSize = 20;

h5 = subplot(3,2,5);
depIdx = 5;
meanRscdmRes = mean(squeeze(rscdmResultsMat(:,depIdx,:)),2);
meanPdcorrRes = mean(squeeze(pdcorrResultsMat(:,depIdx,:)),2);
meanPcorrRes = mean(squeeze(pcorrResultsMat(:,depIdx,:)),2);
meanCmaRes = mean(squeeze(cmaResultsMat(:,depIdx,:)),2);
plot(gammaVec, meanRscdmRes, 'o-.', ...
     gammaVec, meanPdcorrRes, 'x-.', ...
     gammaVec, meanPcorrRes, 'd-.', ...
     gammaVec, meanCmaRes, '+-.');
xlabel('\gamma', 'FontSize', 20); ylabel('\mu[DEP({X,Y}|Z)]', 'FontSize', 20); grid on;
h5.FontSize = 20; 

h6 = subplot(3,2,6);
depIdx = 6;
meanRscdmRes = mean(squeeze(rscdmResultsMat(:,depIdx,:)),2);
meanPdcorrRes = mean(squeeze(pdcorrResultsMat(:,depIdx,:)),2);
meanPcorrRes = mean(squeeze(pcorrResultsMat(:,depIdx,:)),2);
meanCmaRes = mean(squeeze(cmaResultsMat(:,depIdx,:)),2);
plot(gammaVec, meanRscdmRes, 'o-.', ...
     gammaVec, meanPdcorrRes, 'x-.', ...
     gammaVec, meanPcorrRes, 'd-.', ...
     gammaVec, meanCmaRes, '+-.');
xlabel('\gamma', 'FontSize', 20); ylabel('\mu[DEP({X,Y}|Z)]', 'FontSize', 20); grid on;
h6.FontSize = 20; 

legend('RSCDM', 'PDCORR', 'PCORR', 'CMA');  % manually move this using the mouse to a
                                            % good location

%% Conditional Independence Test

clear;
clc;
dbstop if error;

rng(12345);

M = 500;
nsim = 500;

gammaVec = 0:0.1:1;

numDepTypes = 6;
rscdmResultsMat = zeros(length(gammaVec), numDepTypes);
cmaResultsMat = zeros(length(gammaVec), numDepTypes);
hdResultsMat = zeros(length(gammaVec), numDepTypes);
pdcorrResultsMat = zeros(length(gammaVec), numDepTypes);
pcorrResultsMat = zeros(length(gammaVec), numDepTypes);

rscdmResultsVec = zeros(1,nsim);
cmaResultsVec = zeros(1,nsim);
hdResultsVec = zeros(1,nsim);
pdcorrResultsVec = zeros(1,nsim);
pcorrResultsVec = zeros(1,nsim);

dispstat('','init'); % One time only initialization
dispstat(sprintf('Begining the simulation...\n'),'keepthis','timestamp');
for gammaIdx=1:length(gammaVec)
    gamma = gammaVec(gammaIdx);
    dispstat(sprintf('Computing for gamma=%0.02f',gamma),'keepthis', 'timestamp');
    for jj=1:numDepTypes
        for ii=1:nsim
            dispstat(sprintf('Simulating -- %0.02f %%', ii/nsim*100),'timestamp');
            
            X = rand(M,1);
            eps = randn(M,1);
            
            switch(jj)
                case 1
                    Y = gamma*X + (1-gamma)*eps;
                    Z = gamma*X + (1-gamma)*eps;
                case 2
                    Y = gamma*(X-0.5).^2 + (1-gamma)*eps;
                    Z = gamma*(X-0.5).^2 + (1-gamma)*eps;
                case 3
                    Y = gamma*sin(4*pi*X) + (1-gamma)*eps;
                    Z = gamma*cos(4*pi*X) + (1-gamma)*eps;
                case 4
                    Y = gamma*nthroot(X,4) + (1-gamma)*eps;
                    Z = gamma*nthroot(X,4) + (1-gamma)*eps;
                case 5
                    Y = gamma*X + (1-gamma)*eps;
                    Z = gamma*(X-0.5).^2 + (1-gamma)*eps;
                case 6
                    Y = gamma*(X-0.5).^2 + (1-gamma)*eps;
                    Z = gamma*cos(4*pi*X) + (1-gamma)*eps;
            end
            
            rscdmVal = rscdm(Y,Z,X);
            data = struct();
            data.X = Y; data.Y = Z; data.Z = X; cmaVal = cassor(data);
            hdVal = hd(Y,Z,X);
            pdcorrVal = pdcorr_R(Y, Z, X);
            pcorrVal = partialcorr(Y, Z, X);
            
            rscdmResultsVec(ii) = rscdmVal;
            cmaResultsVec(ii) = cmaVal;
            hdResultsVec(ii) = hdVal;
            pdcorrResultsVec(ii) = pdcorrVal;
            pcorrResultsVec(ii) = pcorrVal;
            
        end
        
        rscdmMean = mean(rscdmResultsVec); rscdmSTD = std(rscdmResultsVec);
        cmaMean = mean(cmaResultsVec); cmaSTD = std(cmaResultsVec);
        hdMean = mean(hdResultsVec); hdSTD = std(hdResultsVec);
        pdcorrMean = mean(pdcorrResultsVec); pdcorrSTD = std(pdcorrResultsVec);
        pcorrMean = mean(pcorrResultsVec); pcorrSTD = std(pcorrResultsVec);
        
        rscdmResultsMat(gammaIdx, jj) = rscdmMean;
        cmaResultsMat(gammaIdx, jj) = cmaMean;
        hdResultsMat(gammaIdx, jj) = hdMean;
        pdcorrResultsMat(gammaIdx, jj) = pdcorrMean;
        pcorrResultsMat(gammaIdx, jj) = pcorrMean;
        
    end
end

% save the results
if(ispc)
    save('C:\\Users\\Kiran\\ownCloud\\PhD\\sim_results\\independence\\rscdm_CI.mat');
elseif(ismac)
    save('/Users/Kiran/ownCloud/PhD/sim_results/independence/rscdm_CI.mat');
elseif(isunix)
    save('/home/kiran/ownCloud/PhD/sim_results/independence/rscdm_CI.mat');
end

% plot the dependence metric results versus gamma for each dep type
figure;

h1 = subplot(2,2,1);
plot(gammaVec, rscdmResultsMat(:,1), 'o-.', ...
     gammaVec, pdcorrResultsMat(:,1), 'x-.', ...
     gammaVec, pcorrResultsMat(:,1), 'd-.', ...
     gammaVec, cmaResultsMat(:,1), '+-.', ...
     gammaVec, hdResultsMat(:,1),  '^-');
xlabel('\gamma', 'FontSize', '20'); ylabel('DEP({X,Y}|Z)'); grid on;
 
h2 = subplot(2,2,2);
plot(gammaVec, rscdmResultsMat(:,2), 'o-.', ...
     gammaVec, pdcorrResultsMat(:,2), 'x-.', ...
     gammaVec, pcorrResultsMat(:,2), 'd-.', ...
     gammaVec, cmaResultsMat(:,2), '+-.', ...
     gammaVec, hdResultsMat(:,2),  '^-');
xlabel('\gamma', 'FontSize', '20'); ylabel('DEP({X,Y}|Z)'); grid on;

h3 = subplot(2,2,3);
plot(gammaVec, rscdmResultsMat(:,3), 'o-.', ...
     gammaVec, pdcorrResultsMat(:,3), 'x-.', ...
     gammaVec, pcorrResultsMat(:,3), 'd-.', ...
     gammaVec, cmaResultsMat(:,3), '+-.', ...
     gammaVec, hdResultsMat(:,3),  '^-');
xlabel('\gamma', 'FontSize', '20'); ylabel('DEP({X,Y}|Z)'); grid on;

h4 = subplot(2,2,4);
plot(gammaVec, rscdmResultsMat(:,4), 'o-.', ...
     gammaVec, pdcorrResultsMat(:,4), 'x-.', ...
     gammaVec, pcorrResultsMat(:,4), 'd-.', ...
     gammaVec, cmaResultsMat(:,4), '+-.', ...
     gammaVec, hdResultsMat(:,4),  '^-');
xlabel('\gamma', 'FontSize', '20'); ylabel('DEP({X,Y}|Z)'); grid on;

h5 = subplot(3,2,5);
plot(gammaVec, rscdmResultsMat(:,5), 'o-.', ...
     gammaVec, pdcorrResultsMat(:,5), 'x-.', ...
     gammaVec, pcorrResultsMat(:,5), 'd-.', ...
     gammaVec, cmaResultsMat(:,5), '+-.', ...
     gammaVec, hdResultsMat(:,5),  '^-');
xlabel('\gamma', 'FontSize', 20); ylabel('DEP({X,Y}|Z)', 'FontSize', 20); grid on;

h6 = subplot(3,2,6);
plot(gammaVec, rscdmResultsMat(:,6), 'o-.', ...
     gammaVec, pdcorrResultsMat(:,6), 'x-.', ...
     gammaVec, pcorrResultsMat(:,6), 'd-.', ...
     gammaVec, cmaResultsMat(:,6), '+-.', ...
     gammaVec, hdResultsMat(:,6),  '^-');
xlabel('\gamma', 'FontSize', 20); ylabel('DEP({X,Y}|Z)', 'FontSize', 20); grid on;

legend('RSCDM', 'PDCORR', 'PCORR', 'CMA', 'HD');  % manually move this using the mouse to a
                                                  % good location

%% Do a conditionally independent test

clear;
clc;

rng(123);

M = 500;
noise = 0.01;
alpha = 0.05;

% Generate data from     Y<--X-->Z
x = rand(M,1);
% y = x + noise*randn(M,1); z = x + noise*randn(M,1);
% y = 4*(x-0.5).^2 + noise*randn(M,1); z = 4*(x-0.5).^2 + noise*randn(M,1);
% y = sin(4*pi*x) + noise*randn(M,1); z = cos(4*pi*x) + noise*randn(M,1);
% y = nthroot(x, 4) + noise*randn(M,1); z = nthroot(x, 4) + noise*randn(M,1);
% y = x + noise*randn(M,1); z = 4*(x-0.5).^2 + noise*randn(M,1);
% y = 4*(x-0.5).^2 + noise*randn(M,1); z = cos(4*pi*x) + noise*randn(M,1);

rsdm1 = rsdm(x, y);
rsdm2 = rsdm(x, z);
rsdm3 = rsdm(y, z);
% RSCDM conditions X on Y and Z, and sees how related Y and Z
% are to each other ... in this graphical model, they should be UNRELATED
% after the effect of X is removed... i.e. close to independent.  To see
% why, look at the graphical model, Y indep Z | X according to
% d-separation.  So if we condition upon X (i.e. remove teh effect of X on
% Y and Z separately), then we should get independence.
[rscdmVal, RxAligned, RyAligned] = rscdm(y,z,x);
pdCorr_val = abs(pdcorr_R(y,z,x));
partialCorrVal = abs(partialcorr(y,z,x));
data.X = y; data.Y = z; data.Z = x; cassorVal = abs(cassor(data));
hdVal = abs(hd(y, z, x));

fontSize = 10;

subplot(3,12,1:3);
scatter(x,y); grid on; xlabel('x', 'FontSize', fontSize); ylabel('y', 'FontSize', fontSize); 
title(sprintf('RSDM=%0.2f', rsdm1), 'FontSize', fontSize);

subplot(3,12,5:8);
scatter(y,z); grid on; xlabel('y', 'FontSize', fontSize); ylabel('z', 'FontSize', fontSize); 
title(sprintf('RSDM=%0.2f', rsdm3), 'FontSize', fontSize);

subplot(3,12,10:12);
scatter(x,z); grid on; xlabel('x', 'FontSize', fontSize); ylabel('z', 'FontSize', fontSize); 
title(sprintf('RSDM=%0.2f', rsdm2), 'FontSize', fontSize);

subplot(3,12,13:17);
scatter(1:M, RxAligned); grid on;
xlabel('x', 'FontSize', fontSize); ylabel('r_y', 'FontSize', fontSize);

subplot(3,12,20:24);
scatter(1:M, RyAligned); grid on;
xlabel('x', 'FontSize', fontSize); ylabel('r_z', 'FontSize', fontSize);

subplot(3,12,25:29);
scatter(RxAligned,RyAligned); grid on; 
xlabel('r_y', 'FontSize', fontSize); ylabel('r_z', 'FontSize', fontSize);  
title(sprintf('%0.02f/%0.02f/%0.02f/%0.02f/%0.02f', ...
    rscdmVal, pdCorr_val, partialCorrVal, cassorVal, hdVal), 'FontSize', fontSize);

subplot(3,12,32:36);
scatter(pobs(RxAligned),pobs(RyAligned), 'r'); grid on; 
xlabel('F_{r_y}', 'FontSize', fontSize); ylabel('F_{r_z}', 'FontSize', fontSize);


%% Do a conditionally dependent test

% Tests Conditional Independence w/ RSDM and the residuals processing
% algorithm.

clear;
clc;

% rng(123);

M = 500;
noise = 0;
alpha = 0.05;

% Generate data from     Y-->X<--Z
y = rand(M,1);
z = rand(M,1);
% x = y + z + noise*randn(M,1);
% x = (y-0.5).^2 + (z-0.5).^2 + noise*randn(M,1);
% x = sin(4*pi*y)+cos(4*pi*z) + noise*randn(M,1);
% x = nthroot(y,4)+nthroot(z,4) + noise*randn(M,1);
% x = y + (z-0.5).^2 + noise*randn(M,1);
% x = (y-0.5).^2 + cos(4*pi*z) + noise*randn(M,1);

rsdm1 = rsdm(x, y);
rsdm2 = rsdm(x, z);
rsdm3 = rsdm(y, z);
% In this graphical model, Y and Z are independent of each other, but when
% conditioned upon X, they become dependent.  Refer to the rules of
% d-separation to see why, this is a V-Structure!
[rscdmVal, RxAligned, RyAligned] = rscdm(y,z,x);
pdCorr_val = abs(pdcorr_R(y,z,x));
partialCorrVal = abs(partialcorr(y,z,x));
data.X = y; data.Y = z; data.Z = x; cassorVal = abs(cassor(data));
hdVal = abs(hd(y, z, x));

fontSize = 10;

subplot(3,12,1:3);
scatter(x,y); grid on; xlabel('x', 'FontSize', fontSize); ylabel('y', 'FontSize', fontSize); 
title(sprintf('RSDM=%0.2f', rsdm1), 'FontSize', fontSize);

subplot(3,12,5:8);
scatter(y,z); grid on; xlabel('y', 'FontSize', fontSize); ylabel('z', 'FontSize', fontSize); 
title(sprintf('RSDM=%0.2f', rsdm3), 'FontSize', fontSize);

subplot(3,12,10:12);
scatter(x,z); grid on; xlabel('x', 'FontSize', fontSize); ylabel('z', 'FontSize', fontSize); 
title(sprintf('RSDM=%0.2f', rsdm2), 'FontSize', fontSize);

subplot(3,12,13:17);
scatter(1:M, RxAligned); grid on;
xlabel('x', 'FontSize', fontSize); ylabel('r_y', 'FontSize', fontSize);

subplot(3,12,20:24);
scatter(1:M, RyAligned); grid on;
xlabel('x', 'FontSize', fontSize); ylabel('r_z', 'FontSize', fontSize);

subplot(3,12,25:29);
scatter(RxAligned,RyAligned); grid on; 
xlabel('r_y', 'FontSize', fontSize); ylabel('r_z', 'FontSize', fontSize);  
title(sprintf('%0.02f/%0.02f/%0.02f/%0.02f/%0.02f', ...
    rscdmVal, pdCorr_val, partialCorrVal, cassorVal, hdVal), 'FontSize', fontSize);

subplot(3,12,32:36);
scatter(pobs(RxAligned),pobs(RyAligned), 'r'); grid on; 
xlabel('F_{r_y}', 'FontSize', fontSize); ylabel('F_{r_z}', 'FontSize', fontSize);

%% Characterize null distribution {Y indep Z} | X
