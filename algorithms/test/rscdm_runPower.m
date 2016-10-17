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

rng(12345);

M = 500;
nsim = 500;

numDepTypes = 6;
rscdmResultsMat = zeros(length(gammaVec), numDepTypes);
cmaResultsMat = zeros(length(gammaVec), numDepTypes);
cmiResultsMat = zeros(length(gammaVec), numDepTypes);
hdResultsMat = zeros(length(gammaVec), numDepTypes);
hsncicResultsMat = zeros(length(gammaVec), numDepTypes);

rscdmResultsVec = zeros(1,nsim);
cmaResultsVec = zeros(1,nsim);
cmiResultsVec = zeros(1,nsim);
hdResultsVec = zeros(1,nsim);
hsncicResultsVec = zeros(1,nsim);

gammaVec = 0:0.1:1;
for gammaIdx=1:length(gammaVec)
    gamma = gammaVec(gammaIdx);
    for jj=1:numDepTypes
        for ii=1:nsim
            Y = rand(M,1);
            Z = rand(M,1);
            eps = randn(M,1);
            
            switch(jj)
                case 1
                    X = gamma*(Y+Z) + (1-gamma)*eps;
                case 2
                    X = gamma*((Y-0.5).^2 + (Z-0.5).^2) + (1-gamma)*eps;
                case 3
                    X = gamma*(sin(4*pi*Y)+cos(4*pi*Y)) + (1-gamma)*eps;
                case 4
                    X = gamma*(nthroot(Y,4)+nthroot(Z,4)) + (1-gamma)*eps;
                case 5
                    X = gamma*(Y + (Z-0.5).^2) + (1-gamma)*eps;
                case 6
                    X = gamma*((Y-0.5).^2 + cos(4*pi*Z)) + (1-gamma)*eps;
            end
            
            rscdmVal = rscdm(y,z,x);
            data.X = Y; data.Y = Z; data.Z = X;
            cmaVal = cassor(data);
            cmiVal = cmi(data);
            hdVal = hd(Y,Z,X);
            hsncicVal = hsncic(Y,Z,X);
            
            rscdmResultsVec(ii) = rscdmVal;
            cmaResultsVec(ii) = cmaVal;
            cmiResultsVec(ii) = cmiVal;
            hdResultsVec(ii) = hdVal;
            hsncicResultsVec(ii) = hsncicVal;
            
        end
        
        rscdmMean = mean(rscdmResultsVec); rscdmSTD = std(rscdmResultsVec);
        cmaMean = mean(cmaResultsVec); cmaSTD = std(cmaResultsVec);
        cmiMean = mean(cmiResultsVec); cmiSTD = std(cmiResultsVec);
        hdMean = mean(hdResultsVec); hdSTD = std(hdResultsVec);
        hsncicMean = mean(hsncicResultsVec); hsncicSTD = std(hsncicResultsVec);
        
        rscdmResultsMat(gammaIdx, jj) = rscdmMean;
        cmaResultsMat(gammaIdx, jj) = cmaMean;
        cmiResultsMat(gammaIdx, jj) = cmiMean;
        hdResultsMat(gammaIdx, jj) = hdMean;
        hsncicResultsMat(gammaIdx, jj) = hsncicMean;
        
    end
end

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

h1 = subplot(2,2,1);
plot(gammaVec, rscdmResultsMat(:,1), 'o-.', ...
     gammaVec, cmaResultsMat(:,1), '+-.', ...
     gammaVec, cmiResultsMat(:,1), 'd-.', ...
     gammaVec, hdResultsMat(:,1),  '^-');
xlabel('\gamma', 'FontSize', '20'); ylabel('DEP({X,Y}|Z)'); grid on;
 
h2 = subplot(2,2,2);
plot(gammaVec, rscdmResultsMat(:,2), 'o-.', ...
     gammaVec, cmaResultsMat(:,2), '+-.', ...
     gammaVec, cmiResultsMat(:,2), 'd-.', ...
     gammaVec, hdResultsMat(:,2),  '^-');
xlabel('\gamma', 'FontSize', '20'); ylabel('DEP({X,Y}|Z)'); grid on;

h3 = subplot(2,2,3);
plot(gammaVec, rscdmResultsMat(:,3), 'o-.', ...
     gammaVec, cmaResultsMat(:,3), '+-.', ...
     gammaVec, cmiResultsMat(:,3), 'd-.', ...
     gammaVec, hdResultsMat(:,3),  '^-');
xlabel('\gamma', 'FontSize', '20'); ylabel('DEP({X,Y}|Z)'); grid on;

h4 = subplot(2,2,4);
plot(gammaVec, rscdmResultsMat(:,4), 'o-.', ...
     gammaVec, cmaResultsMat(:,4), '+-.', ...
     gammaVec, cmiResultsMat(:,4), 'd-.', ...
     gammaVec, hdResultsMat(:,4),  '^-');
xlabel('\gamma', 'FontSize', '20'); ylabel('DEP({X,Y}|Z)'); grid on;

%% Conditional Independence Test

clear;
clc;

rng(12345);

M = 500;
nsim = 500;

numDepTypes = 6;
rscdmResultsMat = zeros(length(gammaVec), numDepTypes);
cmaResultsMat = zeros(length(gammaVec), numDepTypes);
cmiResultsMat = zeros(length(gammaVec), numDepTypes);
hdResultsMat = zeros(length(gammaVec), numDepTypes);
hsncicResultsMat = zeros(length(gammaVec), numDepTypes);

rscdmResultsVec = zeros(1,nsim);
cmaResultsVec = zeros(1,nsim);
cmiResultsVec = zeros(1,nsim);
hdResultsVec = zeros(1,nsim);
hsncicResultsVec = zeros(1,nsim);

gammaVec = 0:0.1:1;
for gammaIdx=1:length(gammaVec)
    gamma = gammaVec(gammaIdx);
    for jj=1:numDepTypes
        for ii=1:nsim
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
            
            rscdmVal = rscdm(y,z,x);
            data.X = Y; data.Y = Z; data.Z = X;
            cmaVal = cassor(data);
            cmiVal = cmi(data);
            hdVal = hd(Y,Z,X);
            hsncicVal = hsncic(Y,Z,X);
            
            rscdmResultsVec(ii) = rscdmVal;
            cmaResultsVec(ii) = cmaVal;
            cmiResultsVec(ii) = cmiVal;
            hdResultsVec(ii) = hdVal;
            hsncicResultsVec(ii) = hsncicVal;
            
        end
        
        rscdmMean = mean(rscdmResultsVec); rscdmSTD = std(rscdmResultsVec);
        cmaMean = mean(cmaResultsVec); cmaSTD = std(cmaResultsVec);
        cmiMean = mean(cmiResultsVec); cmiSTD = std(cmiResultsVec);
        hdMean = mean(hdResultsVec); hdSTD = std(hdResultsVec);
        hsncicMean = mean(hsncicResultsVec); hsncicSTD = std(hsncicResultsVec);
        
        rscdmResultsMat(gammaIdx, jj) = rscdmMean;
        cmaResultsMat(gammaIdx, jj) = cmaMean;
        cmiResultsMat(gammaIdx, jj) = cmiMean;
        hdResultsMat(gammaIdx, jj) = hdMean;
        hsncicResultsMat(gammaIdx, jj) = hsncicMean;
        
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
     gammaVec, cmaResultsMat(:,1), '+-.', ...
     gammaVec, cmiResultsMat(:,1), 'd-.', ...
     gammaVec, hdResultsMat(:,1),  '^-');
xlabel('\gamma', 'FontSize', '20'); ylabel('DEP({X,Y}|Z)'); grid on;
 
h2 = subplot(2,2,2);
plot(gammaVec, rscdmResultsMat(:,2), 'o-.', ...
     gammaVec, cmaResultsMat(:,2), '+-.', ...
     gammaVec, cmiResultsMat(:,2), 'd-.', ...
     gammaVec, hdResultsMat(:,2),  '^-');
xlabel('\gamma', 'FontSize', '20'); ylabel('DEP({X,Y}|Z)'); grid on;

h3 = subplot(2,2,3);
plot(gammaVec, rscdmResultsMat(:,3), 'o-.', ...
     gammaVec, cmaResultsMat(:,3), '+-.', ...
     gammaVec, cmiResultsMat(:,3), 'd-.', ...
     gammaVec, hdResultsMat(:,3),  '^-');
xlabel('\gamma', 'FontSize', '20'); ylabel('DEP({X,Y}|Z)'); grid on;

h4 = subplot(2,2,4);
plot(gammaVec, rscdmResultsMat(:,4), 'o-.', ...
     gammaVec, cmaResultsMat(:,4), '+-.', ...
     gammaVec, cmiResultsMat(:,4), 'd-.', ...
     gammaVec, hdResultsMat(:,4),  '^-');
xlabel('\gamma', 'FontSize', '20'); ylabel('DEP({X,Y}|Z)'); grid on;

%% Do a conditionally independent test

clear;
clc;

rng(123);

M = 500;
noise = 0.01;
alpha = 0.05;

% Generate data from     Y<--X-->Z
x = rand(M,1);
% y = 4*(x-0.5).^2 + noise*randn(M,1);
% z = 4*(x-0.5).^2 + noise*randn(M,1);
% y = x + noise*randn(M,1);
% z = x + noise*randn(M,1);
y = sin(2*pi*x) + noise*randn(M,1);
z = 4*(x-0.5).^2 + noise*randn(M,1);

rsdm1 = rsdm(x, y);
rsdm2 = rsdm(x, z);
rsdm3 = rsdm(y, z);
% RSCDM conditions X on Y and Z, and sees how related Y and Z
% are to each other ... in this graphical model, they should be UNRELATED
% after the effect of X is removed... i.e. close to independent.  To see
% why, look at the graphical model, Y indep Z | X according to
% d-separation.  So if we condition upon X (i.e. remove teh effect of X on
% Y and Z separately), then we should get independence.
[rscdmVal, RxStacked, RyStacked, RxPtsStacked, RyPtsStacked] = rscdm(y,z,x);
rscdmVal_Z = sqrt(M-4)*abs(rscdmVal);
testVal = norminv(1-alpha/2);

subplot(3,12,1:3);
scatter(x,y); grid on; xlabel('x', 'FontSize', 20); ylabel('y', 'FontSize', 20); 
title(sprintf('RSDM=%0.2f', rsdm1), 'FontSize', 20);

subplot(3,12,5:8);
scatter(y,z); grid on; xlabel('y', 'FontSize', 20); ylabel('z', 'FontSize', 20); 
title(sprintf('RSDM=%0.2f', rsdm3), 'FontSize', 20);

subplot(3,12,10:12);
scatter(x,z); grid on; xlabel('x', 'FontSize', 20); ylabel('z', 'FontSize', 20); 
title(sprintf('RSDM=%0.2f', rsdm2), 'FontSize', 20);

subplot(3,12,13:17);
scatter(RxPtsStacked, RxStacked); grid on;
xlabel('x', 'FontSize', 20); ylabel('r_y', 'FontSize', 20);

subplot(3,12,20:24);
scatter(RyPtsStacked, RyStacked); grid on;
xlabel('x', 'FontSize', 20); ylabel('r_z', 'FontSize', 20);

subplot(3,12,25:29);
scatter(RxStacked,RyStacked); grid on; 
xlabel('r_y', 'FontSize', 20); ylabel('r_z', 'FontSize', 20);  
title(sprintf('rscdm=%0.2f', rscdmVal), 'FontSize', 20);

subplot(3,12,32:36);
scatter(pobs(RxStacked),pobs(RyStacked), 'r'); grid on; 
xlabel('F_{r_y}', 'FontSize', 20); ylabel('F_{r_z}', 'FontSize', 20);


%% Do a conditionally dependent test

% Tests Conditional Independence w/ RSDM and the residuals processing
% algorithm.

clear;
clc;

rng(123);

M = 500;
noise = 0.5;
alpha = 0.05;

% Generate data from     Y-->X<--Z
y = rand(M,1);
z = rand(M,1);
x = y.^2+z.^3 + noise;      % TODO: fix the bug :D
% x = y + z + noise;


rsdm1 = rsdm(x, y);
rsdm2 = rsdm(x, z);
rsdm3 = rsdm(y, z);
% In this graphical model, Y and Z are independent of each other, but when
% conditioned upon X, they become dependent.  Refer to the rules of
% d-separation to see why, this is a V-Structure!
[rscdmVal, RxStacked, RyStacked, RxPtsStacked, RyPtsStacked] = rscdm(y,z,x);
rscdmVal_Z = sqrt(M-4)*abs(rscdmVal);
testVal = norminv(1-alpha/2);

subplot(3,12,1:3);
scatter(x,y); grid on; xlabel('x', 'FontSize', 20); ylabel('y', 'FontSize', 20); 
title(sprintf('RSDM=%0.2f', rsdm1), 'FontSize', 20);

subplot(3,12,5:8);
scatter(y,z); grid on; xlabel('y', 'FontSize', 20); ylabel('z', 'FontSize', 20); 
title(sprintf('RSDM=%0.2f', rsdm3), 'FontSize', 20);

subplot(3,12,10:12);
scatter(x,z); grid on; xlabel('x', 'FontSize', 20); ylabel('z', 'FontSize', 20); 
title(sprintf('RSDM=%0.2f', rsdm2), 'FontSize', 20);

subplot(3,12,13:17);
scatter(RxPtsStacked, RxStacked); grid on;
xlabel('x', 'FontSize', 20); title('r_y', 'FontSize', 20);

subplot(3,12,20:24);
scatter(RyPtsStacked, RyStacked); grid on;
xlabel('x', 'FontSize', 20); title('r_z', 'FontSize', 20);

subplot(3,12,25:29);
scatter(RxStacked,RyStacked); grid on; xlabel('r_y', 'FontSize', 20); ylabel('r_z', 'FontSize', 20);  
title(sprintf('rscdm=%0.2f', rscdmVal), 'FontSize', 20);

subplot(3,12,32:36);
scatter(pobs(RxStacked),pobs(RyStacked), 'r'); grid on; xlabel('F_{r_y}'); ylabel('F_{r_z}');

%% Characterize null distribution {Y indep Z} | X
