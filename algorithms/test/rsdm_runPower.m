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

%% Generate statistical power curves for our dependency metric using the
% same methodology as Simon & Tibshirani:
% http://statweb.stanford.edu/~tibs/reshef/script.R

clear;
clc;

% WARNING: ENSURE THAT minepy/matlab/ is in the matlab path for MIC to
% work!

rng(1234);
dbstop if error;

nsim_null = 500;   % The number of null datasets we use to estimate our rejection reject regions for an alternative with level 0.05
nsim_alt  = 500;   % Number of alternative datasets we use to estimate our power

num_noise = 30;                    % The number of different noise levels used
noise = 3;                         % A constant to determine the amount of noise

M = 500;                % number of samples
numDepTests = 8;        % the number of different dependency tests we will conduct
                        % TODO: add copula dependencies as well

% Vectors holding the null "correlations" (for pearson, dcor and mic respectively) 
% for each of the nsim null datasets at a given noise level
rsdmNull = zeros(1,nsim_null);
dcorrNull = zeros(1,nsim_null);
ticeNull = zeros(1,nsim_null);
corrNull = zeros(1,nsim_null);
rdcNull  = zeros(1,nsim_null);
cosNull  = zeros(1,nsim_null);
ccorrNull = zeros(1,nsim_null);

rsdmAlt  = zeros(1,nsim_alt);
dcorrAlt = zeros(1,nsim_alt);
ticeAlt = zeros(1,nsim_alt);
corrAlt = zeros(1,nsim_alt);
rdcAlt  = zeros(1,nsim_alt);
cosAlt  = zeros(1,nsim_alt);
ccorrAlt = zeros(1,nsim_alt);

% Arrays holding the estimated power for each of the "correlation" types, 
% for each data type (linear, parabolic, etc...) with each noise level
rsdmPower = zeros(numDepTests, num_noise);
dcorrPower = zeros(numDepTests, num_noise);
ticePower = zeros(numDepTests, num_noise);
corrPower = zeros(numDepTests, num_noise);
rdcPower  = zeros(numDepTests, num_noise);
cosPower  = zeros(numDepTests, num_noise);
ccorrPower = zeros(numDepTests, num_noise);

% Optimal parameters for MICe
mine_c = 15;
mine_alpha = 0.6;

% Optimal parameters for RDC
rdc_k = 20;
rdc_s = 1/6;

% We loop through the noise level and functional form; 
% each time we estimate a null distribution based on the marginals of the data, 
% and then use that null distribution to estimate power

% We use a uniformly distributed x, because in the original paper the 
% authors used the same

% Simon & Tibshirani use xMin=0, xMax=1 for performing their analysis ...
xMin = 0;
xMax = 1;

dispstat('','init'); % One time only initialization
dispstat(sprintf('Begining the simulation...\n'),'keepthis','timestamp');
num_noise_test_min = 1;
num_noise_test_max = 30;
for l=num_noise_test_min:num_noise_test_max
    for typ=1:numDepTests
        dispstat(sprintf('Computing for noise level=%d Dependency Test=%d',l, typ),'keepthis', 'timestamp');
        % simulate data under the null w/ correct marginals
        parfor ii=1:nsim_null
            x = rand(M,1)*(xMax-xMin)+xMin;
            switch(typ)
                case 1
                    % linear
                    y = x + noise*(l/num_noise)*randn(M,1); 
                case 2
                    % parabolic
                    y = 4*(x-.5).^2 + noise*(l/num_noise)*randn(M,1);
                case 3
                    % cubic
                    y = 128*(x-1/3).^3-48*(x-1/3).^3-12*(x-1/3)+10* noise*(l/num_noise)*randn(M,1);
                case 4
                    % low-freq sin
                    y = sin(4*pi*x) + 2*noise*(l/num_noise)*randn(M,1);
                case 5
                    % high-freq sin
                    y = sin(16*pi*x) + noise*(l/num_noise)*randn(M,1);
                case 6
                    % fourth root
                    y = x.^(1/4) + noise*(l/num_noise)*randn(M,1);
                case 7
                    % circle
                    y=(2*binornd(1,0.5,M,1)-1) .* (sqrt(1 - (2*x - 1).^2)) + noise/4*l/num_noise*randn(M,1);
                case 8
                    % step function
                    y = (x > 0.5) + noise*5*l/num_noise*randn(M,1);
                otherwise
                    error('unknown dep type!');
            end
            % resimulate x so we have null scenario
            x = rand(M,1)*(xMax-xMin)+xMin;
            
            % calculate the metrics
            rsdmNull(ii) = rsdm(x, y);
            dcorrNull(ii) = dcorr(x, y);
            % compute MICe
            minestats = mine(x',y',mine_alpha,mine_c,'mic_e');
            ticeNull(ii) = minestats.tic;
            % compute correlation
            corrNull(ii) = corr(x,y);
            % compute RDC
            rdcNull(ii) = rdc(x,y,rdc_k,rdc_s);
            % compute CoS
            cosNull(ii) = cosf(x,y);
            % compute cCorr
            ccorrNull(ii) = cCorr(x,y);
        end
        
        % compute the rejection cutoffs
        rsdm_cut = quantile(rsdmNull, 0.95);
        dcorr_cut = quantile(dcorrNull, 0.95);
        tice_cut = quantile(ticeNull, 0.95);
        corr_cut = quantile(corrNull, 0.95);
        rdc_cut  = quantile(rdcNull, 0.95);
        cos_cut  = quantile(cosNull, 0.95);
        ccorr_cut = quantile(ccorrNull, 0.95);
        
        % resimulate the data under the alternative hypothesis
        parfor ii=1:nsim_alt
            x = rand(M,1)*(xMax-xMin)+xMin;
            switch(typ)
                case 1
                    % linear
                    y = x + noise*(l/num_noise)*randn(M,1); 
                case 2
                    % parabolic
                    y = 4*(x-.5).^2 + noise*(l/num_noise)*randn(M,1);
                case 3
                    % cubic
                    y = 128*(x-1/3).^3-48*(x-1/3).^3-12*(x-1/3) + 10*noise*(l/num_noise)*randn(M,1);
                case 4
                    % low-freq sin
                    y = sin(4*pi*x) + 2*noise*(l/num_noise)*randn(M,1);
                case 5
                    % high-freq sin
                    y = sin(16*pi*x) + noise*(l/num_noise)*randn(M,1);
                case 6
                    % fourth root
                    y = x.^(1/4) + noise*(l/num_noise)*randn(M,1);
                case 7
                    % circle
                    y=(2*binornd(1,0.5,M,1)-1) .* (sqrt(1 - (2*x - 1).^2)) + noise/4*l/num_noise*randn(M,1);
                case 8
                    % step function
                    y = (x > 0.5) + noise*5*l/num_noise*randn(M,1);
                otherwise
                    error('unknown dep type!');
            end
            
            % calculate the metrics
            rsdmAlt(ii) = rsdm(x, y);
            dcorrAlt(ii) = dcorr(x, y);
            % compute MICe
            minestats = mine(x',y',mine_alpha,mine_c,'mic_e');
            ticeAlt(ii) = minestats.tic;
            % compute corr
            corrAlt(ii) = corr(x, y);
            % compute RDC
            rdcAlt(ii) = rdc(x,y,rdc_k,rdc_s);
            % compute CoS
            cosAlt(ii) = cosf(x,y);
            % compute cCorr
            ccorrAlt(ii) = cCorr(x,y);
        end
        
        % compute the power
        rsdmPower(typ, l)   = sum(rsdmAlt > rsdm_cut)/nsim_alt;
        dcorrPower(typ, l)  = sum(dcorrAlt > dcorr_cut)/nsim_alt;
        ticePower(typ, l)   = sum(ticeAlt > tice_cut)/nsim_alt;
        corrPower(typ, l)   = sum(corrAlt > corr_cut)/nsim_alt;
        rdcPower(typ, l)    = sum(rdcAlt > rdc_cut)/nsim_alt;
        cosPower(typ, l)    = sum(cosAlt > cos_cut)/nsim_alt;
        ccorrPower(typ, l)  = sum(ccorrAlt > ccorr_cut)/nsim_alt;
    end
end

% save the data
if(ispc)
    save(sprintf('C:\\Users\\Kiran\\ownCloud\\PhD\\sim_results\\independence\\power_M_%d.mat', M));
elseif(ismac)
    save(sprintf('/Users/Kiran/ownCloud/PhD/sim_results/independence/power_M_%d.mat', M));
else
    save(sprintf('/home/kiran/ownCloud/PhD/sim_results/independence/power_M_%d.mat', M));
end

% inlet plot configuration
M_inlet = 200;
if(M==500)
    inset_bufX = 0.0005; inset_bufY = 0.002;
else
    inset_bufX = 0.15; inset_bufY = 0.26;
end

inset_width = 0.1; inset_height = 0.08;

noiseVec = (num_noise_test_min:num_noise_test_max)/10;
figure;
h1 = subplot(2,2,1);
hh1 = plot(noiseVec, rsdmPower(1,num_noise_test_min:num_noise_test_max), 'o-.', ...
     noiseVec, dcorrPower(1,num_noise_test_min:num_noise_test_max), '+-.', ...
     noiseVec, ticePower(1,num_noise_test_min:num_noise_test_max), 'd-.', ...
     noiseVec, rdcPower(1,num_noise_test_min:num_noise_test_max), 'v-.', ...
     noiseVec, cosPower(1,num_noise_test_min:num_noise_test_max), 's-.', ...
     noiseVec, ccorrPower(1,num_noise_test_min:num_noise_test_max), 'p-.'); 
axis([min(noiseVec) max(noiseVec) 0 1]);
xlabel('Noise Level', 'FontSize', 20); ylabel('Power', 'FontSize', 20); grid on;
h1.FontSize = 20; 
loc_inset = [h1.Position(1)+inset_bufX h1.Position(2)+inset_bufY inset_width inset_height];
ax1 = axes('Position',loc_inset);
tmp1 = linspace(0,1,M_inlet);
tmp2 = tmp1;
plot(tmp1,tmp2, 'k', 'LineWidth', 2);
ax1.Box = 'on'; ax1.XTick = []; ax1.YTick = [];
ax1.XLim = [min(tmp1) max(tmp1)];
ax1.YLim = [min(tmp2) max(tmp2)];
hh1(1).LineWidth = 5; 
hh1(2).LineWidth = 1.5; 
hh1(3).LineWidth = 1.5; 
hh1(4).LineWidth = 1.5; 
hh1(5).LineWidth = 1.5;
hh1(6).LineWidth = 1.5;

h2 = subplot(2,2,2);
hh2 = plot(noiseVec, rsdmPower(2,num_noise_test_min:num_noise_test_max), 'o-.', ...
     noiseVec, dcorrPower(2,num_noise_test_min:num_noise_test_max), '+-.', ...
     noiseVec, ticePower(2,num_noise_test_min:num_noise_test_max), 'd-.', ...
     noiseVec, rdcPower(2,num_noise_test_min:num_noise_test_max), 'v-.', ...
     noiseVec, cosPower(2,num_noise_test_min:num_noise_test_max), 's-.', ...
     noiseVec, ccorrPower(2,num_noise_test_min:num_noise_test_max), 'p-.'); 
axis([min(noiseVec) max(noiseVec) 0 1]);
xlabel('Noise Level', 'FontSize', 20); ylabel('Power', 'FontSize', 20); grid on;
h2.FontSize = 20; 
loc_inset = [h2.Position(1)+inset_bufX h2.Position(2)+inset_bufY inset_width inset_height];
ax2 = axes('Position',loc_inset);
tmp1 = linspace(0,1,M_inlet);
tmp2 = 4*(tmp1-.5).^2;
plot(tmp1,tmp2, 'k', 'LineWidth', 2);
ax2.Box = 'on'; ax2.XTick = []; ax2.YTick = [];
ax2.XLim = [min(tmp1) max(tmp1)];
ax2.YLim = [min(tmp2) max(tmp2)];
hh2(1).LineWidth = 5; 
hh2(2).LineWidth = 1.5; 
hh2(3).LineWidth = 1.5; 
hh2(4).LineWidth = 1.5; 
hh2(5).LineWidth = 1.5;
hh2(6).LineWidth = 1.5;

h3 = subplot(2,2,3); 
hh3 = plot(noiseVec, rsdmPower(3,num_noise_test_min:num_noise_test_max), 'o-.', ...
     noiseVec, dcorrPower(3,num_noise_test_min:num_noise_test_max), '+-.', ...
     noiseVec, ticePower(3,num_noise_test_min:num_noise_test_max), 'd-.', ...
     noiseVec, rdcPower(3,num_noise_test_min:num_noise_test_max), 'v-.', ...
     noiseVec, cosPower(3,num_noise_test_min:num_noise_test_max), 's-.', ...
     noiseVec, ccorrPower(3,num_noise_test_min:num_noise_test_max), 'p-.');   
axis([min(noiseVec) max(noiseVec) 0 1]);
xlabel('Noise Level', 'FontSize', 20); ylabel('Power', 'FontSize', 20); grid on;
h3.FontSize = 20; 
loc_inset = [h3.Position(1)+inset_bufX h3.Position(2)+inset_bufY inset_width inset_height];
ax3 = axes('Position',loc_inset);
tmp1 = linspace(0,1,M_inlet);
tmp2 = 128*(tmp1-1/3).^3-48*(tmp1-1/3).^3-12*(tmp1-1/3);
plot(tmp1,tmp2, 'k', 'LineWidth', 2);
ax3.Box = 'on'; ax3.XTick = []; ax3.YTick = [];
ax3.XLim = [min(tmp1) max(tmp1)];
ax3.YLim = [min(tmp2) max(tmp2)];
hh3(1).LineWidth = 5; 
hh3(2).LineWidth = 1.5; 
hh3(3).LineWidth = 1.5; 
hh3(4).LineWidth = 1.5; 
hh3(5).LineWidth = 1.5;
hh3(6).LineWidth = 1.5;

h4 = subplot(2,2,4); 
hh4 = plot(noiseVec, rsdmPower(4,num_noise_test_min:num_noise_test_max), 'o-.', ...
     noiseVec, dcorrPower(4,num_noise_test_min:num_noise_test_max), '+-.', ...
     noiseVec, ticePower(4,num_noise_test_min:num_noise_test_max), 'd-.', ...
     noiseVec, rdcPower(4,num_noise_test_min:num_noise_test_max), 'v-.', ...
     noiseVec, cosPower(4,num_noise_test_min:num_noise_test_max), 's-.', ...
     noiseVec, ccorrPower(4,num_noise_test_min:num_noise_test_max), 'p-.'); 
axis([min(noiseVec) max(noiseVec) 0 1]);
xlabel('Noise Level', 'FontSize', 20); ylabel('Power', 'FontSize', 20); grid on;
h4.FontSize = 20; 
loc_inset = [h4.Position(1)+inset_bufX h4.Position(2)+inset_bufY inset_width inset_height];
ax4 = axes('Position',loc_inset);
tmp1 = linspace(0,1,M_inlet);
tmp2 = sin(4*pi*tmp1);
plot(tmp1,tmp2, 'k', 'LineWidth', 2);
ax4.Box = 'on'; ax4.XTick = []; ax4.YTick = [];
ax4.XLim = [min(tmp1) max(tmp1)];
ax4.YLim = [min(tmp2) max(tmp2)];
hh4(1).LineWidth = 5; 
hh4(2).LineWidth = 1.5; 
hh4(3).LineWidth = 1.5; 
hh4(4).LineWidth = 1.5; 
hh4(5).LineWidth = 1.5;
hh4(6).LineWidth = 1.5;

figure;
h5 = subplot(2,2,1); 
hh5 = plot(noiseVec, rsdmPower(5,num_noise_test_min:num_noise_test_max), 'o-.', ...
     noiseVec, dcorrPower(5,num_noise_test_min:num_noise_test_max), '+-.', ...
     noiseVec, ticePower(5,num_noise_test_min:num_noise_test_max), 'd-.', ...
     noiseVec, rdcPower(5,num_noise_test_min:num_noise_test_max), 'v-.', ...
     noiseVec, cosPower(5,num_noise_test_min:num_noise_test_max), 's-.', ...
     noiseVec, ccorrPower(5,num_noise_test_min:num_noise_test_max), 'p-.'); 
axis([min(noiseVec) max(noiseVec) 0 1]);
xlabel('Noise Level'); ylabel('Power'); grid on;
h5.FontSize = 20; 
loc_inset = [h5.Position(1)+inset_bufX h5.Position(2)+inset_bufY inset_width inset_height];
ax5 = axes('Position',loc_inset);
tmp1 = linspace(0,1,M_inlet);
tmp2 = sin(16*pi*tmp1);
plot(tmp1,tmp2, 'k', 'LineWidth', 2);
ax5.Box = 'on'; ax5.XTick = []; ax5.YTick = [];
ax5.XLim = [min(tmp1) max(tmp1)];
ax5.YLim = [min(tmp2) max(tmp2)];
hh5(1).LineWidth = 5; 
hh5(2).LineWidth = 1.5; 
hh5(3).LineWidth = 1.5; 
hh5(4).LineWidth = 1.5; 
hh5(5).LineWidth = 1.5;
hh5(6).LineWidth = 1.5;

h6 = subplot(2,2,2); 
hh6 = plot(noiseVec, rsdmPower(6,num_noise_test_min:num_noise_test_max), 'o-.', ...
     noiseVec, dcorrPower(6,num_noise_test_min:num_noise_test_max), '+-.', ...
     noiseVec, ticePower(6,num_noise_test_min:num_noise_test_max), 'd-.', ...
     noiseVec, rdcPower(6,num_noise_test_min:num_noise_test_max), 'v-.', ...
     noiseVec, cosPower(6,num_noise_test_min:num_noise_test_max), 's-.', ...
     noiseVec, ccorrPower(6,num_noise_test_min:num_noise_test_max), 'p-.'); 
axis([min(noiseVec) max(noiseVec) 0 1]);
xlabel('Noise Level'); ylabel('Power'); grid on;
h6.FontSize = 20; 
loc_inset = [h6.Position(1)+inset_bufX h6.Position(2)+inset_bufY inset_width inset_height];
ax6 = axes('Position',loc_inset);
tmp1 = linspace(0,1,M_inlet);
tmp2 = tmp1.^(1/4);
plot(tmp1,tmp2, 'k', 'LineWidth', 2);
ax6.Box = 'on'; ax6.XTick = []; ax6.YTick = [];
ax6.XLim = [min(tmp1) max(tmp1)];
ax6.YLim = [min(tmp2) max(tmp2)];
hh6(1).LineWidth = 5; 
hh6(2).LineWidth = 1.5; 
hh6(3).LineWidth = 1.5; 
hh6(4).LineWidth = 1.5; 
hh6(5).LineWidth = 1.5;
hh6(6).LineWidth = 1.5;

h7 = subplot(2,2,3); 
hh7 = plot(noiseVec, rsdmPower(7,num_noise_test_min:num_noise_test_max), 'o-.', ...
     noiseVec, dcorrPower(7,num_noise_test_min:num_noise_test_max), '+-.', ...
     noiseVec, ticePower(7,num_noise_test_min:num_noise_test_max), 'd-.', ...
     noiseVec, rdcPower(7,num_noise_test_min:num_noise_test_max), 'v-.', ...
     noiseVec, cosPower(7,num_noise_test_min:num_noise_test_max), 's-.', ...
     noiseVec, ccorrPower(7,num_noise_test_min:num_noise_test_max), 'p-.'); 
axis([min(noiseVec) max(noiseVec) 0 1]);
xlabel('Noise Level'); ylabel('Power'); grid on;
h7.FontSize = 20; 
loc_inset = [h7.Position(1)+inset_bufX h7.Position(2)+inset_bufY inset_width inset_height];
ax7 = axes('Position',loc_inset);
tmp1 = linspace(0,1,M_inlet/2);
tmp2 = (sqrt(1 - (2*tmp1 - 1).^2));
tmp3 = -(sqrt(1 - (2*tmp1 - 1).^2));
plot(tmp1,tmp2, 'k', 'LineWidth', 2); hold on;
plot(tmp1,tmp3, 'k', 'LineWidth', 2); 
ax7.Box = 'on'; ax7.XTick = []; ax7.YTick = [];
ax7.XLim = [min(tmp1) max(tmp1)];
ax7.YLim = [min(tmp3) max(tmp2)];
hh7(1).LineWidth = 5; 
hh7(2).LineWidth = 1.5; 
hh7(3).LineWidth = 1.5; 
hh7(4).LineWidth = 1.5; 
hh7(5).LineWidth = 1.5;
hh7(6).LineWidth = 1.5;

h8 = subplot(2,2,4); 
hh8 = plot(noiseVec, rsdmPower(8,num_noise_test_min:num_noise_test_max), 'o-.', ...
     noiseVec, dcorrPower(8,num_noise_test_min:num_noise_test_max), '+-.', ...
     noiseVec, ticePower(8,num_noise_test_min:num_noise_test_max), 'd-.', ...
     noiseVec, rdcPower(8,num_noise_test_min:num_noise_test_max), 'v-.', ...
     noiseVec, cosPower(8,num_noise_test_min:num_noise_test_max), 's-.', ...
     noiseVec, ccorrPower(8,num_noise_test_min:num_noise_test_max), 'p-.'); 
axis([min(noiseVec) max(noiseVec) 0 1]);
h8.FontSize = 20; 
legend('RSDM', 'dCor', 'TIC_e', 'RDC', 'CoS', 'cCor');  % manually move this using the mouse to a
                                                        % good location
xlabel('Noise Level'); ylabel('Power'); grid on;
loc_inset = [h8.Position(1)+inset_bufX h8.Position(2)+inset_bufY inset_width inset_height];
ax8 = axes('Position',loc_inset);
tmp1 = linspace(0,1,M_inlet);
tmp2 = (tmp1 > 0.5);
plot(tmp1,tmp2, 'k', 'LineWidth', 2);
ax8.Box = 'on'; ax8.XTick = []; ax8.YTick = [];
ax8.XLim = [min(tmp1) max(tmp1)];
hh8(1).LineWidth = 5; 
hh8(2).LineWidth = 1.5; 
hh8(3).LineWidth = 1.5; 
hh8(4).LineWidth = 1.5; 
hh8(5).LineWidth = 1.5;
hh8(6).LineWidth = 1.5;

%% Generate curves which show the effect of sample size for statistical 
% power for the various dependency metrics

clear;
clc;

% WARNING: ENSURE THAT minepy/matlab/ is in the matlab path for MIC to
% work!

rng(1234);
dbstop if error;

nsim_null = 500;   % The number of null datasets we use to estimate our rejection reject regions for an alternative with level 0.05
nsim_alt  = 500;   % Number of alternative datasets we use to estimate our power

num_noise = 30;                    % The number of different noise levels used
noise = 3;                         % A constant to determine the amount of noise

M_vec = 25:25:1500;      % number of samples
numDepTests = 8;        % the number of different dependency tests we will conduct
                        % TODO: add copula dependencies as well

% Vectors holding the null "correlations" (for pearson, dcor and mic respectively) 
% for each of the nsim null datasets at a given noise level
rsdmNull = zeros(1,nsim_null);
dcorrNull = zeros(1,nsim_null);
ticeNull = zeros(1,nsim_null);
corrNull = zeros(1,nsim_null);
rdcNull  = zeros(1,nsim_null);
cosNull = zeros(1,nsim_null);
ccorrNull = zeros(1,nsim_null);

rsdmAlt  = zeros(1,nsim_alt);
dcorrAlt = zeros(1,nsim_alt);
ticeAlt = zeros(1,nsim_alt);
corrAlt = zeros(1,nsim_alt);
rdcAlt  = zeros(1,nsim_alt);
cosAlt = zeros(1,nsim_alt);
ccorrAlt = zeros(1,nsim_alt);

% Arrays holding the estimated power for each of the "correlation" types, 
% for each data type (linear, parabolic, etc...) with each noise level
rsdmPower = zeros(numDepTests, num_noise, length(M_vec));
dcorrPower = zeros(numDepTests, num_noise, length(M_vec));
ticePower = zeros(numDepTests, num_noise, length(M_vec));
corrPower = zeros(numDepTests, num_noise, length(M_vec));
rdcPower  = zeros(numDepTests, num_noise, length(M_vec));
cosPower  = zeros(numDepTests, num_noise, length(M_vec));
ccorrPower  = zeros(numDepTests, num_noise, length(M_vec));

% Optimal parameters for MICe
mine_c = 15;
mine_alpha = 0.6;

% Optimal parameters for RDC
rdc_k = 20;
rdc_s = 1/6;

powerThreshold = 0.7;       % we want to see what sample-size is needed by
                            % each dependency metric to acheive 70% 
                            % statistical power

% We loop through the noise level and functional form; 
% each time we estimate a null distribution based on the marginals of the data, 
% and then use that null distribution to estimate power

% We use a uniformly distributed x, because in the original paper the 
% authors used the same

% Simon & Tibshirani use xMin=0, xMax=1 for performing their analysis ...
xMin = 0;
xMax = 1;

dispstat('','init'); % One time only initialization
dispstat(sprintf('Begining the simulation...\n'),'keepthis','timestamp');
num_noise_test_min = 1;
num_noise_test_max = 30;
for m=1:length(M_vec)
    M = M_vec(m);
    for l=num_noise_test_min:num_noise_test_max
        for typ=1:numDepTests
            dispstat(sprintf('M=%d Noise=%d Dependency=%d',M, l, typ),'keepthis', 'timestamp');
            % simulate data under the null w/ correct marginals
            parfor ii=1:nsim_null
                x = rand(M,1)*(xMax-xMin)+xMin;
                switch(typ)
                    case 1
                        % linear
                        y = x + noise*(l/num_noise)*randn(M,1); 
                    case 2
                        % parabolic
                        y = 4*(x-.5).^2 + noise*(l/num_noise)*randn(M,1);
                    case 3
                        % cubic
                        y = 128*(x-1/3).^3-48*(x-1/3).^3-12*(x-1/3)+10* noise*(l/num_noise)*randn(M,1);
                    case 4
                        % low-freq sin
                        y = sin(4*pi*x) + 2*noise*(l/num_noise)*randn(M,1);
                    case 5
                        % high-freq sin
                        y = sin(16*pi*x) + noise*(l/num_noise)*randn(M,1);
                    case 6
                        % fourth root
                        y = x.^(1/4) + noise*(l/num_noise)*randn(M,1);
                    case 7
                        % circle
                        y=(2*binornd(1,0.5,M,1)-1) .* (sqrt(1 - (2*x - 1).^2)) + noise/4*l/num_noise*randn(M,1);
                    case 8
                        % step function
                        y = (x > 0.5) + noise*5*l/num_noise*randn(M,1);
                    otherwise
                        error('unknown dep type!');
                end
                % resimulate x so we have null scenario
                x = rand(M,1)*(xMax-xMin)+xMin;

                % calculate the metrics
                rsdmNull(ii) = rsdm(x, y);
                dcorrNull(ii) = dcorr(x, y);
                % compute MICe
                minestats = mine(x',y',mine_alpha,mine_c,'mic_e');
                ticeNull(ii) = minestats.tic;
                % compute correlation
                corrNull(ii) = corr(x,y);
                % compute RDC
                rdcNull(ii) = rdc(x,y,rdc_k,rdc_s);
                % compute CoS
                cosNull(ii) = cosf(x,y);
                % compute ccorr
                ccorrNull(ii) = cCorr(x,y);
            end

            % compute the rejection cutoffs
            rsdm_cut = quantile(rsdmNull, 0.95);
            dcorr_cut = quantile(dcorrNull, 0.95);
            tice_cut = quantile(ticeNull, 0.95);
            corr_cut = quantile(corrNull, 0.95);
            rdc_cut  = quantile(rdcNull, 0.95);
            cos_cut  = quantile(cosNull, 0.95);
            ccorr_cut = quantile(ccorrNull, 0.95);

            % resimulate the data under the alternative hypothesis
            parfor ii=1:nsim_alt
                x = rand(M,1)*(xMax-xMin)+xMin;
                switch(typ)
                    case 1
                        % linear
                        y = x + noise*(l/num_noise)*randn(M,1); 
                    case 2
                        % parabolic
                        y = 4*(x-.5).^2 + noise*(l/num_noise)*randn(M,1);
                    case 3
                        % cubic
                        y = 128*(x-1/3).^3-48*(x-1/3).^3-12*(x-1/3) + 10*noise*(l/num_noise)*randn(M,1);
                    case 4
                        % low-freq sin
                        y = sin(4*pi*x) + 2*noise*(l/num_noise)*randn(M,1);
                    case 5
                        % high-freq sin
                        y = sin(16*pi*x) + noise*(l/num_noise)*randn(M,1);
                    case 6
                        % fourth root
                        y = x.^(1/4) + noise*(l/num_noise)*randn(M,1);
                    case 7
                        % circle
                        y=(2*binornd(1,0.5,M,1)-1) .* (sqrt(1 - (2*x - 1).^2)) + noise/4*l/num_noise*randn(M,1);
                    case 8
                        % step function
                        y = (x > 0.5) + noise*5*l/num_noise*randn(M,1);
                    otherwise
                        error('unknown dep type!');
                end

                % calculate the metrics
                rsdmAlt(ii) = rsdm(x, y);
                dcorrAlt(ii) = dcorr(x, y);
                % compute MICe
                minestats = mine(x',y',mine_alpha,mine_c,'mic_e');
                ticeAlt(ii) = minestats.tic;
                % compute corr
                corrAlt(ii) = corr(x, y);
                % compute RDC
                rdcAlt(ii) = rdc(x,y,rdc_k,rdc_s);
                % compute CoS
                cosAlt(ii) = cosf(x,y);
                % compute ccorr
                ccorrAlt(ii) = cCorr(x,y);
            end

            % compute the power
            rsdmPower(typ, l, m)   = sum(rsdmAlt > rsdm_cut)/nsim_alt;
            dcorrPower(typ, l, m)  = sum(dcorrAlt > dcorr_cut)/nsim_alt;
            ticePower(typ, l, m)   = sum(ticeAlt > tice_cut)/nsim_alt;
            corrPower(typ, l, m)   = sum(corrAlt > corr_cut)/nsim_alt;
            rdcPower(typ, l, m)    = sum(rdcAlt > rdc_cut)/nsim_alt;
            cosPower(typ, l, m)    = sum(cosAlt > cos_cut)/nsim_alt;
            ccorrPower(typ, l, m)  = sum(ccorrAlt > ccorr_cut)/nsim_alt;
        end
    end
    
    % save intermediate results just in case things crash :(
    if(ispc)
        save('C:\\Users\\Kiran\\ownCloud\\PhD\\sim_results\\independence\\power_all.mat');
    elseif(ismac)
        save('/Users/Kiran/ownCloud/PhD/sim_results/independence/power_all.mat');
    else
        save('/home/kiran/ownCloud/PhD/sim_results/independence/power_all.mat');
    end
end
% save the data
if(ispc)
    save('C:\\Users\\Kiran\\ownCloud\\PhD\\sim_results\\independence\\power_all.mat');
elseif(ismac)
    save('/Users/Kiran/ownCloud/PhD/sim_results/independence/power_all.mat');
else
    save('/home/kiran/ownCloud/PhD/sim_results/independence/power_all.mat');
end

% process these results to search for the minimum M that gives us the
% threshold power
depMetrics = {'rsdm', 'dcorr', 'tice', 'corr', 'rdc', 'cos', 'ccorr'};
sampleSizeAnalysisVec = zeros(numDepTests, length(depMetrics), length(num_noise_test_min:num_noise_test_max));
for depMetricIdx=1:length(depMetrics)
    if(strcmpi(depMetrics{depMetricIdx}, 'rsdm'))
        powerData = rsdmPower;
        rsdmIdx = 1;
    elseif(strcmpi(depMetrics{depMetricIdx}, 'dcorr'))
        powerData = dcorrPower;
        dcorrIdx = 2;
    elseif(strcmpi(depMetrics{depMetricIdx}, 'tice'))
        powerData = ticePower;
        ticeIdx = 3;
    elseif(strcmpi(depMetrics{depMetricIdx}, 'corr'))
        powerData = corrPower;
        corrIdx = 4;
    elseif(strcmpi(depMetrics{depMetricIdx}, 'rdc'))
        powerData = rdcPower;
        rdcIdx = 5;
    elseif(strcmpi(depMetrics{depMetricIdx}, 'cos'))
        powerData = cosPower;
        cosIdx = 6;
    elseif(strcmpi(depMetrics{depMetricIdx}, 'ccorr'))
        powerData = ccorrPower;
        ccorrIdx = 7;
    end
    
    for typ=1:numDepTests
        for l=num_noise_test_min:num_noise_test_max
            for m=1:length(M_vec)
                M = M_vec(m);
                if(powerData(typ, l, m)>powerThreshold)
                    % TODO: we can do some interpolation here, so that we
                    % are not restricted to the boundaries of which tests
                    % were run ...
                    sampleSizeAnalysisVec(typ, depMetricIdx, l) = M;
                    break;
                end
            end
        end
    end
end

% replace 0's w/ NaN's so we don't plot them
sampleSizeAnalysisVec(sampleSizeAnalysisVec==0)=NaN;

% inlet plot configuration
M_inlet = 200;
inset_bufX = 0.0005; inset_bufY = 0.26;
inset_width = 0.1; inset_height = 0.08;

noiseVec = (num_noise_test_min:num_noise_test_max)/10;
figure;
h1 = subplot(2,2,1);
hh1 = plot(noiseVec, squeeze(sampleSizeAnalysisVec(1,rsdmIdx,num_noise_test_min:num_noise_test_max)), 'o-.', ...
     noiseVec, squeeze(sampleSizeAnalysisVec(1,dcorrIdx,num_noise_test_min:num_noise_test_max)), '+-.', ...
     noiseVec, squeeze(sampleSizeAnalysisVec(1,ticeIdx,num_noise_test_min:num_noise_test_max)), 'd-.', ...
     noiseVec, squeeze(sampleSizeAnalysisVec(1,rdcIdx,num_noise_test_min:num_noise_test_max)), 'v-.', ...
     noiseVec, squeeze(sampleSizeAnalysisVec(1,cosIdx,num_noise_test_min:num_noise_test_max)), 's-.', ...
     noiseVec, squeeze(sampleSizeAnalysisVec(1,ccorrIdx,num_noise_test_min:num_noise_test_max)), 'p-.'); 
% axis([min(noiseVec) max(noiseVec) 0 1]);
xlabel('Noise Level', 'FontSize', 20); ylabel('$$M_{min}$$', 'FontSize', 20, 'Interpreter', 'Latex'); grid on;
h1.FontSize = 20; 
loc_inset = [h1.Position(1)+inset_bufX h1.Position(2)+inset_bufY inset_width inset_height];
ax1 = axes('Position',loc_inset);
tmp1 = linspace(0,1,M_inlet);
tmp2 = tmp1;
plot(tmp1,tmp2, 'k', 'LineWidth', 2);
ax1.Box = 'on'; ax1.XTick = []; ax1.YTick = [];
ax1.XLim = [min(tmp1) max(tmp1)];
ax1.YLim = [min(tmp2) max(tmp2)];
hh1(1).LineWidth = 5; 
hh1(2).LineWidth = 1.5; 
hh1(3).LineWidth = 1.5; 
hh1(4).LineWidth = 1.5; 
hh1(5).LineWidth = 1.5;
hh1(6).LineWidth = 1.5;

h2 = subplot(2,2,2);
hh2 = plot(noiseVec, squeeze(sampleSizeAnalysisVec(2,rsdmIdx,num_noise_test_min:num_noise_test_max)), 'o-.', ...
     noiseVec, squeeze(sampleSizeAnalysisVec(2,dcorrIdx,num_noise_test_min:num_noise_test_max)), '+-.', ...
     noiseVec, squeeze(sampleSizeAnalysisVec(2,ticeIdx,num_noise_test_min:num_noise_test_max)), 'd-.', ...
     noiseVec, squeeze(sampleSizeAnalysisVec(2,rdcIdx,num_noise_test_min:num_noise_test_max)), 'v-.', ...
     noiseVec, squeeze(sampleSizeAnalysisVec(2,cosIdx,num_noise_test_min:num_noise_test_max)), 's-.', ...
     noiseVec, squeeze(sampleSizeAnalysisVec(2,ccorrIdx,num_noise_test_min:num_noise_test_max)), 'p-.'); 
% axis([min(noiseVec) max(noiseVec) 0 1]);
xlabel('Noise Level', 'FontSize', 20); ylabel('$$M_{min}$$', 'FontSize', 20, 'Interpreter', 'Latex'); grid on;
h2.FontSize = 20; 
loc_inset = [h2.Position(1)+inset_bufX h2.Position(2)+inset_bufY inset_width inset_height];
ax1 = axes('Position',loc_inset);
tmp1 = linspace(0,1,M_inlet);
tmp2 = 4*(tmp1-.5).^2;
plot(tmp1,tmp2, 'k', 'LineWidth', 2);
ax1.Box = 'on'; ax1.XTick = []; ax1.YTick = [];
ax1.XLim = [min(tmp1) max(tmp1)];
ax1.YLim = [min(tmp2) max(tmp2)];
hh2(1).LineWidth = 5; 
hh2(2).LineWidth = 1.5; 
hh2(3).LineWidth = 1.5; 
hh2(4).LineWidth = 1.5; 
hh2(5).LineWidth = 1.5;
hh2(6).LineWidth = 1.5;

h3 = subplot(2,2,3);
hh3 = plot(noiseVec, squeeze(sampleSizeAnalysisVec(3,rsdmIdx,num_noise_test_min:num_noise_test_max)), 'o-.', ...
     noiseVec, squeeze(sampleSizeAnalysisVec(3,dcorrIdx,num_noise_test_min:num_noise_test_max)), '+-.', ...
     noiseVec, squeeze(sampleSizeAnalysisVec(3,ticeIdx,num_noise_test_min:num_noise_test_max)), 'd-.', ...
     noiseVec, squeeze(sampleSizeAnalysisVec(3,rdcIdx,num_noise_test_min:num_noise_test_max)), 'v-.', ...
     noiseVec, squeeze(sampleSizeAnalysisVec(3,cosIdx,num_noise_test_min:num_noise_test_max)), 's-.', ...
     noiseVec, squeeze(sampleSizeAnalysisVec(3,ccorrIdx,num_noise_test_min:num_noise_test_max)), 'p-.'); 
% axis([min(noiseVec) max(noiseVec) 0 1]);
xlabel('Noise Level', 'FontSize', 20); ylabel('$$M_{min}$$', 'FontSize', 20, 'Interpreter', 'Latex'); grid on;
h3.FontSize = 20; 
loc_inset = [h3.Position(1)+inset_bufX h3.Position(2)+inset_bufY inset_width inset_height];
ax1 = axes('Position',loc_inset);
tmp1 = linspace(0,1,M_inlet);
tmp2 = 128*(tmp1-1/3).^3-48*(tmp1-1/3).^3-12*(tmp1-1/3);
plot(tmp1,tmp2, 'k', 'LineWidth', 2);
ax1.Box = 'on'; ax1.XTick = []; ax1.YTick = [];
ax1.XLim = [min(tmp1) max(tmp1)];
ax1.YLim = [min(tmp2) max(tmp2)];
hh3(1).LineWidth = 5; 
hh3(2).LineWidth = 1.5; 
hh3(3).LineWidth = 1.5; 
hh3(4).LineWidth = 1.5; 
hh3(5).LineWidth = 1.5;
hh3(6).LineWidth = 1.5;

h4 = subplot(2,2,4);
hh4 = plot(noiseVec, squeeze(sampleSizeAnalysisVec(4,rsdmIdx,num_noise_test_min:num_noise_test_max)), 'o-.', ...
     noiseVec, squeeze(sampleSizeAnalysisVec(4,dcorrIdx,num_noise_test_min:num_noise_test_max)), '+-.', ...
     noiseVec, squeeze(sampleSizeAnalysisVec(4,ticeIdx,num_noise_test_min:num_noise_test_max)), 'd-.', ...
     noiseVec, squeeze(sampleSizeAnalysisVec(4,rdcIdx,num_noise_test_min:num_noise_test_max)), 'v-.', ...
     noiseVec, squeeze(sampleSizeAnalysisVec(4,cosIdx,num_noise_test_min:num_noise_test_max)), 's-.', ...
     noiseVec, squeeze(sampleSizeAnalysisVec(4,ccorrIdx,num_noise_test_min:num_noise_test_max)), 'p-.'); 
% axis([min(noiseVec) max(noiseVec) 0 1]);
xlabel('Noise Level', 'FontSize', 20); ylabel('$$M_{min}$$', 'FontSize', 20, 'Interpreter', 'Latex'); grid on;
h4.FontSize = 20; 
loc_inset = [h4.Position(1)+inset_bufX h4.Position(2)+inset_bufY inset_width inset_height];
ax1 = axes('Position',loc_inset);
tmp1 = linspace(0,1,M_inlet);
tmp2 = sin(4*pi*tmp1);
plot(tmp1,tmp2, 'k', 'LineWidth', 2);
ax1.Box = 'on'; ax1.XTick = []; ax1.YTick = [];
ax1.XLim = [min(tmp1) max(tmp1)];
ax1.YLim = [min(tmp2) max(tmp2)];
hh4(1).LineWidth = 5; 
hh4(2).LineWidth = 1.5; 
hh4(3).LineWidth = 1.5; 
hh4(4).LineWidth = 1.5; 
hh4(5).LineWidth = 1.5;
hh4(6).LineWidth = 1.5;

figure;
h5 = subplot(2,2,1);
hh5 = plot(noiseVec, squeeze(sampleSizeAnalysisVec(5,rsdmIdx,num_noise_test_min:num_noise_test_max)), 'o-.', ...
     noiseVec, squeeze(sampleSizeAnalysisVec(5,dcorrIdx,num_noise_test_min:num_noise_test_max)), '+-.', ...
     noiseVec, squeeze(sampleSizeAnalysisVec(5,ticeIdx,num_noise_test_min:num_noise_test_max)), 'd-.', ...
     noiseVec, squeeze(sampleSizeAnalysisVec(5,rdcIdx,num_noise_test_min:num_noise_test_max)), 'v-.', ...
     noiseVec, squeeze(sampleSizeAnalysisVec(5,cosIdx,num_noise_test_min:num_noise_test_max)), 's-.', ...
     noiseVec, squeeze(sampleSizeAnalysisVec(5,ccorrIdx,num_noise_test_min:num_noise_test_max)), 'p-.'); 
% axis([min(noiseVec) max(noiseVec) 0 1]);
xlabel('Noise Level', 'FontSize', 20); ylabel('$$M_{min}$$', 'FontSize', 20, 'Interpreter', 'Latex'); grid on;
h5.FontSize = 20; 
loc_inset = [h5.Position(1)+inset_bufX h5.Position(2)+inset_bufY inset_width inset_height];
ax1 = axes('Position',loc_inset);
tmp1 = linspace(0,1,M_inlet);
tmp2 = sin(16*pi*tmp1);
plot(tmp1,tmp2, 'k', 'LineWidth', 2);
ax1.Box = 'on'; ax1.XTick = []; ax1.YTick = [];
ax1.XLim = [min(tmp1) max(tmp1)];
ax1.YLim = [min(tmp2) max(tmp2)];
hh5(1).LineWidth = 5; 
hh5(2).LineWidth = 1.5; 
hh5(3).LineWidth = 1.5; 
hh5(4).LineWidth = 1.5; 
hh5(5).LineWidth = 1.5;
hh5(6).LineWidth = 1.5;

h6 = subplot(2,2,2);
hh6 = plot(noiseVec, squeeze(sampleSizeAnalysisVec(6,rsdmIdx,num_noise_test_min:num_noise_test_max)), 'o-.', ...
     noiseVec, squeeze(sampleSizeAnalysisVec(6,dcorrIdx,num_noise_test_min:num_noise_test_max)), '+-.', ...
     noiseVec, squeeze(sampleSizeAnalysisVec(6,ticeIdx,num_noise_test_min:num_noise_test_max)), 'd-.', ...
     noiseVec, squeeze(sampleSizeAnalysisVec(6,rdcIdx,num_noise_test_min:num_noise_test_max)), 'v-.', ...
     noiseVec, squeeze(sampleSizeAnalysisVec(6,cosIdx,num_noise_test_min:num_noise_test_max)), 's-.', ...
     noiseVec, squeeze(sampleSizeAnalysisVec(6,ccorrIdx,num_noise_test_min:num_noise_test_max)), 'p-.'); 
% axis([min(noiseVec) max(noiseVec) 0 1]);
xlabel('Noise Level', 'FontSize', 20); ylabel('$$M_{min}$$', 'FontSize', 20, 'Interpreter', 'Latex'); grid on;
h6.FontSize = 20; 
loc_inset = [h6.Position(1)+inset_bufX h6.Position(2)+inset_bufY inset_width inset_height];
ax1 = axes('Position',loc_inset);
tmp1 = linspace(0,1,M_inlet);
tmp2 = tmp1.^(1/4);
plot(tmp1,tmp2, 'k', 'LineWidth', 2);
ax1.Box = 'on'; ax1.XTick = []; ax1.YTick = [];
ax1.XLim = [min(tmp1) max(tmp1)];
ax1.YLim = [min(tmp2) max(tmp2)];
hh6(1).LineWidth = 5; 
hh6(2).LineWidth = 1.5; 
hh6(3).LineWidth = 1.5; 
hh6(4).LineWidth = 1.5; 
hh6(5).LineWidth = 1.5;
hh6(6).LineWidth = 1.5;

h7 = subplot(2,2,3);
hh7 = plot(noiseVec, squeeze(sampleSizeAnalysisVec(7,rsdmIdx,num_noise_test_min:num_noise_test_max)), 'o-.', ...
     noiseVec, squeeze(sampleSizeAnalysisVec(7,dcorrIdx,num_noise_test_min:num_noise_test_max)), '+-.', ...
     noiseVec, squeeze(sampleSizeAnalysisVec(7,ticeIdx,num_noise_test_min:num_noise_test_max)), 'd-.', ...
     noiseVec, squeeze(sampleSizeAnalysisVec(7,rdcIdx,num_noise_test_min:num_noise_test_max)), 'v-.', ...
     noiseVec, squeeze(sampleSizeAnalysisVec(7,cosIdx,num_noise_test_min:num_noise_test_max)), 's-.', ...
     noiseVec, squeeze(sampleSizeAnalysisVec(7,ccorrIdx,num_noise_test_min:num_noise_test_max)), 'p-.'); 
% axis([min(noiseVec) max(noiseVec) 0 1]);
xlabel('Noise Level', 'FontSize', 20); ylabel('$$M_{min}$$', 'FontSize', 20, 'Interpreter', 'Latex'); grid on;
h7.FontSize = 20; 
loc_inset = [h7.Position(1)+inset_bufX h7.Position(2)+inset_bufY inset_width inset_height];
ax1 = axes('Position',loc_inset);
tmp1 = linspace(0,1,M_inlet/2);
tmp2 = (sqrt(1 - (2*tmp1 - 1).^2));
tmp3 = -(sqrt(1 - (2*tmp1 - 1).^2));
plot(tmp1,tmp2, 'k', 'LineWidth', 2); hold on;
plot(tmp1,tmp3, 'k', 'LineWidth', 2); 
ax1.Box = 'on'; ax1.XTick = []; ax1.YTick = [];
ax1.XLim = [min(tmp1) max(tmp1)];
ax1.YLim = [min(tmp3) max(tmp2)];
hh7(1).LineWidth = 5; 
hh7(2).LineWidth = 1.5; 
hh7(3).LineWidth = 1.5; 
hh7(4).LineWidth = 1.5; 
hh7(5).LineWidth = 1.5;
hh7(6).LineWidth = 1.5;

h8 = subplot(2,2,4);
hh8 = plot(noiseVec, squeeze(sampleSizeAnalysisVec(8,rsdmIdx,num_noise_test_min:num_noise_test_max)), 'o-.', ...
     noiseVec, squeeze(sampleSizeAnalysisVec(8,dcorrIdx,num_noise_test_min:num_noise_test_max)), '+-.', ...
     noiseVec, squeeze(sampleSizeAnalysisVec(8,ticeIdx,num_noise_test_min:num_noise_test_max)), 'd-.', ...
     noiseVec, squeeze(sampleSizeAnalysisVec(8,rdcIdx,num_noise_test_min:num_noise_test_max)), 'v-.', ...
     noiseVec, squeeze(sampleSizeAnalysisVec(8,cosIdx,num_noise_test_min:num_noise_test_max)), 's-.', ...
     noiseVec, squeeze(sampleSizeAnalysisVec(8,ccorrIdx,num_noise_test_min:num_noise_test_max)), 'p-.'); 
legend('RSDM', 'dCor', 'TIC_e', 'RDC', 'CoS', 'cCor');  % manually move this using the mouse to a
                                                  % good location
% axis([min(noiseVec) max(noiseVec) 0 1]);
xlabel('Noise Level', 'FontSize', 20); ylabel('$$M_{min}$$', 'FontSize', 20, 'Interpreter', 'Latex'); grid on;
h8.FontSize = 20; 
loc_inset = [h8.Position(1)+inset_bufX h8.Position(2)+inset_bufY inset_width inset_height];
ax1 = axes('Position',loc_inset);
tmp1 = linspace(0,1,M_inlet);
tmp2 = (tmp1 > 0.5);
plot(tmp1,tmp2, 'k', 'LineWidth', 2);
ax1.Box = 'on'; ax1.XTick = []; ax1.YTick = [];
ax1.XLim = [min(tmp1) max(tmp1)];
% ax1.YLim = [min(tmp2) max(tmp2)]; % why is this erroring?
ax1.YLim = [0 1];
hh8(1).LineWidth = 5; 
hh8(2).LineWidth = 1.5; 
hh8(3).LineWidth = 1.5; 
hh8(4).LineWidth = 1.5; 
hh8(5).LineWidth = 1.5;
hh8(6).LineWidth = 1.5;

%% Characterize null distribution (X indep Y) experimentally for \hat{RSDM}
clear;
clc;

rng(1234);

nsim = 1000;
M_vec = 100:100:1000;

xMin = 0; xMax = 1;
yMin = 0; yMax = 1;

FIT_PLOTS = 0;

rsdmNullDistributionResults = zeros(nsim, length(M_vec));

for ii=1:nsim
    parfor jj=1:length(M_vec)
        M = M_vec(jj);
        % create independent x & y
        x = rand(M,1)*(xMax-xMin)+xMin;
        y = rand(M,1)*(yMax-yMin)+yMin;
    
        % compute RSDM
        rsdmNullDistributionResults(ii,jj) = rsdm(x, y);
    end
end

% plot distribution of RSDM under the null distribution 
legendCell = cell(1,length(M_vec));
for ii=1:length(M_vec)
    [f,xi] = ksdensity(rsdmNullDistributionResults(:,ii));
    plot(xi,f); hold on;
    legendCell{ii} = sprintf('M=%d',M_vec(ii));
end
grid on;
legend(legendCell);
title('Distribution of RSDM_{approx}');

D_cell = cell(1,length(M_vec)); 
PD_cell = cell(1,length(M_vec));
idx = 1;
for ii=1:length(M_vec)
    [D, PD] = allfitdist(rsdmNullDistributionResults(:,ii), 'PDF');
    D_cell{idx} = D; 
    PD_cell{idx} = PD;
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

distScores = zeros(4,length(distributions));
for ii=1:length(distributions)
    dist = distributions{ii};
    % find this distribution in the fit and store the BIC, AIC, AICc scores
    % for all M
    NLogL = 0;
    BIC = 0;
    AIC = 0;
    AICc = 0;
    for jj=1:length(M_vec)
        D = D_cell{jj};
        PD = PD_cell{jj};
        
        % find the distribution
        for kk=1:length(PD)
            if(strcmpi(PD{kk}.DistributionName, dist))
                break;
            end
        end
        
        NLogL = NLogL + D(kk).NLogL;
        BIC = BIC + D(kk).BIC;
        AIC = AIC + D(kk).AIC;
        AICc = AICc + D(kk).AICc;
    end
    
    distScores(1,ii) = NLogL;
    distScores(2,ii) = BIC;
    distScores(3,ii) = AIC;
    distScores(4,ii) = AICc;
end

% save the data
if(ispc)
    save('C:\\Users\\Kiran\\ownCloud\\PhD\\sim_results\\independence\\rsdmNullDistribution.mat');
elseif(ismac)
    save('/Users/Kiran/ownCloud/PhD/sim_results/independence/rsdmNullDistribution.mat');
else
    save('/home/kiran/ownCloud/PhD/sim_results/independence/rsdmNullDistribution.mat');
end

% Sort by NLogL
[~,I] = sort(distScores(1,:), 'ascend');
fprintf('NLogL\n');
distributions{I(1)}

% Sort by BIC
[~,I] = sort(distScores(2,:), 'ascend');
fprintf('BIC\n');
distributions{I(1)}

% Sort by AIC
[~,I] = sort(distScores(3,:), 'ascend');
fprintf('AIC\n');
distributions{I(1)}

% Sort by AICc
[~,I] = sort(distScores(4,:), 'ascend');
fprintf('AICc\n');
distributions{I(1)}

%%

if(ispc)
    load('C:\\Users\\Kiran\\ownCloud\\PhD\\sim_results\\independence\\rsdmNullDistribution.mat');
elseif(ismac)
    load('/Users/Kiran/ownCloud/PhD/sim_results/independence/rsdmNullDistribution.mat');
else
    load('/home/kiran/ownCloud/PhD/sim_results/independence/rsdmNullDistribution.mat');
end


% From the above analysis, the Generalized Extreme Value distribution seems 
% to fit best ... Q-Q Plots

% QQ Plot w/ best fit for M=100 and M=1000
pdObjs = cell(1,length(M_vec));
kVec = zeros(1,length(M_vec));
muVec = zeros(1,length(M_vec));
sigmaVec = zeros(1,length(M_vec));
for ii=1:length(M_vec)
    M = M_vec(ii);
    % look for the Inverse Gaussian Distribution in the correct cell array
    D = D_cell(ii);
    PD = PD_cell(ii); PD = PD{1};
    for jj=1:length(PD)
        if(strcmpi('Generalized extreme value', PD{jj}.DistributionName))
            pd = PD{jj};
        end
    end
    pdObjs{ii} = pd;
    kVec(ii) = pd.k;
    muVec(ii) = pd.mu;
    sigmaVec(ii) = pd.sigma;
end

fontSize = 20;

% do the Q-Q plot
pd = pdObjs{5};
h1 = subplot(2,2,1); qqplot(rsdmNullDistributionResults(:,1), pd); grid on;
xlabel(sprintf('Quantiles of GEV(%0.02f,%0.02f, %0.02f)', pd.k, pd.mu, pd.sigma), 'FontSize', 20);
ylabel('Quantiles of Input Samples', 'FontSize', fontSize);
title('M = 500', 'FontSize', fontSize);
h1.FontSize = fontSize;

% plot how mu and lambda change as M goes from 100 --> 1000
h2 = subplot(2,2,2); plot(M_vec, kVec, 'LineWidth', 5);     
grid on; xlabel('M', 'FontSize', fontSize); 
ylabel('k', 'FontSize', fontSize);
h2.FontSize = fontSize;

h3 = subplot(2,2,3); plot(M_vec, muVec, 'LineWidth', 5);   
grid on; xlabel('M', 'FontSize', fontSize); 
ylabel('\mu', 'FontSize', fontSize);
h3.FontSize = fontSize;

h4 = subplot(2,2,4); plot(M_vec, sigmaVec, 'LineWidth', 5);
grid on; xlabel('M', 'FontSize', fontSize); 
ylabel('\lambda', 'FontSize', fontSize);
h4.FontSize = fontSize;

%% Generate curves which show the effect of sample size for statistical for only CoS and cCorr
% TODO: merge these results back in w/ the other power vs sample-size
% calculations ...

clear;
clc;

% WARNING: ENSURE THAT minepy/matlab/ is in the matlab path for MIC to
% work!

rng(1234);
dbstop if error;

nsim_null = 500;   % The number of null datasets we use to estimate our rejection reject regions for an alternative with level 0.05
nsim_alt  = 500;   % Number of alternative datasets we use to estimate our power

num_noise = 30;                    % The number of different noise levels used
noise = 3;                         % A constant to determine the amount of noise

M_vec = 25:25:1500;      % number of samples
numDepTests = 8;        % the number of different dependency tests we will conduct
                        % TODO: add copula dependencies as well

% Optimal parameters for MICe
mine_c = 15;
mine_alpha = 0.6;
                        
% Vectors holding the null "correlations" (for pearson, dcor and mic respectively) 
% for each of the nsim null datasets at a given noise level
cosNull = zeros(1,nsim_null);
ccorrNull = zeros(1,nsim_null);
ticeNull = zeros(1,nsim_null);

cosAlt = zeros(1,nsim_alt);
ccorrAlt = zeros(1,nsim_alt);
ticeAlt = zeros(1,nsim_alt);

% Arrays holding the estimated power for each of the "correlation" types, 
% for each data type (linear, parabolic, etc...) with each noise level
cosPower  = zeros(numDepTests, num_noise, length(M_vec));
ccorrPower  = zeros(numDepTests, num_noise, length(M_vec));
ticePower = zeros(numDepTests, num_noise, length(M_vec));

% We loop through the noise level and functional form; 
% each time we estimate a null distribution based on the marginals of the data, 
% and then use that null distribution to estimate power

% We use a uniformly distributed x, because in the original paper the 
% authors used the same

% Simon & Tibshirani use xMin=0, xMax=1 for performing their analysis ...
xMin = 0;
xMax = 1;

dispstat('','init'); % One time only initialization
dispstat(sprintf('Begining the simulation...\n'),'keepthis','timestamp');
num_noise_test_min = 1;
num_noise_test_max = 30;
% for m=1:length(M_vec)
for m=1:30
    M = M_vec(m);
    for l=num_noise_test_min:num_noise_test_max
        for typ=1:numDepTests
            dispstat(sprintf('M=%d Noise=%d Dependency=%d',M, l, typ),'keepthis', 'timestamp');
            % simulate data under the null w/ correct marginals
            parfor ii=1:nsim_null
                x = rand(M,1)*(xMax-xMin)+xMin;
                switch(typ)
                    case 1
                        % linear
                        y = x + noise*(l/num_noise)*randn(M,1); 
                    case 2
                        % parabolic
                        y = 4*(x-.5).^2 + noise*(l/num_noise)*randn(M,1);
                    case 3
                        % cubic
                        y = 128*(x-1/3).^3-48*(x-1/3).^3-12*(x-1/3)+10* noise*(l/num_noise)*randn(M,1);
                    case 4
                        % low-freq sin
                        y = sin(4*pi*x) + 2*noise*(l/num_noise)*randn(M,1);
                    case 5
                        % high-freq sin
                        y = sin(16*pi*x) + noise*(l/num_noise)*randn(M,1);
                    case 6
                        % fourth root
                        y = x.^(1/4) + noise*(l/num_noise)*randn(M,1);
                    case 7
                        % circle
                        y=(2*binornd(1,0.5,M,1)-1) .* (sqrt(1 - (2*x - 1).^2)) + noise/4*l/num_noise*randn(M,1);
                    case 8
                        % step function
                        y = (x > 0.5) + noise*5*l/num_noise*randn(M,1);
                    otherwise
                        error('unknown dep type!');
                end
                % resimulate x so we have null scenario
                x = rand(M,1)*(xMax-xMin)+xMin;

                % calculate the metrics
                % compute TICe
                minestats = mine(x',y',mine_alpha,mine_c,'mic_e');
                ticeNull(ii) = minestats.tic;
                % compute CoS
                cosNull(ii) = cosf(x,y);
                % compute ccorr
                ccorrNull(ii) = cCorr(x,y);
            end

            % compute the rejection cutoffs
            cos_cut  = quantile(cosNull, 0.95);
            ccorr_cut = quantile(ccorrNull, 0.95);
            tice_cut = quantile(ticeNull, 0.95);

            % resimulate the data under the alternative hypothesis
            parfor ii=1:nsim_alt
                x = rand(M,1)*(xMax-xMin)+xMin;
                switch(typ)
                    case 1
                        % linear
                        y = x + noise*(l/num_noise)*randn(M,1); 
                    case 2
                        % parabolic
                        y = 4*(x-.5).^2 + noise*(l/num_noise)*randn(M,1);
                    case 3
                        % cubic
                        y = 128*(x-1/3).^3-48*(x-1/3).^3-12*(x-1/3) + 10*noise*(l/num_noise)*randn(M,1);
                    case 4
                        % low-freq sin
                        y = sin(4*pi*x) + 2*noise*(l/num_noise)*randn(M,1);
                    case 5
                        % high-freq sin
                        y = sin(16*pi*x) + noise*(l/num_noise)*randn(M,1);
                    case 6
                        % fourth root
                        y = x.^(1/4) + noise*(l/num_noise)*randn(M,1);
                    case 7
                        % circle
                        y=(2*binornd(1,0.5,M,1)-1) .* (sqrt(1 - (2*x - 1).^2)) + noise/4*l/num_noise*randn(M,1);
                    case 8
                        % step function
                        y = (x > 0.5) + noise*5*l/num_noise*randn(M,1);
                    otherwise
                        error('unknown dep type!');
                end

                % calculate the metrics
                % compute TICe
                minestats = mine(x',y',mine_alpha,mine_c,'mic_e');
                ticeAlt(ii) = minestats.tic;
                % compute CoS
                cosAlt(ii) = cosf(x,y);
                % compute ccorr
                ccorrAlt(ii) = cCorr(x,y);
            end

            % compute the power
            cosPower(typ, l, m)    = sum(cosAlt > cos_cut)/nsim_alt;
            ccorrPower(typ, l, m)  = sum(ccorrAlt > ccorr_cut)/nsim_alt;
            ticePower(typ, l, m)   = sum(ticeAlt > tice_cut)/nsim_alt;
        end
    end
    
    % save intermediate results just in case things crash :(
    if(ispc)
        save('C:\\Users\\Kiran\\ownCloud\\PhD\\sim_results\\independence\\rsdmPower_CoS_cCorr_ticE_M_25_750.mat');
    elseif(ismac)
        save('/Users/Kiran/ownCloud/PhD/sim_results/independence/rsdmPower_CoS_cCorr_ticE_M_25_750.mat');
    else
        save('/home/kiran/ownCloud/PhD/sim_results/independence/rsdmPower_CoS_cCorr_ticE_M_25_750.mat');
    end
end
% save the data
if(ispc)
    save('C:\\Users\\Kiran\\ownCloud\\PhD\\sim_results\\independence\\rsdmPower_CoS_cCorr_ticE_M_25_750.mat');
elseif(ismac)
    save('/Users/Kiran/ownCloud/PhD/sim_results/independence/rsdmPower_CoS_cCorr_ticE_M_25_750.mat');
else
    save('/home/kiran/ownCloud/PhD/sim_results/independence/rsdmPower_CoS_cCorr_ticE_M_25_750.mat');
end