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

% try to get 8 workers :D
myCluster = parcluster('local');
myCluster.NumWorkers = 8;
saveProfile(myCluster);

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
rsdm1Null = zeros(1,nsim_null);
rsdm2Null = zeros(1,nsim_null);
rsdm4Null = zeros(1,nsim_null);
dcorrNull = zeros(1,nsim_null);
miceNull = zeros(1,nsim_null);
corrNull = zeros(1,nsim_null);
rdcNull  = zeros(1,nsim_null);

rsdm1Alt  = zeros(1,nsim_alt);
rsdm2Alt  = zeros(1,nsim_alt);
rsdm4Alt  = zeros(1,nsim_alt);
dcorrAlt = zeros(1,nsim_alt);
miceAlt = zeros(1,nsim_alt);
corrAlt = zeros(1,nsim_alt);
rdcAlt  = zeros(1,nsim_alt);

% Arrays holding the estimated power for each of the "correlation" types, 
% for each data type (linear, parabolic, etc...) with each noise level
rsdm1Power = zeros(numDepTests, num_noise);
rsdm2Power = zeros(numDepTests, num_noise);
rsdm4Power = zeros(numDepTests, num_noise);
dcorrPower = zeros(numDepTests, num_noise);
micePower = zeros(numDepTests, num_noise);
corrPower = zeros(numDepTests, num_noise);
rdcPower  = zeros(numDepTests, num_noise);

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
            rsdm1Null(ii) = rsdm_1(x, y);
            rsdm2Null(ii) = rsdm_2(x, y);
            rsdm4Null(ii) = rsdm_4(x, y);
            dcorrNull(ii) = dcorr(x, y);
            % compute MICe
            minestats = mine(x',y',mine_alpha,mine_c,'mic_e');
            miceNull(ii) = minestats.mic;
            % compute correlation
            corrNull(ii) = corr(x,y);
            % compute RDC
            rdcNull(ii) = rdc(x,y,rdc_k,rdc_s);
        end
        
        % compute the rejection cutoffs
        rsdm1_cut = quantile(rsdm1Null, 0.95);
        rsdm2_cut = quantile(rsdm2Null, 0.95);
        rsdm4_cut = quantile(rsdm4Null, 0.95);
        dcorr_cut = quantile(dcorrNull, 0.95);
        mice_cut = quantile(miceNull, 0.95);
        corr_cut = quantile(corrNull, 0.95);
        rdc_cut  = quantile(rdcNull, 0.95);
        
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
            rsdm1Alt(ii) = rsdm_1(x, y);
            rsdm2Alt(ii) = rsdm_2(x, y);
            rsdm4Alt(ii) = rsdm_4(x, y);
            dcorrAlt(ii) = dcorr(x, y);
            % compute MICe
            minestats = mine(x',y',mine_alpha,mine_c,'mic_e');
            miceAlt(ii) = minestats.mic;
            % compute corr
            corrAlt(ii) = corr(x, y);
            % compute RDC
            rdcAlt(ii) = rdc(x,y,rdc_k,rdc_s);
        end
        
        % compute the power
        rsdm1Power(typ, l)  = sum(rsdm1Alt > rsdm1_cut)/nsim_alt;
        rsdm2Power(typ, l)  = sum(rsdm2Alt > rsdm2_cut)/nsim_alt;
        rsdm4Power(typ, l)  = sum(rsdm4Alt > rsdm4_cut)/nsim_alt;
        dcorrPower(typ, l)  = sum(dcorrAlt > dcorr_cut)/nsim_alt;
        micePower(typ, l)   = sum(miceAlt > mice_cut)/nsim_alt;
        corrPower(typ, l)   = sum(corrAlt > corr_cut)/nsim_alt;
        rdcPower(typ, l)    = sum(rdcAlt > rdc_cut)/nsim_alt;
    end
end

% save the data
if(ispc)
    save(sprintf('C:\\Users\\Kiran\\ownCloud\\PhD\\sim_results\\independence\\rsdmParametricsPower_M_%d.mat', M));
elseif(ismac)
    save(sprintf('/Users/Kiran/ownCloud/PhD/sim_results/independence/rsdmParametricsPower_M_%d.mat', M));
else
    save(sprintf('/home/kiran/ownCloud/PhD/sim_results/independence/rsdmParametricsPower_M_%d.mat', M));
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
hh1 = plot(noiseVec, rsdm1Power(1,num_noise_test_min:num_noise_test_max), 'o-.', ...
     noiseVec, rsdm2Power(1,num_noise_test_min:num_noise_test_max), 's-.', ...
     noiseVec, rsdm4Power(1,num_noise_test_min:num_noise_test_max), 'h-.', ...
     noiseVec, dcorrPower(1,num_noise_test_min:num_noise_test_max), '+-.', ...
     noiseVec, micePower(1,num_noise_test_min:num_noise_test_max), 'd-.', ...
     noiseVec, corrPower(1,num_noise_test_min:num_noise_test_max), '^-.', ...
     noiseVec, rdcPower(1,num_noise_test_min:num_noise_test_max), 'v-.'); 
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

h2 = subplot(2,2,2);
hh2 = plot(noiseVec, rsdm1Power(2,num_noise_test_min:num_noise_test_max), 'o-.', ...
     noiseVec, rsdm2Power(2,num_noise_test_min:num_noise_test_max), 's-.', ...
     noiseVec, rsdm4Power(2,num_noise_test_min:num_noise_test_max), 'h-.', ...
     noiseVec, dcorrPower(2,num_noise_test_min:num_noise_test_max), '+-.', ...
     noiseVec, micePower(2,num_noise_test_min:num_noise_test_max), 'd-.', ...
     noiseVec, corrPower(2,num_noise_test_min:num_noise_test_max), '^-.', ...
     noiseVec, rdcPower(2,num_noise_test_min:num_noise_test_max), 'v-.'); 
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

h3 = subplot(2,2,3); 
hh3 = plot(noiseVec, rsdm1Power(3,num_noise_test_min:num_noise_test_max), 'o-.', ...
     noiseVec, rsdm2Power(3,num_noise_test_min:num_noise_test_max), 's-.', ...
     noiseVec, rsdm4Power(3,num_noise_test_min:num_noise_test_max), 'h-.', ...
     noiseVec, dcorrPower(3,num_noise_test_min:num_noise_test_max), '+-.', ...
     noiseVec, micePower(3,num_noise_test_min:num_noise_test_max), 'd-.', ...
     noiseVec, corrPower(3,num_noise_test_min:num_noise_test_max), '^-.', ...
     noiseVec, rdcPower(3,num_noise_test_min:num_noise_test_max), 'v-.');  
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

h4 = subplot(2,2,4); 
hh4 = plot(noiseVec, rsdm1Power(4,num_noise_test_min:num_noise_test_max), 'o-.', ...
     noiseVec, rsdm2Power(4,num_noise_test_min:num_noise_test_max), 's-.', ...
     noiseVec, rsdm4Power(4,num_noise_test_min:num_noise_test_max), 'h-.', ...
     noiseVec, dcorrPower(4,num_noise_test_min:num_noise_test_max), '+-.', ...
     noiseVec, micePower(4,num_noise_test_min:num_noise_test_max), 'd-.', ...
     noiseVec, corrPower(4,num_noise_test_min:num_noise_test_max), '^-.', ...
     noiseVec, rdcPower(4,num_noise_test_min:num_noise_test_max), 'v-.');  
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

figure;
h5 = subplot(2,2,1); 
hh5 = plot(noiseVec, rsdm1Power(5,num_noise_test_min:num_noise_test_max), 'o-.', ...
     noiseVec, rsdm2Power(5,num_noise_test_min:num_noise_test_max), 's-.', ...
     noiseVec, rsdm4Power(5,num_noise_test_min:num_noise_test_max), 'h-.', ...
     noiseVec, dcorrPower(5,num_noise_test_min:num_noise_test_max), '+-.', ...
     noiseVec, micePower(5,num_noise_test_min:num_noise_test_max), 'd-.', ...
     noiseVec, corrPower(5,num_noise_test_min:num_noise_test_max), '^-.', ...
     noiseVec, rdcPower(5,num_noise_test_min:num_noise_test_max), 'v-.');  
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

h6 = subplot(2,2,2); 
hh6 = plot(noiseVec, rsdm1Power(6,num_noise_test_min:num_noise_test_max), 'o-.', ...
     noiseVec, rsdm2Power(6,num_noise_test_min:num_noise_test_max), 's-.', ...
     noiseVec, rsdm4Power(6,num_noise_test_min:num_noise_test_max), 'h-.', ...
     noiseVec, dcorrPower(6,num_noise_test_min:num_noise_test_max), '+-.', ...
     noiseVec, micePower(6,num_noise_test_min:num_noise_test_max), 'd-.', ...
     noiseVec, corrPower(6,num_noise_test_min:num_noise_test_max), '^-.', ...
     noiseVec, rdcPower(6,num_noise_test_min:num_noise_test_max), 'v-.'); 
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

h7 = subplot(2,2,3); 
hh7 = plot(noiseVec, rsdm1Power(7,num_noise_test_min:num_noise_test_max), 'o-.', ...
     noiseVec, rsdm2Power(7,num_noise_test_min:num_noise_test_max), 's-.', ...
     noiseVec, rsdm4Power(7,num_noise_test_min:num_noise_test_max), 'h-.', ...
     noiseVec, dcorrPower(7,num_noise_test_min:num_noise_test_max), '+-.', ...
     noiseVec, micePower(7,num_noise_test_min:num_noise_test_max), 'd-.', ...
     noiseVec, corrPower(7,num_noise_test_min:num_noise_test_max), '^-.', ...
     noiseVec, rdcPower(7,num_noise_test_min:num_noise_test_max), 'v-.'); 
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

h8 = subplot(2,2,4); 
hh8 = plot(noiseVec, rsdm1Power(8,num_noise_test_min:num_noise_test_max), 'o-.', ...
     noiseVec, rsdm2Power(8,num_noise_test_min:num_noise_test_max), 's-.', ...
     noiseVec, rsdm4Power(8,num_noise_test_min:num_noise_test_max), 'h-.', ...
     noiseVec, dcorrPower(8,num_noise_test_min:num_noise_test_max), '+-.', ...
     noiseVec, micePower(8,num_noise_test_min:num_noise_test_max), 'd-.', ...
     noiseVec, corrPower(8,num_noise_test_min:num_noise_test_max), '^-.', ...
     noiseVec, rdcPower(8,num_noise_test_min:num_noise_test_max), 'v-.');  
axis([min(noiseVec) max(noiseVec) 0 1]);
h8.FontSize = 20; 
legend('RSDM_1', 'RSDM_2', 'dcorr', 'MIC_e', 'corr', 'RDC');  % manually move this using the mouse to a
                                                  % good location
xlabel('Noise Level'); ylabel('Power'); grid on;
loc_inset = [h8.Position(1)+inset_bufX h8.Position(2)+inset_bufY inset_width inset_height];
ax8 = axes('Position',loc_inset);
tmp1 = linspace(0,1,M_inlet);
tmp2 = (tmp1 > 0.5);
plot(tmp1,tmp2, 'k', 'LineWidth', 2);
ax8.Box = 'on'; ax8.XTick = []; ax8.YTick = [];
ax8.XLim = [min(tmp1) max(tmp1)];