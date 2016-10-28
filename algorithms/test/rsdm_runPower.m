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
miceNull = zeros(1,nsim_null);
corrNull = zeros(1,nsim_null);
rdcNull  = zeros(1,nsim_null);

rsdmAlt  = zeros(1,nsim_alt);
dcorrAlt = zeros(1,nsim_alt);
miceAlt = zeros(1,nsim_alt);
corrAlt = zeros(1,nsim_alt);
rdcAlt  = zeros(1,nsim_alt);

% Arrays holding the estimated power for each of the "correlation" types, 
% for each data type (linear, parabolic, etc...) with each noise level
rsdmPower = zeros(numDepTests, num_noise);
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
            rsdmNull(ii) = rsdm(x, y);
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
        rsdm_cut = quantile(rsdmNull, 0.95);
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
            rsdmAlt(ii) = rsdm(x, y);
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
        rsdmPower(typ, l)  = sum(rsdmAlt > rsdm_cut)/nsim_alt;
        dcorrPower(typ, l)  = sum(dcorrAlt > dcorr_cut)/nsim_alt;
        micePower(typ, l)   = sum(miceAlt > mice_cut)/nsim_alt;
        corrPower(typ, l)   = sum(corrAlt > corr_cut)/nsim_alt;
        rdcPower(typ, l)    = sum(rdcAlt > rdc_cut)/nsim_alt;
    end
end

% save the data
if(ispc)
    save(sprintf('C:\\Users\\Kiran\\ownCloud\\PhD\\sim_results\\independence\\rsdmPower_M_%d.mat', M));
elseif(ismac)
    save(sprintf('/Users/Kiran/ownCloud/PhD/sim_results/independence/rsdmPower_M_%d.mat', M));
else
    save(sprintf('/home/kiran/ownCloud/PhD/sim_results/independence/rsdmPower_M_%d.mat', M));
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
hh2 = plot(noiseVec, rsdmPower(2,num_noise_test_min:num_noise_test_max), 'o-.', ...
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
hh3 = plot(noiseVec, rsdmPower(3,num_noise_test_min:num_noise_test_max), 'o-.', ...
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
hh4 = plot(noiseVec, rsdmPower(4,num_noise_test_min:num_noise_test_max), 'o-.', ...
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
hh5 = plot(noiseVec, rsdmPower(5,num_noise_test_min:num_noise_test_max), 'o-.', ...
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
hh6 = plot(noiseVec, rsdmPower(6,num_noise_test_min:num_noise_test_max), 'o-.', ...
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
hh7 = plot(noiseVec, rsdmPower(7,num_noise_test_min:num_noise_test_max), 'o-.', ...
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
hh8 = plot(noiseVec, rsdmPower(8,num_noise_test_min:num_noise_test_max), 'o-.', ...
     noiseVec, dcorrPower(8,num_noise_test_min:num_noise_test_max), '+-.', ...
     noiseVec, micePower(8,num_noise_test_min:num_noise_test_max), 'd-.', ...
     noiseVec, corrPower(8,num_noise_test_min:num_noise_test_max), '^-.', ...
     noiseVec, rdcPower(8,num_noise_test_min:num_noise_test_max), 'v-.');  
axis([min(noiseVec) max(noiseVec) 0 1]);
h8.FontSize = 20; 
legend('RSDM', 'dcorr', 'MIC_e', 'corr', 'RDC');  % manually move this using the mouse to a
                                                  % good location
xlabel('Noise Level'); ylabel('Power'); grid on;
loc_inset = [h8.Position(1)+inset_bufX h8.Position(2)+inset_bufY inset_width inset_height];
ax8 = axes('Position',loc_inset);
tmp1 = linspace(0,1,M_inlet);
tmp2 = (tmp1 > 0.5);
plot(tmp1,tmp2, 'k', 'LineWidth', 2);
ax8.Box = 'on'; ax8.XTick = []; ax8.YTick = [];
ax8.XLim = [min(tmp1) max(tmp1)];

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
miceNull = zeros(1,nsim_null);
corrNull = zeros(1,nsim_null);
rdcNull  = zeros(1,nsim_null);

rsdmAlt  = zeros(1,nsim_alt);
dcorrAlt = zeros(1,nsim_alt);
miceAlt = zeros(1,nsim_alt);
corrAlt = zeros(1,nsim_alt);
rdcAlt  = zeros(1,nsim_alt);

% Arrays holding the estimated power for each of the "correlation" types, 
% for each data type (linear, parabolic, etc...) with each noise level
rsdmPower = zeros(numDepTests, num_noise, length(M_vec));
dcorrPower = zeros(numDepTests, num_noise, length(M_vec));
micePower = zeros(numDepTests, num_noise, length(M_vec));
corrPower = zeros(numDepTests, num_noise, length(M_vec));
rdcPower  = zeros(numDepTests, num_noise, length(M_vec));

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
                miceNull(ii) = minestats.mic;
                % compute correlation
                corrNull(ii) = corr(x,y);
                % compute RDC
                rdcNull(ii) = rdc(x,y,rdc_k,rdc_s);
            end

            % compute the rejection cutoffs
            rsdm_cut = quantile(rsdmNull, 0.95);
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
                rsdmAlt(ii) = rsdm(x, y);
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
            rsdmPower(typ, l, m)  = sum(rsdmAlt > rsdm_cut)/nsim_alt;
            dcorrPower(typ, l, m)  = sum(dcorrAlt > dcorr_cut)/nsim_alt;
            micePower(typ, l, m)   = sum(miceAlt > mice_cut)/nsim_alt;
            corrPower(typ, l, m)   = sum(corrAlt > corr_cut)/nsim_alt;
            rdcPower(typ, l, m)    = sum(rdcAlt > rdc_cut)/nsim_alt;
        end
    end
end
% save the data
if(ispc)
    save('C:\\Users\\Kiran\\ownCloud\\PhD\\sim_results\\independence\\rsdmPower.mat');
elseif(ismac)
    save('/Users/Kiran/ownCloud/PhD/sim_results/independence/rsdmPower.mat');
else
    save('/home/kiran/ownCloud/PhD/sim_results/independence/rsdmPower.mat');
end

% process these results to search for the minimum M that gives us the
% threshold power
depMetrics = {'rsdm', 'dcorr', 'mice', 'corr', 'rdc'};
sampleSizeAnalysisVec = zeros(numDepTests, length(depMetrics), length(num_noise_test_min:num_noise_test_max));
for depMetricIdx=1:length(depMetrics)
    if(strcmpi(depMetrics{depMetricIdx}, 'rsdm'))
        powerData = rsdmPower;
    elseif(strcmpi(depMetrics{depMetricIdx}, 'dcorr'))
        powerData = dcorrPower;
    elseif(strcmpi(depMetrics{depMetricIdx}, 'mice'))
        powerData = micePower;
    elseif(strcmpi(depMetrics{depMetricIdx}, 'corr'))
        powerData = corrPower;
    elseif(strcmpi(depMetrics{depMetricIdx}, 'rdc'))
        powerData = rdcPower;
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

% define indices to get to the correct dep metric easily
rsdmIdx = 1;
dcorrIdx = 2;
miceIdx = 3;
corrIdx = 4;
rdcIdx = 5;

% inlet plot configuration
M_inlet = 200;
inset_bufX = 0.0005; inset_bufY = 0.26;
inset_width = 0.1; inset_height = 0.08;

noiseVec = (num_noise_test_min:num_noise_test_max)/10;
figure;
h1 = subplot(2,2,1);
hh1 = plot(noiseVec, squeeze(sampleSizeAnalysisVec(1,rsdmIdx,num_noise_test_min:num_noise_test_max)), 'o-.', ...
     noiseVec, squeeze(sampleSizeAnalysisVec(1,dcorrIdx,num_noise_test_min:num_noise_test_max)), '+-.', ...
     noiseVec, squeeze(sampleSizeAnalysisVec(1,miceIdx,num_noise_test_min:num_noise_test_max)), 'd-.', ...
     noiseVec, squeeze(sampleSizeAnalysisVec(1,corrIdx,num_noise_test_min:num_noise_test_max)), '^-.', ...
     noiseVec, squeeze(sampleSizeAnalysisVec(1,rdcIdx,num_noise_test_min:num_noise_test_max)), 'v-.'); 
% axis([min(noiseVec) max(noiseVec) 0 1]);
xlabel('Noise Level', 'FontSize', 20); ylabel('min(Sample Size)', 'FontSize', 20); grid on;
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
hh2 = plot(noiseVec, squeeze(sampleSizeAnalysisVec(2,rsdmIdx,num_noise_test_min:num_noise_test_max)), 'o-.', ...
     noiseVec, squeeze(sampleSizeAnalysisVec(2,dcorrIdx,num_noise_test_min:num_noise_test_max)), '+-.', ...
     noiseVec, squeeze(sampleSizeAnalysisVec(2,miceIdx,num_noise_test_min:num_noise_test_max)), 'd-.', ...
     noiseVec, squeeze(sampleSizeAnalysisVec(2,corrIdx,num_noise_test_min:num_noise_test_max)), '^-.', ...
     noiseVec, squeeze(sampleSizeAnalysisVec(2,rdcIdx,num_noise_test_min:num_noise_test_max)), 'v-.'); 
% axis([min(noiseVec) max(noiseVec) 0 1]);
xlabel('Noise Level', 'FontSize', 20); ylabel('min(Sample Size)', 'FontSize', 20); grid on;
h2.FontSize = 20; 
loc_inset = [h2.Position(1)+inset_bufX h2.Position(2)+inset_bufY inset_width inset_height];
ax1 = axes('Position',loc_inset);
tmp1 = linspace(0,1,M_inlet);
tmp2 = 4*(tmp1-.5).^2;
plot(tmp1,tmp2, 'k', 'LineWidth', 2);
ax1.Box = 'on'; ax1.XTick = []; ax1.YTick = [];
ax1.XLim = [min(tmp1) max(tmp1)];
ax1.YLim = [min(tmp2) max(tmp2)];

h3 = subplot(2,2,3);
hh3 = plot(noiseVec, squeeze(sampleSizeAnalysisVec(3,rsdmIdx,num_noise_test_min:num_noise_test_max)), 'o-.', ...
     noiseVec, squeeze(sampleSizeAnalysisVec(3,dcorrIdx,num_noise_test_min:num_noise_test_max)), '+-.', ...
     noiseVec, squeeze(sampleSizeAnalysisVec(3,miceIdx,num_noise_test_min:num_noise_test_max)), 'd-.', ...
     noiseVec, squeeze(sampleSizeAnalysisVec(3,corrIdx,num_noise_test_min:num_noise_test_max)), '^-.', ...
     noiseVec, squeeze(sampleSizeAnalysisVec(3,rdcIdx,num_noise_test_min:num_noise_test_max)), 'v-.'); 
% axis([min(noiseVec) max(noiseVec) 0 1]);
xlabel('Noise Level', 'FontSize', 20); ylabel('min(Sample Size)', 'FontSize', 20); grid on;
h3.FontSize = 20; 
loc_inset = [h3.Position(1)+inset_bufX h3.Position(2)+inset_bufY inset_width inset_height];
ax1 = axes('Position',loc_inset);
tmp1 = linspace(0,1,M_inlet);
tmp2 = 128*(tmp1-1/3).^3-48*(tmp1-1/3).^3-12*(tmp1-1/3);
plot(tmp1,tmp2, 'k', 'LineWidth', 2);
ax1.Box = 'on'; ax1.XTick = []; ax1.YTick = [];
ax1.XLim = [min(tmp1) max(tmp1)];
ax1.YLim = [min(tmp2) max(tmp2)];

h4 = subplot(2,2,4);
hh4 = plot(noiseVec, squeeze(sampleSizeAnalysisVec(4,rsdmIdx,num_noise_test_min:num_noise_test_max)), 'o-.', ...
     noiseVec, squeeze(sampleSizeAnalysisVec(4,dcorrIdx,num_noise_test_min:num_noise_test_max)), '+-.', ...
     noiseVec, squeeze(sampleSizeAnalysisVec(4,miceIdx,num_noise_test_min:num_noise_test_max)), 'd-.', ...
     noiseVec, squeeze(sampleSizeAnalysisVec(4,corrIdx,num_noise_test_min:num_noise_test_max)), '^-.', ...
     noiseVec, squeeze(sampleSizeAnalysisVec(4,rdcIdx,num_noise_test_min:num_noise_test_max)), 'v-.'); 
% axis([min(noiseVec) max(noiseVec) 0 1]);
xlabel('Noise Level', 'FontSize', 20); ylabel('min(Sample Size)', 'FontSize', 20); grid on;
h4.FontSize = 20; 
loc_inset = [h4.Position(1)+inset_bufX h4.Position(2)+inset_bufY inset_width inset_height];
ax1 = axes('Position',loc_inset);
tmp1 = linspace(0,1,M_inlet);
tmp2 = sin(4*pi*tmp1);
plot(tmp1,tmp2, 'k', 'LineWidth', 2);
ax1.Box = 'on'; ax1.XTick = []; ax1.YTick = [];
ax1.XLim = [min(tmp1) max(tmp1)];
ax1.YLim = [min(tmp2) max(tmp2)];

figure;
h5 = subplot(2,2,1);
hh5 = plot(noiseVec, squeeze(sampleSizeAnalysisVec(5,rsdmIdx,num_noise_test_min:num_noise_test_max)), 'o-.', ...
     noiseVec, squeeze(sampleSizeAnalysisVec(5,dcorrIdx,num_noise_test_min:num_noise_test_max)), '+-.', ...
     noiseVec, squeeze(sampleSizeAnalysisVec(5,miceIdx,num_noise_test_min:num_noise_test_max)), 'd-.', ...
     noiseVec, squeeze(sampleSizeAnalysisVec(5,corrIdx,num_noise_test_min:num_noise_test_max)), '^-.', ...
     noiseVec, squeeze(sampleSizeAnalysisVec(5,rdcIdx,num_noise_test_min:num_noise_test_max)), 'v-.'); 
% axis([min(noiseVec) max(noiseVec) 0 1]);
xlabel('Noise Level', 'FontSize', 20); ylabel('min(Sample Size)', 'FontSize', 20); grid on;
h5.FontSize = 20; 
loc_inset = [h5.Position(1)+inset_bufX h5.Position(2)+inset_bufY inset_width inset_height];
ax1 = axes('Position',loc_inset);
tmp1 = linspace(0,1,M_inlet);
tmp2 = sin(16*pi*tmp1);
plot(tmp1,tmp2, 'k', 'LineWidth', 2);
ax1.Box = 'on'; ax1.XTick = []; ax1.YTick = [];
ax1.XLim = [min(tmp1) max(tmp1)];
ax1.YLim = [min(tmp2) max(tmp2)];

h6 = subplot(2,2,2);
hh6 = plot(noiseVec, squeeze(sampleSizeAnalysisVec(6,rsdmIdx,num_noise_test_min:num_noise_test_max)), 'o-.', ...
     noiseVec, squeeze(sampleSizeAnalysisVec(6,dcorrIdx,num_noise_test_min:num_noise_test_max)), '+-.', ...
     noiseVec, squeeze(sampleSizeAnalysisVec(6,miceIdx,num_noise_test_min:num_noise_test_max)), 'd-.', ...
     noiseVec, squeeze(sampleSizeAnalysisVec(6,corrIdx,num_noise_test_min:num_noise_test_max)), '^-.', ...
     noiseVec, squeeze(sampleSizeAnalysisVec(6,rdcIdx,num_noise_test_min:num_noise_test_max)), 'v-.'); 
% axis([min(noiseVec) max(noiseVec) 0 1]);
xlabel('Noise Level', 'FontSize', 20); ylabel('min(Sample Size)', 'FontSize', 20); grid on;
h6.FontSize = 20; 
loc_inset = [h6.Position(1)+inset_bufX h6.Position(2)+inset_bufY inset_width inset_height];
ax1 = axes('Position',loc_inset);
tmp1 = linspace(0,1,M_inlet);
tmp2 = tmp1.^(1/4);
plot(tmp1,tmp2, 'k', 'LineWidth', 2);
ax1.Box = 'on'; ax1.XTick = []; ax1.YTick = [];
ax1.XLim = [min(tmp1) max(tmp1)];
ax1.YLim = [min(tmp2) max(tmp2)];

h7 = subplot(2,2,3);
hh7 = plot(noiseVec, squeeze(sampleSizeAnalysisVec(7,rsdmIdx,num_noise_test_min:num_noise_test_max)), 'o-.', ...
     noiseVec, squeeze(sampleSizeAnalysisVec(7,dcorrIdx,num_noise_test_min:num_noise_test_max)), '+-.', ...
     noiseVec, squeeze(sampleSizeAnalysisVec(7,miceIdx,num_noise_test_min:num_noise_test_max)), 'd-.', ...
     noiseVec, squeeze(sampleSizeAnalysisVec(7,corrIdx,num_noise_test_min:num_noise_test_max)), '^-.', ...
     noiseVec, squeeze(sampleSizeAnalysisVec(7,rdcIdx,num_noise_test_min:num_noise_test_max)), 'v-.'); 
% axis([min(noiseVec) max(noiseVec) 0 1]);
xlabel('Noise Level', 'FontSize', 20); ylabel('min(Sample Size)', 'FontSize', 20); grid on;
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

h8 = subplot(2,2,4);
hh8 = plot(noiseVec, squeeze(sampleSizeAnalysisVec(8,rsdmIdx,num_noise_test_min:num_noise_test_max)), 'o-.', ...
     noiseVec, squeeze(sampleSizeAnalysisVec(8,dcorrIdx,num_noise_test_min:num_noise_test_max)), '+-.', ...
     noiseVec, squeeze(sampleSizeAnalysisVec(8,miceIdx,num_noise_test_min:num_noise_test_max)), 'd-.', ...
     noiseVec, squeeze(sampleSizeAnalysisVec(8,corrIdx,num_noise_test_min:num_noise_test_max)), '^-.', ...
     noiseVec, squeeze(sampleSizeAnalysisVec(8,rdcIdx,num_noise_test_min:num_noise_test_max)), 'v-.'); 
legend('RSDM', 'dcorr', 'MIC_e', 'corr', 'RDC');  % manually move this using the mouse to a
                                                  % good location
% axis([min(noiseVec) max(noiseVec) 0 1]);
xlabel('Noise Level', 'FontSize', 20); ylabel('min(Sample Size)', 'FontSize', 20); grid on;
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

% save the data
if(ispc)
    save('C:\\Users\\Kiran\\ownCloud\\PhD\\sim_results\\independence\\rsdmNullDistribution.mat');
elseif(ismac)
    save('/Users/Kiran/ownCloud/PhD/sim_results/independence/rsdmNullDistribution.mat');
else
    save('/home/kiran/ownCloud/PhD/sim_results/independence/rsdmNullDistribution.mat');
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

% From the above analysis, the Inverse Gaussian distribution seems to fit
% best ... Do Q-Q Plots of the Inverse Gaussian Distribution Fit

% QQ Plot w/ best fit for M=100 and M=1000
pdObjs = cell(1,length(M_vec));
muVec = zeros(1,length(M_vec));
lambdaVec = zeros(1,length(M_vec));
for ii=1:length(M_vec)
    M = M_vec(ii);
    % look for the Inverse Gaussian Distribution in the correct cell array
    D = D_cell(ii);
    PD = PD_cell(ii); PD = PD{1};
    for jj=1:length(PD)
        if(strcmpi('Inverse Gaussian', PD{jj}.DistributionName))
            pd = PD{jj};
        end
    end
    pdObjs{ii} = pd;
    muVec(ii) = pd.mu;
    lambdaVec(ii) = pd.lambda;
end

% do the Q-Q plot
pd = pdObjs{1};
subplot(2,2,1); qqplot(rsdmNullDistributionResults(:,1), pd); grid on;
xlabel(sprintf('Quantiles of IG(%0.02f,%0.02f)', pd.mu, pd.lambda), 'FontSize', 20);
ylabel('Quantiles of Input Samples', 'FontSize', 20);
title('M = 100', 'FontSize', 24);

pd = pdObjs{10};
subplot(2,2,3); qqplot(rsdmNullDistributionResults(:,1), pd); grid on;
xlabel(sprintf('Quantiles of IG(%0.02f,%0.02f)', pd.mu, pd.lambda), 'FontSize', 20);
ylabel('Quantiles of Input Samples', 'FontSize', 20);
title('M = 1000', 'FontSize', 24);

% plot how mu and lambda change as M goes from 100 --> 1000
subplot(2,2,2); plot(M_vec, muVec); grid on; xlabel('M', 'FontSize', 20); ylabel('\mu', 'FontSize', 20); grid on;
subplot(2,2,4); plot(M_vec, lambdaVec); grid on; xlabel('M', 'FontSize', 20); ylabel('\lambda', 'FontSize', 20); grid on;