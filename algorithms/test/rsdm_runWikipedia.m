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

%% Generates the Wikipedia Graph here: 
% https://en.wikipedia.org/wiki/Correlation_and_dependence#/media/File:Correlation_examples2.svg
% Adapted from: https://github.com/lopezpaz/randomized_dependence_coefficient/blob/master/code/run_wikipedia.r
clear;
clc;

n = 800;
rhos = [1.0, -1.0, 0.8, 0.4, 0.0, -0.4, -0.8];

nCases = 18;

% configure RSDM
rsdm1_minscanincr = 0.025;
rsdm1_diffthresh = 100;
rsdm1_alpha = 0.08;

rsdm2_minscanincr = 0.05;
rsdm2_diffthresh = 120;
rsdm2_alpha = 0.1;

rsdm4_minscanincr = 0.025;
rsdm4_diffthresh = 160;
rsdm4_alpha = 0.8;

% configure RDC
rdc_k = 20;
rdc_s = 1/6;

% configure MIC_e
mine_c = 15;
mine_alpha = 0.6;

rng(2);

% do all the gaussian ones
xy = cell(1,nCases);
xyIdx = 1;

x = rand(n, 1)*2-1; % uniform between +/- 1
y = rand(n, 1)*2-1;
xy{xyIdx} = [x y]; xyIdx = xyIdx + 1;

rsdm1Bias = rsdm_1(x, y, rsdm1_minscanincr, rsdm1_diffthresh, rsdm1_alpha);
rsdm2Bias = rsdm_2(x, y, rsdm2_minscanincr, rsdm2_diffthresh, rsdm2_alpha);
rsdm4Bias = rsdm_4(x, y, rsdm4_minscanincr, rsdm4_diffthresh, rsdm4_alpha);
rdcBias = rdc(x,y,rdc_k,rdc_s);
dcorrBias = dcorr(x, y);
corrBias = corr(x,y);
minestats = mine(x',y',mine_alpha,mine_c,'mic_e');
mic_e_bias = minestats.mic;

for rho=rhos
    xy{xyIdx} = mvnrnd([0 0], [1 rho; rho 1], n);
    xyIdx = xyIdx + 1;
end
% do the "others" now
spread = 0;

x = rand(n,1)*2-1;
y = rand(n,1)*2-1;
t = -pi/8;
rotationMat = [cos(t) sin(t); -sin(t) cos(t)];
xy{xyIdx} = [x y]*rotationMat; xyIdx = xyIdx + 1;

x = rand(n,1)*2-1;
y = rand(n,1)*2-1;
t = -pi/4;
rotationMat = [cos(t) sin(t); -sin(t) cos(t)];
xy{xyIdx} = [x y]*rotationMat; xyIdx = xyIdx + 1;

x = rand(n,1)*2-1;
y = 2*x.^2 + (rand(n,1)*2-1)*spread;
xy{xyIdx} = [x y]; xyIdx = xyIdx + 1;

x = rand(n,1)*2-1;
y = (x.^2 + rand(n,1)/2 * spread) .* (2*binornd(1,0.5,n,1)-1);
xy{xyIdx} = [x y]; xyIdx = xyIdx + 1;

xx = (rand(n,1)*2-1);
y = cos(xx*pi) + randn(n,1)*1/8 * spread;
x = sin(xx*pi) + randn(n,1)*1/8 * spread;
xy{xyIdx} = [x y]; xyIdx = xyIdx + 1;

spread = 0.1;
xy1 = mvnrnd([3 3], [1 spread; spread 1], n/4);
xy2 = mvnrnd([-3 3], [1 spread; spread 1], n/4);
xy3 = mvnrnd([3 -3], [1 spread; spread 1], n/4);
xy4 = mvnrnd([-3 -3], [1 spread; spread 1], n/4);
xy{xyIdx} = [xy1; xy2; xy3; xy4]; xyIdx = xyIdx + 1;

x = rand(n,1);
y = sin(4*pi*x);
xy{xyIdx} = [x y]; xyIdx = xyIdx + 1;

x = rand(n,1);
y = exp(x);
xy{xyIdx} = [x y]; xyIdx = xyIdx + 1;

x = rand(n,1)*3;
y = exp(-x).*cos(pi*x);
% x = rand(n,1);
% y = 128*(x-1/3).^3-48*(x-1/3).^3-12*(x-1/3);
xy{xyIdx} = [x y]; xyIdx = xyIdx + 1;

x = rand(n,1);
y = x.^(1/4);
xy{xyIdx} = [x y]; xyIdx = xyIdx + 1;

%% Plot with barcharts and inset
clc;
figure;
nCases = 9;
dataIdx = [2 1 3 4 6 8 13 12 15];

inset_width = 0.075; inset_height = 0.075;
inset_bufX = 0.05; inset_bufY = 0.02;

for ii=1:nCases
    h = subplot(3,3,ii);
    h.Position
    data = xy{dataIdx(ii)};
    x = data(:,1); y = data(:,2);
    
    rsdm1Val = rsdm_1(x, y, rsdm1_minscanincr, rsdm1_diffthresh, rsdm1_alpha);
    rsdm2Val = rsdm_2(x, y, rsdm2_minscanincr, rsdm2_diffthresh, rsdm2_alpha);
    rsdm4Val = rsdm_4(x, y, rsdm4_minscanincr, rsdm4_diffthresh, rsdm4_alpha);
    rdcVal = rdc(x,y,rdc_k,rdc_s);
    minestats = mine(x',y',mine_alpha,mine_c,'mic_e');
    mic_e_val = minestats.mic;
    dcorrVal = dcorr(x,y);
    corrVal = corr(x,y);
    
    rsdm1Print = rsdm1Val-rsdm1Bias;
    rsdm2Print = rsdm2Val-rsdm2Bias;
    rsdm4Print = rsdm4Val-rsdm4Bias;
    rdcPrint = rdcVal-rdcBias;
    dcorrPrint = dcorrVal-dcorrBias;
    mic_e_print = mic_e_val-mic_e_bias;
    corrPrint = abs(corrVal) - corrBias;
    
    % show the bias for the independence plot
    if(dataIdx(ii)==1)
        rsdm1Print = rsdm1Bias;
        rsdm2Print = rsdm2Bias;
        rsdm4Print = rsdm4Bias;
        rdcPrint = rdcBias;
        dcorrPrint = dcorrBias;
        mic_e_print = mic_e_bias;
        corrPrint = corrBias;
    end
    
    b = bar([rsdm1Print rsdm2Print rsdm4Print rdcPrint dcorrPrint mic_e_print corrPrint]);
    b.BarWidth = 0.6;
    ylim([0 1])
    Labels = {'RSDM1', 'RSDM2', 'RSDM4', 'RDC', 'MIC_e', 'dCorr', '|corr|'};
    set(gca, 'XTick', 1:7, 'XTickLabel', Labels, 'FontSize', 20);
    rotateXLabels( gca(), 80 )
    
%     loc_inset = [h.Position(1)+h.Position(3)-inset_bufX h.Position(2)+h.Position(4)-inset_bufY inset_width inset_height];
    loc_inset = [h.Position(1)-inset_bufX h.Position(2)-inset_bufY inset_width inset_height];
    ax = axes('Position',loc_inset);
    scatter(x,y, 'k');
    ax.Box = 'on'; ax.XTick = []; ax.YTick = [];
    
end

%% Plot w/ title showing the metric
% Plotting section
figure; 
nCases = 9;
handles = cell(1,nCases);
dataIdx = [1 2 3 4 6 8 5 7 14];
for ii=1:nCases   
    h = subplot(3,3,ii);
    data = xy{dataIdx(ii)};
    x = data(:,1); y = data(:,2);
    scatter(x, y);
    handles{ii} = h;
    
    % remove axis labels
    set(h,'XTick',[],'YTick',[],'XColor','w','YColor','w','box','off')
    
    % compute the dependency and put as title
    rsdm1Val = rsdm(x, y, rsdm1_minscanincr, rsdm1_diffthresh, rsdm1_alpha);
    rdcVal = rdc(x,y,rdc_k,rdc_s);
    minestats = mine(x',y',mine_alpha,mine_c,'mic_e');
    mic_e_val = minestats.mic;
    dcorrVal = dcorr(x,y);
    corrVal = corr(x,y);

    if(ii==1)
        rsdm1Print = rsdm1Val;
        rdcPrint = rdcVal;
        dcorrPrint = dcorrVal;
        mic_e_print = mic_e_val;
        corrPrint = corrVal;
    else
        rsdm1Print = rsdm1Val-rsdm1Bias;
        rdcPrint = rdcVal-rdcBias;
        dcorrPrint = mic_e_val-mic_e_bias;
        mic_e_print = mic_e_val-mic_e_bias;
        corrPrint = corrVal - corrBias;
    end
    
    title({[strcat('\fontsize{22} {\color{blue}', sprintf('%1.2f}|',rsdm1Print)), ...
           strcat('{\color{red}', sprintf('%1.2f}|', dcorrPrint)), ...       
           strcat('{\color{orange}', sprintf('%1.2f}|', mic_e_print))]; ...
           [strcat('{\color{magenta}', sprintf('%1.2f}|', corrPrint)), ...
            strcat('{\color{darkGreen}', sprintf('%1.2f}', rdcPrint))]});
end

figure; 
nCases = 9;
handles = cell(1,nCases);
dataIdx = [9 10 11 12 13 15 16 17 18];
for ii=1:nCases   
    h = subplot(3,3,ii);
    data = xy{dataIdx(ii)};
    x = data(:,1); y = data(:,2);
    scatter(x, y);
    handles{ii} = h;
    
    % remove axis labels
    set(h,'XTick',[],'YTick',[],'XColor','w','YColor','w','box','off')
    
    % compute the dependency and put as title
    rsdm1Val = rsdm(x, y, rsdm1_minscanincr, rsdm1_diffthresh, rsdm1_alpha);
    rdcVal = rdc(x,y,rdc_k,rdc_s);
    minestats = mine(x',y',mine_alpha,mine_c,'mic_e');
    mic_e_val = minestats.mic;
    dcorrVal = dcorr(x,y);
    corrVal = corr(x,y);

    if(ii==1)
        rsdm1Print = rsdm1Val;
        rdcPrint = rdcVal;
        dcorrPrint = dcorrVal;
        mic_e_print = mic_e_val;
        corrPrint = corrVal;
    else
        rsdm1Print = rsdm1Val-rsdm1Bias;
        rdcPrint = rdcVal-rdcBias;
        dcorrPrint = mic_e_val-mic_e_bias;
        mic_e_print = mic_e_val-mic_e_bias;
        corrPrint = corrVal - corrBias;
    end
    
    title({[strcat('\fontsize{22} {\color{blue}', sprintf('%1.2f}|',rsdm1Print)), ...
           strcat('{\color{red}', sprintf('%1.2f}|', dcorrPrint)), ...       
           strcat('{\color{orange}', sprintf('%1.2f}|', mic_e_print))]; ...
           [strcat('{\color{magenta}', sprintf('%1.2f}|', corrPrint)), ...
            strcat('{\color{darkGreen}', sprintf('%1.2f}', rdcPrint))]});
end

% % reduce border of subplots
% tightfig(f);
% subplotsqueeze(f,1.1);
