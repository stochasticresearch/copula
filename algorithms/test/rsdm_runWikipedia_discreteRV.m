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
nCases = 9;

% configure RSDM
rsdm_minscanincr = 0.025;
rsdm_diffthresh = 100;
rsdm_alpha = 0.08;

% configure RDC
rdc_k = 20;
rdc_s = 1/6;

% configure MIC_e
mine_c = 15;
mine_alpha = 0.6;

numDiscreteIntervals = 5;

rng(2);

% do all the gaussian ones
xy = cell(1,nCases);

x = unidrnd(numDiscreteIntervals, n, 1);
y = unidrnd(numDiscreteIntervals, n, 1);
x = discretizeRv(x, numDiscreteIntervals)';
y = discretizeRv(y, numDiscreteIntervals)';
xy{1} = [x y];

rsdmBias = rsdm(x, y, rsdm_minscanincr, rsdm_diffthresh, rsdm_alpha);
rdcBias = rdc(x,y,rdc_k,rdc_s);
dcorrBias = dcorr(x, y);
corrBias = corr(x,y);
minestats = mine(x',y',mine_alpha,mine_c,'mic_e');
mic_e_bias = minestats.mic;

% data = mvnrnd([0 0], [1 1; 1 1], n);
% x = data(:,1); y = data(:,2);
x = unidrnd(numDiscreteIntervals, n, 1);
y = x;
xy{2} = [x y];

% data = mvnrnd([0 0], [1 -1; -1 1], n);
% x = data(:,1); y = data(:,2);
x = unidrnd(numDiscreteIntervals, n, 1);
y = -x;
xy{3} = [x y];

x = rand(n,1)*2-1;
x = discretizeRv(x, numDiscreteIntervals)';
y = x.^3;
xy{4} = [x y];

x = rand(n,1)*2-1;
x = discretizeRv(x, numDiscreteIntervals)';
y = 2*x.^2;
xy{5} = [x y];

x = rand(n,1);
x = discretizeRv(x, numDiscreteIntervals)';
y = sin(2*pi*x);
xy{6} = [x y];

x = rand(n,1);
x = discretizeRv(x, numDiscreteIntervals)';
y = exp(x);
xy{7} = [x y];

x = rand(n,1)*3;
x = discretizeRv(x, numDiscreteIntervals)';
y = exp(-x).*cos(pi*x);
xy{8} = [x y];

x = rand(n,1);
x = discretizeRv(x, numDiscreteIntervals)';
y = x.^(1/4);
xy{9} = [x y];

%% Plot with barcharts and inset
clc;
figure;
nCases = 9;
dataIdx = [2 1 3 4 5 6 7 8 9];

inset_width = 0.075; inset_height = 0.075;
inset_bufX = 0.05; inset_bufY = 0.02;

for ii=1:nCases
    h = subplot(3,3,ii);
    h.Position
    data = xy{dataIdx(ii)};
    x = data(:,1); y = data(:,2);
    
    rsdmVal = rsdm(x, y, rsdm_minscanincr, rsdm_diffthresh, rsdm_alpha);
    rdcVal = rdc(x,y,rdc_k,rdc_s);
    minestats = mine(x',y',mine_alpha,mine_c,'mic_e');
    mic_e_val = minestats.mic;
    dcorrVal = dcorr(x,y);
    corrVal = corr(x,y);
    
    rsdmPrint = rsdmVal-rsdmBias;
    rdcPrint = rdcVal-rdcBias;
    dcorrPrint = dcorrVal-dcorrBias;
    mic_e_print = mic_e_val-mic_e_bias;
    corrPrint = abs(corrVal) - corrBias;
    
    % show the bias for the independence plot
    if(dataIdx(ii)==1)
        rsdmPrint = rsdmBias;
        rdcPrint = rdcBias;
        dcorrPrint = dcorrBias;
        mic_e_print = mic_e_bias;
        corrPrint = corrBias;
    end
    
    b = bar([rsdmPrint rdcPrint dcorrPrint mic_e_print corrPrint]);
    b.BarWidth = 0.6;
    ylim([0 1])
    Labels = {'RSDM', 'RDC', 'MIC_e', 'dCorr', '|corr|'};
    set(gca, 'XTick', 1:5, 'XTickLabel', Labels, 'FontSize', 20);
    rotateXLabels( gca(), 80 )
    
%     loc_inset = [h.Position(1)+h.Position(3)-inset_bufX h.Position(2)+h.Position(4)-inset_bufY inset_width inset_height];
    loc_inset = [h.Position(1)-inset_bufX h.Position(2)-inset_bufY inset_width inset_height];
    ax = axes('Position',loc_inset);
    scatter(x,y, 'k', 'filled');
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
    rsdmVal = rsdm(x, y, rsdm_minscanincr, rsdm_diffthresh, rsdm_alpha);
    rdcVal = rdc(x,y,rdc_k,rdc_s);
    minestats = mine(x',y',mine_alpha,mine_c,'mic_e');
    mic_e_val = minestats.mic;
    dcorrVal = dcorr(x,y);
    corrVal = corr(x,y);

    if(ii==1)
        rsdmPrint = rsdmVal;
        rdcPrint = rdcVal;
        dcorrPrint = dcorrVal;
        mic_e_print = mic_e_val;
        corrPrint = corrVal;
    else
        rsdmPrint = rsdmVal-rsdmBias;
        rdcPrint = rdcVal-rdcBias;
        dcorrPrint = mic_e_val-mic_e_bias;
        mic_e_print = mic_e_val-mic_e_bias;
        corrPrint = corrVal - corrBias;
    end
    
    title({[strcat('\fontsize{22} {\color{blue}', sprintf('%1.2f}|',rsdmPrint)), ...
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
    rsdmVal = rsdm(x, y, rsdm_minscanincr, rsdm_diffthresh, rsdm_alpha);
    rdcVal = rdc(x,y,rdc_k,rdc_s);
    minestats = mine(x',y',mine_alpha,mine_c,'mic_e');
    mic_e_val = minestats.mic;
    dcorrVal = dcorr(x,y);
    corrVal = corr(x,y);

    if(ii==1)
        rsdmPrint = rsdmVal;
        rdcPrint = rdcVal;
        dcorrPrint = dcorrVal;
        mic_e_print = mic_e_val;
        corrPrint = corrVal;
    else
        rsdmPrint = rsdmVal-rsdmBias;
        rdcPrint = rdcVal-rdcBias;
        dcorrPrint = mic_e_val-mic_e_bias;
        mic_e_print = mic_e_val-mic_e_bias;
        corrPrint = corrVal - corrBias;
    end
    
    title({[strcat('\fontsize{22} {\color{blue}', sprintf('%1.2f}|',rsdmPrint)), ...
           strcat('{\color{red}', sprintf('%1.2f}|', dcorrPrint)), ...       
           strcat('{\color{orange}', sprintf('%1.2f}|', mic_e_print))]; ...
           [strcat('{\color{magenta}', sprintf('%1.2f}|', corrPrint)), ...
            strcat('{\color{darkGreen}', sprintf('%1.2f}', rdcPrint))]});
end

% % reduce border of subplots
% tightfig(f);
% subplotsqueeze(f,1.1);
