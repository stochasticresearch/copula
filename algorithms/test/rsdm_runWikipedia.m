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
rhos = [0.8, 0.4, 0.0, -0.4, -0.8];

nCases = 18;
numDiscreteIntervals = 6;

% configure RDC
rdc_k = 20;
rdc_s = 1/6;

% configure MIC_e
mine_c = 15;
mine_alpha = 0.6;

rng(2);

% do all the gaussian ones
xy = cell(1,nCases);
xxyy = cell(1,nCases);
xyIdx = 1;

xContinuous = rand(n, 1)*2-1; % uniform between +/- 1
yContinuous = rand(n, 1)*2-1;
xx = discretizeRv(xContinuous,numDiscreteIntervals)';
yy = discretizeRv(yContinuous,numDiscreteIntervals)';
xy{xyIdx} = [xContinuous yContinuous]; xxyy{xyIdx} = [xx yy];
xyIdx = xyIdx + 1;

rsdmBias = rsdm(xContinuous, yContinuous);
rdcBias = rdc(xContinuous,yContinuous,rdc_k,rdc_s);
dcorrBias = dcorr(xContinuous, yContinuous);
corrBias = corr(xContinuous,yContinuous);
minestatsContinuous = mine(xContinuous',yContinuous',mine_alpha,mine_c,'mic_e');
mic_e_bias = minestatsContinuous.mic;

xContinuous = rand(n,1)*2-1;
yContinuous = xContinuous;
dat = [xContinuous yContinuous];
xy{xyIdx} = dat; xxyy{xyIdx} = [discretizeRv(dat(:,1),numDiscreteIntervals)' discretizeRv(dat(:,2),numDiscreteIntervals)'];
xyIdx = xyIdx + 1;

xContinuous = rand(n,1)*2-1;
yContinuous = -1*xContinuous;
dat = [xContinuous yContinuous];
xy{xyIdx} = dat; xxyy{xyIdx} = [discretizeRv(dat(:,1),numDiscreteIntervals)' discretizeRv(dat(:,2),numDiscreteIntervals)'];
xyIdx = xyIdx + 1;

for rho=rhos
    dat = mvnrnd([0 0], [1 rho; rho 1], n);
    xy{xyIdx} = dat; xxyy{xyIdx} = [discretizeRv(dat(:,1),numDiscreteIntervals)' discretizeRv(dat(:,2),numDiscreteIntervals)'];
    xyIdx = xyIdx + 1;
end

% do the "others" now
spread = 0;

xContinuous = rand(n,1)*2-1;
yContinuous = rand(n,1)*2-1;
t = -pi/8;
rotationMat = [cos(t) sin(t); -sin(t) cos(t)];
dat = [xContinuous yContinuous]*rotationMat; 
xy{xyIdx} = dat; xxyy{xyIdx} = [discretizeRv(dat(:,1),numDiscreteIntervals)' discretizeRv(dat(:,2),numDiscreteIntervals)'];
xyIdx = xyIdx + 1;

xContinuous = rand(n,1)*2-1;
yContinuous = rand(n,1)*2-1;
t = -pi/4;
rotationMat = [cos(t) sin(t); -sin(t) cos(t)];
dat = [xContinuous yContinuous]*rotationMat; xyIdx = xyIdx + 1;
xy{xyIdx} = dat; xxyy{xyIdx} = [discretizeRv(dat(:,1),numDiscreteIntervals)' discretizeRv(dat(:,2),numDiscreteIntervals)'];
xyIdx = xyIdx + 1;

xContinuous = rand(n,1)*2-1;
yContinuous = 2*xContinuous.^2 + (rand(n,1)*2-1)*spread;
dat = [xContinuous yContinuous];
xy{xyIdx} = dat; xxyy{xyIdx} = [discretizeRv(dat(:,1),numDiscreteIntervals)' discretizeRv(dat(:,2),numDiscreteIntervals)'];
xyIdx = xyIdx + 1;

xContinuous = rand(n,1)*2-1;
yContinuous = (xContinuous.^2 + rand(n,1)/2 * spread) .* (2*binornd(1,0.5,n,1)-1);
dat = [xContinuous yContinuous];
xy{xyIdx} = dat; xxyy{xyIdx} = [discretizeRv(dat(:,1),numDiscreteIntervals)' discretizeRv(dat(:,2),numDiscreteIntervals)'];
xyIdx = xyIdx + 1;

xx = (rand(n,1)*2-1);
yContinuous = cos(xx*pi) + randn(n,1)*1/8 * spread;
xContinuous = sin(xx*pi) + randn(n,1)*1/8 * spread;
dat = [xContinuous yContinuous];
xy{xyIdx} = dat; xxyy{xyIdx} = [discretizeRv(dat(:,1),numDiscreteIntervals)' discretizeRv(dat(:,2),numDiscreteIntervals)'];
xyIdx = xyIdx + 1;

spread = 0.1;
xy1 = mvnrnd([3 3], [1 spread; spread 1], n/4);
xy2 = mvnrnd([-3 3], [1 spread; spread 1], n/4);
xy3 = mvnrnd([3 -3], [1 spread; spread 1], n/4);
xy4 = mvnrnd([-3 -3], [1 spread; spread 1], n/4);
dat = [xy1; xy2; xy3; xy4];
xy{xyIdx} = dat; xxyy{xyIdx} = dat;     % TODO: fix discrete!

xContinuous = rand(n,1);
yContinuous = sin(4*pi*xContinuous);
dat = [xContinuous yContinuous];
xy{xyIdx} = dat; xxyy{xyIdx} = [discretizeRv(dat(:,1),numDiscreteIntervals)' discretizeRv(dat(:,2),numDiscreteIntervals)'];
xyIdx = xyIdx + 1;

xContinuous = rand(n,1);
yContinuous = exp(xContinuous);
dat = [xContinuous yContinuous];
xy{xyIdx} = dat; xxyy{xyIdx} = [discretizeRv(dat(:,1),numDiscreteIntervals)' discretizeRv(dat(:,2),numDiscreteIntervals)'];
xyIdx = xyIdx + 1;

xContinuous = rand(n,1)*3;
yContinuous = exp(-xContinuous).*cos(pi*xContinuous);
% x = rand(n,1);
% y = 128*(x-1/3).^3-48*(x-1/3).^3-12*(x-1/3);
dat = [xContinuous yContinuous];
xy{xyIdx} = dat; xxyy{xyIdx} = [discretizeRv(dat(:,1),numDiscreteIntervals)' discretizeRv(dat(:,2),numDiscreteIntervals)'];
xyIdx = xyIdx + 1;

xContinuous = rand(n,1);
yContinuous = xContinuous.^(1/4);
dat = [xContinuous yContinuous];
xy{xyIdx} = dat; xxyy{xyIdx} = [discretizeRv(dat(:,1),numDiscreteIntervals)' discretizeRv(dat(:,2),numDiscreteIntervals)'];
xyIdx = xyIdx + 1;

%% Plot with barcharts and inset
clc;
figure;
nCases = 9;
dataIdx = [2 1 3 4 6 8 13 12 15];

inset_width = 0.075; inset_height = 0.075;
% inset_bufX = 0.05; inset_bufY = 0.02; % for bottom left inset
inset_bufX = 0.17; inset_bufY = 0.15;   % for top right inset

for ii=1:nCases
    h = subplot(3,3,ii);

    dataContinuous = xy{dataIdx(ii)};
    xContinuous = dataContinuous(:,1); yContinuous = dataContinuous(:,2);
    dataDiscrete = xxyy{dataIdx(ii)};
    xDiscrete = dataDiscrete(:,1); yDiscrete = dataDiscrete(:,2);
    
    rsdmValContinuous = rsdm(xContinuous, yContinuous);
    rdcValContinuous = rdc(xContinuous,yContinuous,rdc_k,rdc_s);
    minestatsContinuous = mine(xContinuous',yContinuous',mine_alpha,mine_c,'mic_e');
    mic_e_valContinuous = minestatsContinuous.mic;
    dcorrValContinuous = dcorr(xContinuous,yContinuous);
    corrValContinuous = corr(xContinuous,yContinuous);
    cosValContinuous = cosdv(xContinuous,yContinuous);
    ccorrValContinuous = cCorr(xContinuous,yContinuous);
    
    rsdmValDiscrete = rsdm(xDiscrete, yDiscrete);
    rdcValDiscrete = rdc(xDiscrete,yDiscrete,rdc_k,rdc_s);
    minestatsDiscrete = mine(xDiscrete',yDiscrete',mine_alpha,mine_c,'mic_e');
    mic_e_valDiscrete = minestatsDiscrete.mic;
    dcorrValDiscrete = dcorr(xDiscrete,yDiscrete);
    corrValDiscrete = corr(xDiscrete,yDiscrete);
%     cosValDiscrete = cosdv(xDiscrete,yDiscrete);
    cosValDiscrete = 0;     % glitches b/c it is only meant for continuous
%     ccorrValDiscrete = cCorr(xDiscrete,yDiscrete);
    ccorrValDiscrete = 0;   % glitches b/c it is only meant for discrete
    
    barVals = [rsdmValContinuous rsdmValDiscrete; ...
               rdcValContinuous rdcValDiscrete; ...
               dcorrValContinuous dcorrValDiscrete; ...
               mic_e_valContinuous mic_e_valDiscrete; ...
               cosValContinuous cosValDiscrete; ...
               ccorrValContinuous ccorrValDiscrete];
    
    b = bar(barVals);
    b(1).BarWidth = 1;
    ylim([0 1])
    Labels = {'RSDM', 'RDC', 'MIC_e', 'dCor', 'CoS', 'cCor'};
    set(gca, 'XTick', 1:6, 'XTickLabel', Labels, 'FontSize', 20);
    rotateXLabels( gca(), 80 )
   
    if(ii==nCases)    
        legend({'Continuous', 'Discrete'});
    end
    
    loc_inset = [h.Position(1)+inset_bufX h.Position(2)+inset_bufY inset_width inset_height];
    ax = axes('Position',loc_inset);
    scatter(xContinuous,yContinuous, 'k');
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
    dataContinuous = xy{dataIdx(ii)};
    xContinuous = dataContinuous(:,1); yContinuous = dataContinuous(:,2);
    scatter(xContinuous, yContinuous);
    handles{ii} = h;
    
    % remove axis labels
    set(h,'XTick',[],'YTick',[],'XColor','w','YColor','w','box','off')
    
    % compute the dependency and put as title
    rsdmValContinuous = rsdm(xContinuous, yContinuous, rsdm1_minscanincr, rsdm1_diffthresh, rsdm1_alpha);
    rdcValContinuous = rdc(xContinuous,yContinuous,rdc_k,rdc_s);
    minestatsContinuous = mine(xContinuous',yContinuous',mine_alpha,mine_c,'mic_e');
    mic_e_valContinuous = minestatsContinuous.mic;
    dcorrValContinuous = dcorr(xContinuous,yContinuous);
    corrValContinuous = corr(xContinuous,yContinuous);

    if(ii==1)
        rsdmPrint = rsdmValContinuous;
        rdcPrint = rdcValContinuous;
        dcorrPrint = dcorrValContinuous;
        mic_e_print = mic_e_valContinuous;
        corrPrint = corrValContinuous;
    else
        rsdmPrint = rsdmValContinuous-rsdmBias;
        rdcPrint = rdcValContinuous-rdcBias;
        dcorrPrint = mic_e_valContinuous-mic_e_bias;
        mic_e_print = mic_e_valContinuous-mic_e_bias;
        corrPrint = corrValContinuous - corrBias;
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
    dataContinuous = xy{dataIdx(ii)};
    xContinuous = dataContinuous(:,1); yContinuous = dataContinuous(:,2);
    scatter(xContinuous, yContinuous);
    handles{ii} = h;
    
    % remove axis labels
    set(h,'XTick',[],'YTick',[],'XColor','w','YColor','w','box','off')
    
    % compute the dependency and put as title
    rsdmValContinuous = rsdm(xContinuous, yContinuous, rsdm1_minscanincr, rsdm1_diffthresh, rsdm1_alpha);
    rdcValContinuous = rdc(xContinuous,yContinuous,rdc_k,rdc_s);
    minestatsContinuous = mine(xContinuous',yContinuous',mine_alpha,mine_c,'mic_e');
    mic_e_valContinuous = minestatsContinuous.mic;
    dcorrValContinuous = dcorr(xContinuous,yContinuous);
    corrValContinuous = corr(xContinuous,yContinuous);

    if(ii==1)
        rsdmPrint = rsdmValContinuous;
        rdcPrint = rdcValContinuous;
        dcorrPrint = dcorrValContinuous;
        mic_e_print = mic_e_valContinuous;
        corrPrint = corrValContinuous;
    else
        rsdmPrint = rsdmValContinuous-rsdmBias;
        rdcPrint = rdcValContinuous-rdcBias;
        dcorrPrint = mic_e_valContinuous-mic_e_bias;
        mic_e_print = mic_e_valContinuous-mic_e_bias;
        corrPrint = corrValContinuous - corrBias;
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
