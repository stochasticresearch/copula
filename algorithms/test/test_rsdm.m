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

%% 
clear;
clc;

M = 1000;

x = rand(M,1);
y = x;

% x = rand(M,1);
% y = sin(2*pi*x);

% x = rand(M,1)*2-1;
% y = x.^2;

% x = rand(M,1);
% y=(2*binornd(1,0.5,M,1)-1).* (sqrt(1 - (2*x - 1).^2));

[rsdmMetric, rsdmResidual, residAssocIdxs, rsdmRectangleCfg] = rsdm(x,y);

%% Understand how different tau's work w/ discrete and hybrid functional dependencies
clear;
clc;
close all;

% rng(12345);

M = 500;

numDiscreteIntervals = 4;

% optimal parameters for MICe
mine_c = 15;
mine_alpha = 0.6;

% Optimal parameters for RDC
rdc_k = 20;
rdc_s = 1/6;

% Strictly monotonic
x = rand(M,1);
% y = x.^3;
y = x;
xx = discretizeRv(x,numDiscreteIntervals)';
yy = discretizeRv(y,numDiscreteIntervals)';
% yyy = xx.^3;        % both x and y are discrete here
yyy = xx;

% some required but not used configuration parameters
alpha_dontCare = 0.05;
wantplot_dontCare = 0;

% calculate the different kind's of tau for each of these
% continuous/discrete/hybrid scenarios
% we have 5 measures -- tau, tau_b, tau_cj, mic_e, rdc
dep_xC_yC = zeros(1,5);
dep_xC_yD = zeros(1,5);
dep_xD_yC = zeros(1,5);
dep_xD_yD = zeros(1,5);

tau = corr(x,y,'type','kendall');
taub = ktaub([x y], alpha_dontCare, wantplot_dontCare);
tau_hat = ktauhat(x, y);
minestats = mine(x',y',mine_alpha,mine_c,'mic_e');
rdcVal = rdc(x,y,rdc_k,rdc_s);
dep_xC_yC(1:5) = [tau taub tau_hat minestats.mic rdcVal];

tau = corr(x,yy,'type','kendall');
taub = ktaub([x yy], alpha_dontCare, wantplot_dontCare);
tau_hat = ktauhat(x, yy);
minestats = mine(x',yy',mine_alpha,mine_c,'mic_e');
rdcVal = rdc(x,yy,rdc_k,rdc_s);
dep_xC_yD(1:5) = [tau taub tau_hat minestats.mic rdcVal];

tau = corr(xx,y,'type','kendall');
taub = ktaub([xx y], alpha_dontCare, wantplot_dontCare);
tau_hat = ktauhat(xx, y);
minestats = mine(xx',y',mine_alpha,mine_c,'mic_e');
rdcVal = rdc(xx,y,rdc_k,rdc_s);
dep_xD_yC(1:5) = [tau taub tau_hat minestats.mic rdcVal];

tau = corr(xx,yyy,'type','kendall');
taub = ktaub([xx yyy], alpha_dontCare, wantplot_dontCare);
tau_hat = ktauhat(xx, yyy);
minestats = mine(xx',yyy',mine_alpha,mine_c,'mic_e');
rdcVal = rdc(xx,yyy,rdc_k,rdc_s);
dep_xD_yD(1:5) = [tau taub tau_hat minestats.mic rdcVal];

figure;
subplot(2,2,1);
scatter(x,y); grid on;
title(sprintf('$\\tau=%0.02f \\  \\tau_b=%0.02f \\ \\hat{\\tau}=%0.02f \\ MIC_e=%0.02f \\ RDC=%0.02f$', ...
    dep_xC_yC(1), dep_xC_yC(2), dep_xC_yC(3), dep_xC_yC(4), dep_xC_yC(5) ), 'Interpreter', 'Latex');

subplot(2,2,2);
scatter(x,yy); grid on;
title(sprintf('$\\tau=%0.02f \\ \\tau_b=%0.02f \\ \\hat{\\tau}=%0.02f \\ MIC_e=%0.02f \\  RDC=%0.02f$', ...
    dep_xC_yD(1), dep_xC_yD(2), dep_xC_yD(3), dep_xC_yD(4), dep_xC_yD(5) ), 'Interpreter', 'Latex');

subplot(2,2,3);
scatter(xx,y); grid on;
title(sprintf('$\\tau=%0.02f \\ \\tau_b=%0.02f \\ \\hat{\\tau}=%0.02f \\ MIC_e=%0.02f \\ RDC=%0.02f$', ...
    dep_xD_yC(1), dep_xD_yC(2), dep_xD_yC(3), dep_xD_yC(4), dep_xD_yC(5) ), 'Interpreter', 'Latex');

subplot(2,2,4);
scatter(xx,yyy); grid on;
title(sprintf('$\\tau=%0.02f \\ \\tau_b=%0.02f \\ \\hat{\\tau}=%0.02f \\ MIC_e=%0.02f \\ RDC=%0.02f$', ...
    dep_xD_yD(1), dep_xD_yD(2), dep_xD_yD(3), dep_xD_yD(4), dep_xD_yD(5) ), 'Interpreter', 'Latex');

% Strictly counter-monotonic
% Strictly monotonic
x = rand(M,1);
% y = -x.^3;
y = -x;
xx = discretizeRv(x,numDiscreteIntervals)';
yy = discretizeRv(y,numDiscreteIntervals)';
% yyy = -xx.^3;        % both x and y are discrete here
yyy = -xx;

% some required but not used configuration parameters
alpha_dontCare = 0.05;
wantplot_dontCare = 0;

% calculate the different kind's of tau for each of these
% continuous/discrete/hybrid scenarios
% we have 5 measures -- tau, tau_b, tau_cj, mic_e, rdc
dep_xC_yC = zeros(1,5);
dep_xC_yD = zeros(1,5);
dep_xD_yC = zeros(1,5);
dep_xD_yD = zeros(1,5);

tau = corr(x,y,'type','kendall');
taub = ktaub([x y], alpha_dontCare, wantplot_dontCare);
tau_hat = ktauhat(x, y);
minestats = mine(x',y',mine_alpha,mine_c,'mic_e');
rdcVal = rdc(x,y,rdc_k,rdc_s);
dep_xC_yC(1:5) = [tau taub tau_hat minestats.mic rdcVal];

tau = corr(x,yy,'type','kendall');
taub = ktaub([x yy], alpha_dontCare, wantplot_dontCare);
tau_hat = ktauhat(x, yy);
minestats = mine(x',yy',mine_alpha,mine_c,'mic_e');
rdcVal = rdc(x,yy,rdc_k,rdc_s);
dep_xC_yD(1:5) = [tau taub tau_hat minestats.mic rdcVal];

tau = corr(xx,y,'type','kendall');
taub = ktaub([xx y], alpha_dontCare, wantplot_dontCare);
tau_hat = ktauhat(xx, y);
minestats = mine(xx',y',mine_alpha,mine_c,'mic_e');
rdcVal = rdc(xx,y,rdc_k,rdc_s);
dep_xD_yC(1:5) = [tau taub tau_hat minestats.mic rdcVal];

tau = corr(xx,yyy,'type','kendall');
taub = ktaub([xx yyy], alpha_dontCare, wantplot_dontCare);
tau_hat = ktauhat(xx, yyy);
minestats = mine(xx',yyy',mine_alpha,mine_c,'mic_e');
rdcVal = rdc(xx,yyy,rdc_k,rdc_s);
dep_xD_yD(1:5) = [tau taub tau_hat minestats.mic rdcVal];

figure;
subplot(2,2,1);
scatter(x,y); grid on;
title(sprintf('$\\tau=%0.02f \\ \\tau_b=%0.02f \\ \\hat{\\tau}=%0.02f \\ MIC_e=%0.02f \\ RDC=%0.02f$', ...
    dep_xC_yC(1), dep_xC_yC(2), dep_xC_yC(3), dep_xC_yC(4), dep_xC_yC(5) ), 'Interpreter', 'Latex');

subplot(2,2,2);
scatter(x,yy); grid on;
title(sprintf('$\\tau=%0.02f \\ \\tau_b=%0.02f \\ \\hat{\\tau}=%0.02f \\ MIC_e=%0.02f \\ RDC=%0.02f$', ...
    dep_xC_yD(1), dep_xC_yD(2), dep_xC_yD(3), dep_xC_yD(4), dep_xC_yD(5) ), 'Interpreter', 'Latex');

subplot(2,2,3);
scatter(xx,y); grid on;
title(sprintf('$\\tau=%0.02f \\ \\tau_b=%0.02f \\ \\hat{\\tau}=%0.02f \\ MIC_e=%0.02f \\ RDC=%0.02f$', ...
    dep_xD_yC(1), dep_xD_yC(2), dep_xD_yC(3), dep_xD_yC(4), dep_xD_yC(5) ), 'Interpreter', 'Latex');

subplot(2,2,4);
scatter(xx,yyy); grid on;
title(sprintf('$\\tau=%0.02f \\ \\tau_b=%0.02f \\ \\hat{\\tau}=%0.02f \\ MIC_e=%0.02f \\ RDC=%0.02f$', ...
    dep_xD_yD(1), dep_xD_yD(2), dep_xD_yD(3), dep_xD_yD(4), dep_xD_yD(5) ), 'Interpreter', 'Latex');


%% Understand how RSDM works w/ discrete function dependencies (TODO)

%% Understand how scanfordep handles discrete and hybrid data (TODO)