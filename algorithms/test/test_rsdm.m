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

M = 500;

numDiscreteIntervals = 10;

% Strictly monotonic
x = rand(M,1);
y = x.^3;
xx = discretizeRv(x,numDiscreteIntervals)';
yy = discretizeRv(y,numDiscreteIntervals)';
yyy = xx.^3;        % both x and y are discrete here

% some required but not used configuration parameters
alpha_dontCare = 0.05;
wantplot_dontCare = 0;

% calculate the different kind's of tau for each of these
% continuous/discrete/hybrid scenarios
% we have 3 measures -- tau, tau_b, tau_cj
taus_xC_yC = zeros(1,3);
taus_xC_yD = zeros(1,3);
taus_xD_yC = zeros(1,3);
taus_xD_yD = zeros(1,3);

[~, tau] = ktaub([x y], alpha_dontCare, wantplot_dontCare);
taus_xC_yC(1) = tau;

[taub, tau] = ktaub([x yy], alpha_dontCare, wantplot_dontCare);
taucj = ktaucj(x, yy);
taus_xC_yD(1:3) = [tau taub taucj];

[taub, tau] = ktaub([xx y], alpha_dontCare, wantplot_dontCare);
taucj = ktaucj(xx, y);
taus_xD_yC(1:3) = [tau taub taucj];

[taub, tau] = ktaub([xx yyy], alpha_dontCare, wantplot_dontCare);
taucj = ktaucj(xx, yyy);
taus_xD_yD(1:3) = [tau taub taucj];

figure;
subplot(2,2,1);
scatter(x,y); grid on;
title(sprintf('\\tau=%0.02f', taus_xC_yC(1)) );

subplot(2,2,2);
scatter(x,yy); grid on;
title(sprintf('\\tau=%0.02f \\tau_b=%0.02f \\tau_{cj}=%0.02f', ...
    taus_xC_yD(1), taus_xC_yD(2), taus_xC_yD(3) ));

subplot(2,2,3);
scatter(xx,y); grid on;
title(sprintf('\\tau=%0.02f \\tau_b=%0.02f \\tau_{cj}=%0.02f', ...
    taus_xD_yC(1), taus_xD_yC(2), taus_xD_yC(3) ));

subplot(2,2,4);
scatter(xx,yyy); grid on;
title(sprintf('\\tau=%0.02f \\tau_b=%0.02f \\tau_{cj}=%0.02f', ...
    taus_xD_yD(1), taus_xD_yD(2), taus_xD_yD(3) ));

% Strictly counter-monotonic
% Strictly monotonic
x = rand(M,1);
y = -x.^3;
xx = discretizeRv(x,numDiscreteIntervals)';
yy = discretizeRv(y,numDiscreteIntervals)';
yyy = -xx.^3;        % both x and y are discrete here

% some required but not used configuration parameters
alpha_dontCare = 0.05;
wantplot_dontCare = 0;

% calculate the different kind's of tau for each of these
% continuous/discrete/hybrid scenarios
% we have 3 measures -- tau, tau_b, tau_cj
taus_xC_yC = zeros(1,3);
taus_xC_yD = zeros(1,3);
taus_xD_yC = zeros(1,3);
taus_xD_yD = zeros(1,3);

[~, tau] = ktaub([x y], alpha_dontCare, wantplot_dontCare);
taus_xC_yC(1) = tau;

[taub, tau] = ktaub([x yy], alpha_dontCare, wantplot_dontCare);
taucj = ktaucj(x, yy);
taus_xC_yD(1:3) = [tau taub taucj];

[taub, tau] = ktaub([xx y], alpha_dontCare, wantplot_dontCare);
taucj = ktaucj(xx, y);
taus_xD_yC(1:3) = [tau taub taucj];

[taub, tau] = ktaub([xx yyy], alpha_dontCare, wantplot_dontCare);
taucj = ktaucj(xx, yyy);
taus_xD_yD(1:3) = [tau taub taucj];

figure;
subplot(2,2,1);
scatter(x,y); grid on;
title(sprintf('\\tau=%0.02f', taus_xC_yC(1)) );

subplot(2,2,2);
scatter(x,yy); grid on;
title(sprintf('\\tau=%0.02f \\tau_b=%0.02f \\tau_{cj}=%0.02f', ...
    taus_xC_yD(1), taus_xC_yD(2), taus_xC_yD(3) ));

subplot(2,2,3);
scatter(xx,y); grid on;
title(sprintf('\\tau=%0.02f \\tau_b=%0.02f \\tau_{cj}=%0.02f', ...
    taus_xD_yC(1), taus_xD_yC(2), taus_xD_yC(3) ));

subplot(2,2,4);
scatter(xx,yyy); grid on;
title(sprintf('\\tau=%0.02f \\tau_b=%0.02f \\tau_{cj}=%0.02f', ...
    taus_xD_yD(1), taus_xD_yD(2), taus_xD_yD(3) ));


%% Understand how RSDM works w/ discrete function dependencies (TODO)

%% Understand how scanfordep handles discrete and hybrid data (TODO)