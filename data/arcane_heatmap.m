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
%%  process a heatmap of the test ARCANE gene expression data
clear;
clc;

% load the arcane data --
% Each row contains a microarray experiment and each column contains a gene
X = csvread('arcane.data',1,1);
R_rsdm = pairrsdm(X);
R_rdc = pairrdc(X);
R_mice = pairmice(X);

% load the results from data for ARACNE algorithm
arcaneResults = csvread('arcane.res', 1, 1);

% save the data
if(ispc)
    save('C:\\Users\\Kiran\\ownCloud\\PhD\\sim_results\\independence\\arcaneRSDM.mat');
elseif(ismac)
    save('/Users/Kiran/ownCloud/PhD/sim_results/independence/arcaneRSDM.mat');
else
    save('/home/kiran/ownCloud/PhD/sim_results/independence/arcaneRSDM.mat');
end

%% Load the data for plotting
if(ispc)
    load('C:\\Users\\Kiran\\ownCloud\\PhD\\sim_results\\independence\\arcaneRSDM.mat');
elseif(ismac)
    load('/Users/Kiran/ownCloud/PhD/sim_results/independence/arcaneRSDM.mat');
else
    load('/home/kiran/ownCloud/PhD/sim_results/independence/arcaneRSDM.mat');
end

I = (arcaneResults==0);
% R_rsdm(I) = 0;
% R_rdc(I) = 0;
% R_mice(I) = 0;
R_rsdm(R_rsdm==0) = eps;
R_rdc(R_rdc==0) = eps;
R_mice(R_mice==0) = eps;

colormap('jet');

subplot(2,2,1);
im = imagesc(arcaneResults);        % draw image and scale colormap to values range
im.AlphaData = .8;
colorbar;          % show color scale
title('ARACNe');

subplot(2,2,2);
im = imagesc(R_rsdm./R_rdc);        % draw image and scale colormap to values range
im.AlphaData = .8;
colorbar;          % show color scale
title('RSDM/RDC');

subplot(2,2,3);
im = imagesc(R_rsdm./R_mice);        % draw image and scale colormap to values range
im.AlphaData = .8;
colorbar;          % show color scale
title('RSDM/MICe');

subplot(2,2,4);
im = imagesc(R_rdc./R_mice);        % draw image and scale colormap to values range
im.AlphaData = .8;
colorbar;          % show color scale
title('RDC/MICe');
