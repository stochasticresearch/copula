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
%%  plot a heatmap of the test ARCANE gene expression data
clear;
clc;

% load the arcane data --
% Each row contains a microarray experiment and each column contains a gene
X = csvread('arcane.data',1,1);
R = pairrsdm(X);

% save the data
if(ispc)
    save('C:\\Users\\Kiran\\ownCloud\\PhD\\sim_results\\independence\\arcaneRSDM.mat');
elseif(ismac)
    save('/Users/Kiran/ownCloud/PhD/sim_results/independence/arcaneRSDM.mat');
else
    save('/home/kiran/ownCloud/PhD/sim_results/independence/arcaneRSDM.mat');
end

% load the results from data for ARACNE algorithm
arcaneResults = csvread('arcane.res', 1, 1);

subplot(1,2,1);
myheatmap(R)
title('RSDM');
subplot(1,2,2);
myheatmap(arcaneResults);