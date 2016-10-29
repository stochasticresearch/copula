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

M = size(X,1);
N = size(R_rsdm,2);

% compute the structure using an aracne style algorithm, but instead of the
% data processing inequality, we use conditional (in)dependence as a
% surrogate
alpha = 0.95;       % significance level for independence testing
edges_rsdm = R_rsdm;
for ii=1:N
    for jj=1:N
        rsdm_val = R_rsdm(ii,jj);
%         pval = rsdm_pval(rsdm_val, M);
%         thresh = alpa*3;        % TODO; query the GEV distribution for thresh
%         if(pval > thres)
%             rsdm_indep = 1;
%         else
%             rsdm_indep = 0;
%         end
        
        % TODO: a p-value based test!
        if(rsdm_val < 0.1)
            rsdm_indep = 1;
            edges_rsdm(ii,jj) = 0; edges_rsdm(jj,ii) = 0;
        else
            rsdm_indep = 0;
        end
        if(~rsdm_indep)
            for kk=1:N
                if(kk~=ii && kk~=jj)
                    rscdm_val = rscdm(X(:,ii), X(:,jj), X(:,kk));
                    % test for the conditional independence scenario
                    if(rsdm_val > rscdm_val)    % TODO: need a p-value based test
                        edges_rsdm(ii,jj) = 0; edges_rsdm(jj,ii) = 0;
                    end
                end
            end
        end
    end
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
