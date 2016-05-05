%**************************************************************************
%* 
%* Copyright (C) 2016  Kiran Karra <kiran.karra@gmail.com>
%*
%* This program is free software: you can redistribute it and/or modify
%* it under the terms of the GNU General Public License as published by
%* the Free Software Foundation, either version 3 of the License, or
%* (at your option) any later version.
%*
%* This program is distributed in the hope that it will be useful,
%* but WITHOUT ANY WARRANTY; without even the implied warranty of
%* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%* GNU General Public License for more details.
%*
%* You should have received a copy of the GNU General Public License
%* along with this program.  If not, see <http://www.gnu.org/licenses/>
%* 
%**************************************************************************

% A script to test the model selection algorithm

%% 2-D tests
clear;
clc;

Rho = [1 0.4; 0.4 1];
alpha = 3;
M = 1000;

U_Gaussian = copularnd('Gaussian', Rho, M);
U_Frank = frankcopularnd(M, 2, alpha);
U_Gumbel = gumbelcopularnd(M, 2, alpha);
U_Clayton = claytoncopularnd(M, 2, alpha);

fprintf('----------------------------------------------------\n');
[modelType, modelParams] = copulamodelselect(U_Gaussian)
fprintf('----------------------------------------------------\n');
[modelType, modelParams] = copulamodelselect(U_Frank)
fprintf('----------------------------------------------------\n');
[modelType, modelParams] = copulamodelselect(U_Gumbel)
fprintf('----------------------------------------------------\n');
[modelType, modelParams] = copulamodelselect(U_Clayton)
fprintf('----------------------------------------------------\n');

%% 3-D tests