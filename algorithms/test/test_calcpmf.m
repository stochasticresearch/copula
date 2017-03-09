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
%* along with this program.  If not, see <http://www.gnu.org/licenses/>.
%* 
%**************************************************************************


%% Test 1 -- simple all discrete
clear;
clc;

X = [1 1; ...
     1 1; ...
     2 2; ...
     3 3; ...
     4 3; ...
     1 1; ...
     2 2; ...
     3 3; ...
     4 3; ...
     1 1; ...
     2 2; ...
     3 3; ...
     4 3];
X(:,1) = X(:,1)/4;
X(:,2) = X(:,2)/3;
[h, mapping] = calcpmf( X );
h
celldisp(mapping)

%% Test 1 - functional discrete
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

X = [xx yyy];
[h, mapping] = calcpmf( X );
h
mapping{h~=0}
[tau_plus, tau_minus] = carleybounds(h);
scatter(xx, yyy); grid on;
title(sprintf('\\tau_{-}=%0.02f \\tau_{+}=%0.02f', tau_plus, tau_minus));

%% Test 2 -- discrete and continuous


%% Test 3 -- continuous and discrete


%% Test 4 -- continuous and continuous
