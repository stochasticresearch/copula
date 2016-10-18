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
y = x + 1*randn(M,1);

% x = rand(M,1);
% y = sin(4*pi*x);

% x = rand(M,1)*2-1;
% y = x.^2;

% x = rand(M,1);
% y=(2*binornd(1,0.5,M,1)-1).* (sqrt(1 - (2*x - 1).^2));

[metric, resid, residAssocIdxs, residAssocPts, rectangleCfg] = rsdm_2(x,y);

%% Understand how RSDM works w/ discrete function dependencies (TODO)

%% Understand how scanfordep handles discrete and hybrid data (TODO)