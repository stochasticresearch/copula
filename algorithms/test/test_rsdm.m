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

%% some basic rsdm_s tests before running parametric power tests
clear;
clc;

M = 500;
x = rand(M,1);

y = x; 
lin = rsdm(x,y); 
profile on
lin_s = rsdm_s(x,y);
profile viewer
y = (x-0.5).^2; quad = rsdm(x,y); quad_s = rsdm_s(x,y);
y = 128*(x-1/3).^3-48*(x-1/3).^3-12*(x-1/3); cubi = rsdm(x,y); cubi_s = rsdm_s(x,y);
y = sin(4*pi*x); sine1 = rsdm(x,y); sine1_s = rsdm_s(x,y);
y = sin(16*pi*x); sine2 = rsdm(x,y); sine2_s = rsdm_s(x,y);
y = x.^(1/4); fr = rsdm(x,y); fr_s = rsdm_s(x,y);
y = (2*binornd(1,0.5,M,1)-1) .* (sqrt(1 - (2*x - 1).^2)); cir = rsdm(x,y); cir_s = rsdm_s(x,y);
y = (x > 0.5); step = rsdm(x,y); step_s = rsdm_s(x,y);

fprintf('Linear -- Original=%0.02f, Streaming=%0.02f\n', lin, lin_s);
fprintf('Quadratic -- Original=%0.02f, Streaming=%0.02f\n', quad, quad_s);
fprintf('Cubic -- Original=%0.02f, Streaming=%0.02f\n', cubi, cubi_s);
fprintf('Sine 1 -- Original=%0.02f, Streaming=%0.02f\n', sine1, sine1_s);
fprintf('Sine 2 -- Original=%0.02f, Streaming=%0.02f\n', sine2, sine2_s);
fprintf('Fourth Root -- Original=%0.02f, Streaming=%0.02f\n', fr, fr_s);
fprintf('Circle -- Original=%0.02f, Streaming=%0.02f\n', cir, cir_s);
fprintf('Step -- Original=%0.02f, Streaming=%0.02f\n', step, step_s);

%% Test dependency metric of RSDM-2 for some interesting dependencies...

clear;
clc;

rng(123);

M = 500;
noise = 0.1;
alpha = 0.05;

% Generate data from     Y<--X-->Z
x = rand(M,1);
y = sin(2*pi*x) + noise*randn(M,1); z = 4*(x-0.5).^2 + noise*randn(M,1);

rsdm_val = rsdm_s(y,z)
rdc_val = rdc(y,z)

%% Test residual alignment of RSDM_2
clear;
clc;
close all;

M = 1000;

x = rand(M,1);
y = x;% + 1*randn(M,1);
[metric, resid, residAssocIdxs, residAssocPts, rectangleCfg] = rsdm_2(x,y);
fprintf('linear metric = %0.02f\n', metric);

x = rand(M,1);
y = sin(4*pi*x);
[metric, resid, residAssocIdxs, residAssocPts, rectangleCfg] = rsdm_2(x,y);
fprintf('sin metric = %0.02f\n', metric);

x = rand(M,1)*2-1;
y = x.^2;
[metric, resid, residAssocIdxs, residAssocPts, rectangleCfg] = rsdm_2(x,y);
fprintf('quadratic metric = %0.02f\n', metric);

x = rand(M,1);
y=(2*binornd(1,0.5,M,1)-1).* (sqrt(1 - (2*x - 1).^2));
[metric, resid, residAssocIdxs, residAssocPts, rectangleCfg] = rsdm_2(x,y);
fprintf('circular metric = %0.02f\n', metric);

rng(6);     % expected orientation
x = rand(M,1);
y = (2*binornd(1,0.5,M,1)-1).*sqrt(x);
[metric, resid, residAssocIdxs, residAssocPts, rectangleCfg] = rsdm_2(x,y);
fprintf('half-function metric = %0.02f\n', metric);

rng(7);     % this produces the weird orientation bug :D
x = rand(M,1);
y = (2*binornd(1,0.5,M,1)-1).*sqrt(x);
[metric, resid, residAssocIdxs, residAssocPts, rectangleCfg] = rsdm_2(x,y);
fprintf('half-function metric = %0.02f\n', metric);

x = rand(M,1);
y = 128*(x-1/3).^3-48*(x-1/3).^3-12*(x-1/3);
[metric, resid, residAssocIdxs, residAssocPts, rectangleCfg] = rsdm_2(x,y);
fprintf('cubic metric = %0.02f\n', metric);

x = rand(M,1);
y = x.^(1/4);
[metric, resid, residAssocIdxs, residAssocPts, rectangleCfg] = rsdm_2(x,y);
fprintf('fourth root metric = %0.02f\n', metric);

%% Understand how RSDM works w/ discrete function dependencies (TODO)

%% Understand how scanfordep handles discrete and hybrid data (TODO)