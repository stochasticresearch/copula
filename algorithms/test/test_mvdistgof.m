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

%% TEST 1
% Test the eqdistetest w/ R implementation
clear;
clc;

u = csvread('../energy-r/testfiles/eqdisttest1_input1.csv');
v = csvread('../energy-r/testfiles/eqdisttest1_input2.csv');

nperm = 200;

X = [u;v];
sizes = [size(u,1) size(v,1)];
p1 = mvdistgof(X, sizes, nperm);
fprintf('p = %f -- should be close to 0\n', p1);

X = [u;u];
sizes = [size(u,1) size(u,1)];
p2 = mvdistgof(X, sizes, nperm);
fprintf('p = %f -- should be close to 1\n', p2);

X = [v;v];
sizes = [size(v,1) size(v,1)];
p3 = mvdistgof(X, sizes, nperm);
fprintf('p = %f -- should be close to 1\n', p3);

%% TEST 2

u = csvread('../energy-r/testfiles/eqdisttest2_input1.csv');
v = csvread('../energy-r/testfiles/eqdisttest2_input2.csv');

nperm = 200;

X = [u;v];
sizes = [size(u,1) size(v,1)];
p1 = mvdistgof(X, sizes, nperm);
fprintf('p = %f -- should be close to 0\n', p1);

X = [u;u];
sizes = [size(u,1) size(u,1)];
p2 = mvdistgof(X, sizes, nperm);
fprintf('p = %f -- should be close to 1\n', p2);

X = [v;v];
sizes = [size(v,1) size(v,1)];
p3 = mvdistgof(X, sizes, nperm);
fprintf('p = %f -- should be close to 1\n', p3);