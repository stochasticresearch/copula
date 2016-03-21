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
%**************************************************************************

% Script which compares the distribution of numbers generated from
% rF01Joe in R and the Matlab port

clear;
clc;

% Test 1 -> theta0 = 0.4, theta1 = 3, rej = 10, approx = 100
V0 = csvread('../../testfiles/rF01Frank_input1.csv');
y_R = csvread('../../testfiles/rF01Frank_output1.csv');
y_Matlab = F01Frankrnd(V0, 0.4, 3, 10, 100);
figure; qqplot(y_R, y_Matlab); title('Test 1')

% Test 2 -> theta0 = 0.8, theta1 = 8, rej = 50, approx = 200
V0 = csvread('../../testfiles/rF01Frank_input2.csv');
y_R = csvread('../../testfiles/rF01Frank_output2.csv');
y_Matlab = F01Frankrnd(V0, 0.8, 8, 50, 200);
figure; qqplot(y_R, y_Matlab); title('Test 2')
