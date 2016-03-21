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

% Test 1 -> alpha = 0.3
V0 = csvread('../../testfiles/rF01Joe_input1.csv');
y_R = csvread('../../testfiles/rF01Joe_output1.csv');
y_Matlab = F01Joernd(V0, 0.3, 100000);
figure; qqplot(y_R, y_Matlab); title('Test 1')

% Test 2 --> alpha = 0.8
V0 = csvread('../../testfiles/rF01Joe_input2.csv');
y_R = csvread('../../testfiles/rF01Joe_output2.csv');
y_Matlab = F01Joernd(V0, 0.8, 100000);
figure; qqplot(y_R, y_Matlab); title('Test 2')
