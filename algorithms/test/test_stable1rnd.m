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
% rstable1 and stable1rnd in R and the Matlab port 

clear;
clc;

% Test1, alpha = 0.8, beta = 1, gamma=1, delta=0
y_R = csvread('../../testfiles/rstable1_output1.csv');
y_Matlab = stable1rnd(1000, 0.8, 1, 0);
figure; qqplot(log(y_R), log(y_Matlab)); title('Test 1')
 
% Test2, alpha = 0.25, beta = 1, gamma=1, delta=0
y_R = csvread('../../testfiles/rstable1_output2.csv');
y_Matlab = stable1rnd(1000, 0.25, 1, 0);
figure; qqplot(log(y_R), log(y_Matlab)); title('Test 2')