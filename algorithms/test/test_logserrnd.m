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
% rlogseries in R and the Matlab port


clear;
clc;
close all;

% Test 1 -> alpha = 0.75
y_R = csvread('../../testfiles/rlogseries_output1.csv');
y_Matlab = logserrnd(length(y_R), 0.75);
figure;
qqplot(log(y_R), log(y_Matlab)); title('Test 1 - \alpha = 0.75'); grid on

% Test 1 -> alpha = 0.25
y_R = csvread('../../testfiles/rlogseries_output2.csv');
y_Matlab = logserrnd(length(y_R), 0.25);
figure;
qqplot(log(y_R), log(y_Matlab)); title('Test 2 - \alpha=0.25'); grid on

% Test 1 -> alpha = 0.99
y_R = csvread('../../testfiles/rlogseries_output3.csv');
y_Matlab = logserrnd(length(y_R), 0.99);
figure;
qqplot(log(y_R), log(y_Matlab)); title('Test 3 - \alpha=0.99'); grid on