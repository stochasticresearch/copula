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
% rSibuya in R and the Matlab port

% we do Q-Q plots of the log of the random variates because values are
% bunched up in smaller values, and log smoothes this out. Is this valid?

clear;
clc;
close all;

% Test 1 -> alpha = 0.3
y_R = csvread('../../testfiles/rSibuya_output1.csv');
y_Matlab = sibuyarnd(length(y_R), 0.3, gamma(0.7));
figure;
qqplot(log(y_R), log(y_Matlab));

% Test 1 -> alpha = 0.9
y_R = csvread('../../testfiles/rSibuya_output2.csv');
y_Matlab = sibuyarnd(length(y_R), 0.9, gamma(0.1));
figure;
qqplot(log(y_R), log(y_Matlab));
