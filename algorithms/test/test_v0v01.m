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
% V0 and V01 in R and the Matlab port for the three Archimedean copulas

clear;
clc;

% Clayton copula, theta0 = 0.5, theta1 = 1.7
V0_R = csvread('../../testfiles/clayton_v0.csv');
V01_R = csvread('../../testfiles/clayton_v01.csv');
V0_matlab = claytonV0rnd(1000, 0.5);
V01_matlab = claytonV01rnd(V0_R, 0.5, 1.7);
figure; 
subplot(1,2,1); qqplot(V0_R, V0_matlab); title('V0 - Clayton')
subplot(1,2,2); qqplot(V01_R, V01_matlab); title('V01 - Clayton')

% Gumbel copula, theta0 = 1.1, theta1 = 2.55
V0_R = csvread('../../testfiles/gumbel_v0.csv');
V01_R = csvread('../../testfiles/gumbel_v01.csv');
V0_matlab = gumbelV0rnd(1000, 1.1);
V01_matlab = gumbelV01rnd(V0_R, 1.1, 2.55);
figure; 
subplot(1,2,1); qqplot(V0_R, V0_matlab); title('V0 - Gumbel')
subplot(1,2,2); qqplot(V01_R, V01_matlab); title('V01 - Gumbel')

% Frank copula, theta0 = 0.5, theta1 = 5.0
V0_R = csvread('../../testfiles/frank_v0.cfsv');
V01_R = csvread('../../testfiles/frank_v01.csv');
V0_matlab = frankV0rnd(1000, 1.1);
V01_matlab = frankV01rnd(V0_R, 1.1, 2.55);
figure; 
subplot(1,2,1); qqplot(V0_R, V0_matlab); title('V0 - Frank')
subplot(1,2,2); qqplot(V01_R, V01_matlab); title('V01 - Frank')

