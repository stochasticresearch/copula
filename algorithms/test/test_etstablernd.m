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
% retstable in R and the Matlab port

clear;
clc;

% Test 1 -> alpha = 0.5, h = 1
V0 = csvread('../../testfiles/retstable_input1.csv');
St_R = csvread('../../testfiles/retstable_output1.csv');
St_Matlab = etstablernd(V0, 1, 0.5);
figure; qqplot(St_R, St_Matlab); title('Test 1')

% Test 2 --> alpha = 1, h = 2
V0 = csvread('../../testfiles/retstable_input2.csv');
St_R = csvread('../../testfiles/retstable_output2.csv');
St_Matlab = etstablernd(V0, 2, 1);
figure; qqplot(St_R, St_Matlab); title('Test 2')

% Test 3 --> alpha = 0.1, h = 4
V0 = csvread('../../testfiles/retstable_input3.csv');
St_R = csvread('../../testfiles/retstable_output3.csv');
St_Matlab = etstablernd(V0, 4, 0.1);
figure; qqplot(St_R, St_Matlab); title('Test 3')
