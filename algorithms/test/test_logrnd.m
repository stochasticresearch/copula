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
% rlog in R and the Matlab port

clear;
clc;

% Test 1 -> p = 0.2
y_R = csvread('../../testfiles/rlog_output1.csv');
y_Matlab = logrnd(length(y_R), 0.2);
subplot(1,2,1); hist(y_R);
subplot(1,2,2); hist(y_Matlab);