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

% Test the mvdistgof against various examples
%% compare MVNRND against MVTRND
clear;
clc;

SIGMA = [1 0.8;0.8 1];
X = mvnrnd([0 0], SIGMA, 100);
Y = mvtrnd(SIGMA, 3, 100);

alpha = 0.05;
nperm = 100;
[h, p] = mvdistgof(X, Y, alpha, nperm)
