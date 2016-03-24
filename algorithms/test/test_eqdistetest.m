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

% Test the eqdistetest w/ R implementation
clear;
clc;

u = csvread('/home/kiran/stochasticresearch/energy-r/testfiles/test1_input1.csv');
v = csvread('/home/kiran/stochasticresearch/energy-r/testfiles/test1_input2.csv');

nperm = 200;
p = eqdistetest(u,v,nperm)