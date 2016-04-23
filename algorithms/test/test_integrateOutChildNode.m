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

% Script to see if we can just integrate out the child node to get the
% copula for the parents with the total copula that we already calculated

clear;
clc;

K = 25;

copulaType = 'Gaussian'; Rho = [1 .4 .2; .4 1 -.8; .2 -.8 1];

u = linspace(0,1,K);
[U1,U2,U3] = ndgrid(u);

c = copulapdf(copulaType, [U1(:) U2(:) U3(:)], Rho);
c = reshape(c,K,K,K);

% integrate out the first dimension
c_integratedout = squeeze(sum(c,3));
[U1,U2] = ndgrid(u);
c_manufactured = copulapdf(copulaType, [U1(:) U2(:)], Rho([1 2],[1 2]));
c_manufactured = reshape(c_manufactured, K, K);

subplot(1,2,1); surf(c_integratedout); grid on;
subplot(1,2,2); surf(c_manufactured); grid on;

% Comments -- it generally seems to work, but the scaling factor is off,
% not sure how to correct the scaling factor for the pdf yet ...