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
%
% simulation to show how to generate linked gaussian rv's.  In the example
% below, U_3d is a Gaussian copula w/ a certain correlation matrix.  We'd
% like to generate U_2d, which is linked to U_3d in the sense that it
% shares one of the variables, but the 2nd variable in U_2d should be
% generated in such a way that the relationship between the linked variable
% from U_3d and the 2nd variable in U_2d have the correlation matrix given
% in Rho2.  i.e. this is a BN w/ the following structure:
%      A    B
%     /  \ /
%    C    D
% where all the arrows point downward.  Rho1 represents the dependency
% structure between A,B,D and Rho2 represents the dependency structure
% between A,C.  Note that we use 'upper' for the cholesky decomposition,
% not 'lower'.  Lower is not correct in this situation ... as the testing
% below shows

clear;
clc;

Rho1 = [1 0 .2; 0 1 -.8; .2 -.8 1];
Rho2 = [1 -0.6; -0.6 1];

U_3d = copularnd('Gaussian', Rho1, 1000);
copulafit('Gaussian', U_3d)

u1 = U_3d(:,1);
x1 = norminv(u1, 0, 1);

fprintf('Using upper option -- Rho2 fitted matrix --->\n');
U = chol(Rho2,'upper');
x2 = [x1 normrnd(0, 1, length(u1), 1)]*U;
U_2d_2 = normcdf(x2(:,2));
U_2d = [u1 U_2d_2];
copulafit('Gaussian', U_2d)

fprintf('Using lower option (notice this is incorrect!) -- Rho2 fitted matrix --->\n');
L = chol(Rho2,'lower');
x2 = [x1 normrnd(0, 1, length(u1), 1)]*L;
U_2d_2 = normcdf(x2(:,2));
U_2d = [u1 U_2d_2];
copulafit('Gaussian', U_2d)

