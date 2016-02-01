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

% Compares the distribution of observations with known marginal
% distributions versus pseudo-observations computed with empirical marginal
% distributions

clear;
clc;

% TODO: maybe monte-carlo?

% generate dependent random variables w/ gaussian copula
M = 1000;
alpha = 5;
u = copularnd('Frank', alpha, M);

mu = 2; a = 3; b = 2;
x = [expinv(u(:,1),mu), gaminv(u(:,2), a, b)];

numECDFPts = 100;
U_pseudoObs = pseudoobs(x, 'ecdf', numECDFPts);
U_rank = pseudoobs(x);

subplot(3,3,1); scatter(x(:,1),x(:,2)); grid on; xlabel('X'); ylabel('Y');
subplot(3,3,2); scatter(U_pseudoObs(:,1),U_pseudoObs(:,2)); grid on; xlabel('U_{pseudo}'); ylabel('V_{pseudo}');
subplot(3,3,3); scatter(U_rank(:,1),U_rank(:,2)); grid on; xlabel('U_{rank}'); ylabel('V_{rank}');

nHistBins = 25;
subplot(3,3,4); histogram(x(:,1),nHistBins); grid on; title('f_X(x)')
subplot(3,3,5); histogram(U_pseudoObs(:,1),nHistBins); grid on; title('f_{U_{pseudo}}(u)')
subplot(3,3,6); histogram(U_rank(:,1),nHistBins); grid on; title('f_{U_{rank}}(u)')

subplot(3,3,7); histogram(x(:,2),nHistBins); grid on; title('f_Y(y)')
subplot(3,3,8); histogram(U_pseudoObs(:,2),nHistBins); grid on; title('f_{V_{pseudo}}(v)')
subplot(3,3,9); histogram(U_rank(:,2),nHistBins); grid on; title('f_{V_{rank}}(v)')