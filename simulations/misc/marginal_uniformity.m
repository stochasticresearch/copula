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
U_ecdf = pseudoobs(x, 'ecdf', numECDFPts);
U_rank = pseudoobs(x);

figure;
subplot(3,3,1); scatter(x(:,1),x(:,2)); grid on; xlabel('X'); ylabel('Y');
subplot(3,3,2); scatter(U_ecdf(:,1),U_ecdf(:,2)); grid on; xlabel('U_{pseudo}'); ylabel('V_{pseudo}');
subplot(3,3,3); scatter(U_rank(:,1),U_rank(:,2)); grid on; xlabel('U_{rank}'); ylabel('V_{rank}');

nHistBins = 25;
subplot(3,3,4); histogram(x(:,1),nHistBins); grid on; title('f_X(x)')
subplot(3,3,5); histogram(U_ecdf(:,1),nHistBins); grid on; title('f_{U_{pseudo}}(u)')
subplot(3,3,6); histogram(U_rank(:,1),nHistBins); grid on; title('f_{U_{rank}}(u)')

subplot(3,3,7); histogram(x(:,2),nHistBins); grid on; title('f_Y(y)')
subplot(3,3,8); histogram(U_ecdf(:,2),nHistBins); grid on; title('f_{V_{pseudo}}(v)')
subplot(3,3,9); histogram(U_rank(:,2),nHistBins); grid on; title('f_{V_{rank}}(v)')

% Show the impact of using pseudo-observations versus  assuming marginal 
% distributions are known for a given point in the hypercube
u_v = [0.2 0.3];      % known point in hypercube
true_value = copulapdf('Frank', u_v, alpha);
nsim = 200;
h = 0.05;
K = 50;
c_ecdf = zeros(1,nsim); c_rank = zeros(1,nsim);
for ii=1:nsim
    
    fprintf('Sim %d\n', ii);
    
    % generate the data
    M = 1000;
    alpha = 5;
    u = copularnd('Frank', alpha, M);

    mu = 2; a = 3; b = 2;
    x = [expinv(u(:,1),mu), gaminv(u(:,2), a, b)];

    numECDFPts = 100;
    U_ecdf = pseudoobs(x, 'ecdf', numECDFPts);
    U_rank = pseudoobs(x);
    
    c_ecdf(ii) = empcopula_val(empcopulapdf(U_ecdf, h, K, 'betak'), u_v);
    c_rank(ii) = empcopula_val(empcopulapdf(U_rank, h, K, 'betak'), u_v);
end

figure;
numKsDensityPts = 50; pts = linspace( min( min(c_ecdf), min(c_rank) ), max( max(c_ecdf), max(c_rank) ), numKsDensityPts );
[f_ecdf,xi] = ksdensity(c_ecdf, pts);
[f_rank] = ksdensity(c_rank, pts);
plot(pts, f_ecdf, 'b--');  hold on; grid on;
plot(pts, f_rank, 'r*');
stem(true_value, max( max(f_ecdf), max(f_rank) ) + 1);
var_ecdf = var(c_ecdf); var_rank = var(c_rank);
legend( sprintf('Known Marginals Estimate VAR=%.4f', var_ecdf), ...
        sprintf('Pseudo-Samples Estimate VAR=%.4f', var_rank), ...
        'True Value' )
title('Performance of Estimators of C, Frank Copula, \alpha=5, [u,v]=[0.2,0.3]')