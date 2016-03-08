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

% a test script for testing the multivariate (N>2) versions of generating
% random variates from the Frank, Gumbel, and Clayton copulas

%% Test N = 2
M = 1000;
N = 2;

alpha = 5;

U_clayton = claytoncopularnd(M,N,alpha);
U_frank   = frankcopularnd(M,N,alpha);
U_gumbel  = gumbelcopularnd(M,N,alpha);

U_clayton_matlab = copularnd('Clayton', alpha, M);
U_frank_matlab = copularnd('Frank', alpha, M);
U_gumbel_matlab = copularnd('Gumbel', alpha, M);

figure;
subplot(1,2,1); scatter(U_clayton(:,1), U_clayton(:,2)); grid on; title('Clayton Copula')
subplot(1,2,2); scatter(U_clayton_matlab(:,1), U_clayton_matlab(:,2)); grid on; title('Clayton Copula - Matlab')

figure;
subplot(1,2,1); scatter(U_frank(:,1), U_frank(:,2)); grid on; title('Frank Copula')
subplot(1,2,2); scatter(U_frank_matlab(:,1), U_frank_matlab(:,2)); grid on; title('Frank Copula - Matlab')

figure;
subplot(1,2,1); scatter(U_gumbel(:,1), U_gumbel(:,2)); grid on; title('Gumbel Copula')
subplot(1,2,2); scatter(U_gumbel_matlab(:,1), U_gumbel_matlab(:,2)); grid on; title('Gumbel Copula - Matlab')


%% Test N = 3
M = 1000;
N = 3;

alpha = 5;

U_clayton = claytoncopularnd(M,N,alpha);
U_frank   = frankcopularnd(M,N,alpha);
U_gumbel  = gumbelcopularnd(M,N,alpha);

figure;
plotmatrix(U_clayton);
title('Multivariate Clayton Copula')

figure;
plotmatrix(U_frank);
title('Multivariate Frank Copula')

figure;
plotmatrix(U_gumbel);
title('Multivariate Gumbel Copula')
