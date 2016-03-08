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

% Tests the depcopularnd script

%% Test Gaussian copula CI test with partial correlation
clear;
clc;

Rho1 = [1 0.5; 0.5 1];
Rho2 = [1 -0.3; -0.3 1];
M = 1000;

U_init = copularnd('Gaussian', Rho1, M);        % generates [Z X]
U_dep = depcopularnd(U_init(:,1), 2, 'Gaussian', Rho2); % input [Z] to generate [Z Y]

Z1 = U_init(:,1); Z2 = U_dep(:,1);
X = U_init(:,2);
Y = U_dep(:,2);

fprintf('********** Gaussian RV Testing - X indep Y | Z **********\n');
fprintf('rho(X,Y|Z1) = %f\n', partialcorr(X,Y,Z1));
fprintf('rho(X,Y) = %f\n', corr(X,Y));
fprintf('rho(X,Z1) = %f\n', corr(X,Z1));
fprintf('rho(Y,Z1) = %f\n', corr(Y,Z1));

fprintf('rho(X,Y|Z2) = %f\n', partialcorr(X,Y,Z2));
fprintf('rho(X,Z2) = %f\n', corr(X,Z2));
fprintf('rho(Y,Z2) = %f\n', corr(Y,Z2));

fprintf('\nrho(Z1,Z2)=%f\n', corr(Z1,Z2));

Rho3 = [1 0.5 0.3; 0.5 1 -0.2; 0.3 -0.2 1];
U = mvnrnd([0 0 0], Rho3, M);
X = U(:,1); Y = U(:,2); Z1 = U(:,3);

fprintf('********** Gaussian RV Testing - X !indep Y | Z **********\n');
fprintf('rho(X,Y|Z) = %f\n', partialcorr(X,Y,Z1));
fprintf('rho(X,Y) = %f\n', corr(X,Y));
fprintf('rho(X,Z) = %f\n', corr(X,Z1));
fprintf('rho(Y,Z) = %f\n', corr(Y,Z1));

%% Test the Clayton Copula
clear;
clc;

alpha = 5;
M = 1000;
N = 2;
[U_init] = claytoncopularnd(M, N, alpha);
U_dep = depcopularnd(U_init(:,1), N, 'Clayton', alpha+3);

Z1 = U_dep(:,1); Z2 = U_init(:,1);
X = U_init(:,2);
Y = U_dep(:,2);

fprintf('********** CLAYTON Dependency Testing - X indep Y | Z ******\n');
fprintf('rho_s(Z1,Z2)=%f\n', corr(Z1,Z2,'type','Spearman'));

fprintf('rho_s(X,Y|Z1)=%f\n', partialcorr(X, Y, Z1, 'Type', 'Spearman'));
fprintf('rho_s(X,Y)=%f\n', corr(X, Y, 'Type', 'Spearman'));
fprintf('rho_s(X,Z1)=%f\n', corr(X, Z1, 'Type', 'Spearman'));
fprintf('rho_s(Y,Z1)=%f\n', corr(Y, Z1, 'Type', 'Spearman'));

fprintf('rho_s(X,Y|Z2)=%f\n', partialcorr(X, Y, Z2, 'Type', 'Spearman'));
fprintf('rho_s(X,Z2)=%f\n', corr(X, Z2, 'Type', 'Spearman'));
fprintf('rho_s(X,Z2)=%f\n', corr(Y, Z2, 'Type', 'Spearman'));

%% Test the Frank Copula
clear;
clc;

alpha = 5;
M = 1000;
N = 2;
[U_init] = frankcopularnd(M, N, alpha);
U_dep = depcopularnd(U_init(:,1), N, 'Frank', alpha-3);

Z1 = U_dep(:,1); Z2 = U_init(:,1);
X = U_init(:,2);
Y = U_dep(:,2);

fprintf('********** FRANK Dependency Testing - X indep Y | Z ******\n');
fprintf('rho_s(Z1,Z2)=%f\n', corr(Z1,Z2,'type','Spearman'));

fprintf('rho_s(X,Y|Z1)=%f\n', partialcorr(X, Y, Z1, 'Type', 'Spearman'));
fprintf('rho_s(X,Y)=%f\n', corr(X, Y, 'Type', 'Spearman'));
fprintf('rho_s(X,Z1)=%f\n', corr(X, Z1, 'Type', 'Spearman'));
fprintf('rho_s(Y,Z1)=%f\n', corr(Y, Z1, 'Type', 'Spearman'));

fprintf('rho_s(X,Y|Z2)=%f\n', partialcorr(X, Y, Z2, 'Type', 'Spearman'));
fprintf('rho_s(X,Z2)=%f\n', corr(X, Z2, 'Type', 'Spearman'));
fprintf('rho_s(X,Z2)=%f\n', corr(Y, Z2, 'Type', 'Spearman'));