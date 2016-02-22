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

% Script which characterizes the error between the multivariate Archimedean
% copula PDF calculation and the results from R

clear;
clc;

alpha = 4;

% read input test file
u = csvread('../../testfiles/mvArchimedeanCopula_input.csv');
y_claytonR = csvread('../../testfiles/claytonPdf3D_output.csv');
y_frankR = csvread('../../testfiles/frankPdf3D_output.csv');
y_gumbelR = csvread('../../testfiles/gumbelPdf3D_output.csv');

y_claytonMlab = claytoncopulapdf(u, alpha);
y_frankMlab = frankcopulapdf(u, alpha);
y_gumbelMlab = gumbelcopulapdf(u, alpha);


mse_clayton_3D = sum((y_claytonR-y_claytonMlab).^2);
mse_frank_3D = sum((y_frankR-y_frankMlab).^2);
mse_gumbel_3D = sum((y_gumbelR-y_gumbelMlab).^2);

mse_clayton_3D_str = strrep(sprintf('%8.2E',mse_clayton_3D),'E-0','E-');
mse_frank_3D_str = strrep(sprintf('%8.2E',mse_frank_3D),'E-0','E-');
mse_gumbel_3D_str = strrep(sprintf('%8.2E',mse_gumbel_3D),'E-0','E-');

% plot jsut for sanity checks
subplot(1,3,1);
plot(1:length(y_claytonR), y_claytonR, 'b*', 1:length(y_claytonR), y_claytonMlab, 'r-');
grid on; title(sprintf('Clayton Copula | MSE=%s', mse_clayton_3D_str)); legend('R', 'Matlab'); ylabel('MSE');

subplot(1,3,2);
plot(1:length(y_frankR), y_frankR, 'b*', 1:length(y_frankR), y_frankMlab, 'r-');
grid on; title(sprintf('Frank Copula | MSE=%s', mse_frank_3D_str)); legend('R', 'Matlab'); ylabel('MSE');

subplot(1,3,3);
plot(1:length(y_gumbelR), y_gumbelR, 'b*', 1:length(y_gumbelR), y_gumbelMlab, 'r-');
grid on; title(sprintf('Gumbel Copula | MSE=%s', mse_gumbel_3D_str)); legend('R', 'Matlab'); ylabel('MSE');
