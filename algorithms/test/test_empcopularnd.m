%**************************************************************************
%*                                                                        *
%* Copyright (C) 2016  Kiran Karra <kiran.karra@gmail.com>                *
%*                                                                        *
%* This program is free software: you can redistribute it and/or modify   *
%* it under the terms of the GNU General Public License as published by   *
%* the Free Software Foundation, either version 3 of the License, or      *
%* (at your option) any later version.                                    *
%*                                                                        *
%* This program is distributed in the hope that it will be useful,        *
%* but WITHOUT ANY WARRANTY; without even the implied warranty of         *
%* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          *
%* GNU General Public License for more details.                           *
%*                                                                        *
%* You should have received a copy of the GNU General Public License      *
%* along with this program.  If not, see <http://www.gnu.org/licenses/>.  *
%*                                                                        *
%**************************************************************************

% Test empcopularnd against matlab generation
clear;
clc;
close all;

% Generate the copula pdf
K = 100;
M = 1000;
uu = linspace(0,1,K);
[U1,U2] = meshgrid(uu,uu);
alpha = 10;
copType = 'Frank';
c_density = copulapdf(copType,[U1(:) U2(:)],alpha);
c_density = reshape(c_density,K,K);
uu = copularnd(copType,alpha,M);
vv = empcopularnd(c_density,M);
figure(1);
subplot(1,2,1);
scatter(uu(:,1),uu(:,2)); title('FRANK COPULA Built-in')
subplot(1,2,2);
scatter(vv(:,1),vv(:,2)); title('FRANK COPULA Generated')
pause;

copType = 'Gumbel';
c_density = copulapdf(copType,[U1(:) U2(:)],alpha);
c_density = reshape(c_density,K,K);
uu = copularnd(copType,alpha,M);
vv = empcopularnd(c_density,M);
figure(1);
subplot(1,2,1);
scatter(uu(:,1),uu(:,2)); title('GUMBEL COPULA Built-in')
subplot(1,2,2);
scatter(vv(:,1),vv(:,2)); title('GUMBEL COPULA Generated')
pause;

copType = 'Clayton';
c_density = copulapdf(copType,[U1(:) U2(:)],alpha);
c_density = reshape(c_density,K,K);
uu = copularnd(copType,alpha,M);
vv = empcopularnd(c_density,M);
figure(1);
subplot(1,2,1);
scatter(uu(:,1),uu(:,2)); title('CLAYTON COPULA Built-in')
subplot(1,2,2);
scatter(vv(:,1),vv(:,2)); title('CLAYTON COPULA Generated')
pause;

copType = 'Gaussian';
Rho = [1 0.8; 0.8 1];
c_density = copulapdf(copType,[U1(:) U2(:)],Rho);
c_density = reshape(c_density,K,K);
uu = copularnd(copType,Rho,M);
vv = empcopularnd(c_density,M);
figure(1);
subplot(1,2,1);
scatter(uu(:,1),uu(:,2)); title('Gaussian COPULA Built-in')
subplot(1,2,2);
scatter(vv(:,1),vv(:,2)); title('Gaussian COPULA Generated')
pause;

%% Try 3-D tests
M = 1000;
Rho = [1 .4 .2; .4 1 -.8; .2 -.8 1];
uu = linspace(0,1,K);
[U1,U2,U3] = ndgrid(uu);
c_density = copulapdf('Gaussian',[U1(:) U2(:) U3(:)],Rho);
c_density = reshape(c_density,[K,K,K]);
vv = empcopularnd(c_density,M);

U = copularnd('Gaussian',Rho,M);
h1 = subplot(1,2,1);
plot3(U(:,1),U(:,2),U(:,3),'.')
grid on
view([-55, 15])
xlabel('U1')
ylabel('U2')
zlabel('U3')
title('Built-in')

h2 = subplot(1,2,2);
plot3(vv(:,1),vv(:,2),vv(:,3),'.')
grid on
view([-55, 15])
xlabel('U1')
ylabel('U2')
zlabel('U3')
title('Generated')

hlink = linkprop([h1,h2],{'CameraPosition','CameraUpVector'});
rotate3d on

%% 
% Test full integration as follows:
%  1.) Generate samples from a copula using copularnd
%  2.) Compute the empirical copula density using empcopulapdf
%  3.) Generate samples from the empirical copula density
%  4.) Compare the outputs of Steps #1 and #3
clear;
clc;

% random sample generation parameters
M = 1000;
alpha = 10;
copType = 'Frank';

% copula density estimation parameters
h = 0.05; 
K = 250;

% generate "true" random samples
uu = copularnd(copType,alpha,M);

% generate samples by first computing empirical density
c_density_hat = empcopulapdf(uu, h, K, 'betak');
uu_gen = empcopularnd(c_density_hat, M);

subplot(1,2,1); scatter(uu(:,1),uu(:,2)); grid on; title(sprintf('Built-in %s Samples', copType))
subplot(1,2,2); scatter(uu_gen(:,1),uu_gen(:,2)); grid on; title(sprintf('Generated %s Samples', copType))
pause;

% generate "true" random samples
copType = 'Gumbel';
uu = copularnd(copType,alpha,M);

% generate samples by first computing empirical density
c_density_hat = empcopulapdf(uu, h, K, 'betak');
uu_gen = empcopularnd(c_density_hat, M);

subplot(1,2,1); scatter(uu(:,1),uu(:,2)); grid on; title(sprintf('Built-in %s Samples', copType))
subplot(1,2,2); scatter(uu_gen(:,1),uu_gen(:,2)); grid on; title(sprintf('Generated %s Samples', copType))
pause;

% generate "true" random samples
copType = 'Clayton';
uu = copularnd(copType,alpha,M);

% generate samples by first computing empirical density
c_density_hat = empcopulapdf(uu, h, K, 'betak');
uu_gen = empcopularnd(c_density_hat, M);

subplot(1,2,1); scatter(uu(:,1),uu(:,2)); grid on; title(sprintf('Built-in %s Samples', copType))
subplot(1,2,2); scatter(uu_gen(:,1),uu_gen(:,2)); grid on; title(sprintf('Generated %s Samples', copType))
pause;

copType = 'Gaussian';
Rho = [1 0.8; 0.8 1];
uu = copularnd(copType,Rho,M);

% generate samples by first computing empirical density
c_density_hat = empcopulapdf(uu, h, K, 'betak');
uu_gen = empcopularnd(c_density_hat, M);

subplot(1,2,1); scatter(uu(:,1),uu(:,2)); grid on; title(sprintf('Built-in %s Samples', copType))
subplot(1,2,2); scatter(uu_gen(:,1),uu_gen(:,2)); grid on; title(sprintf('Generated %s Samples', copType))
pause;

%% 
M = 1000;
copType = 'Gaussian';
Rho = [1 .4 .2; .4 1 -.8; .2 -.8 1];
uu = copularnd(copType,Rho,M);
h = .01; K = 50;
c_density_hat = empcopulapdf(uu, h, K, 'betak');
uu_gen = empcopularnd(c_density_hat, M);

h1 = subplot(1,2,1);
plot3(uu(:,1),uu(:,2),uu(:,3),'.')
grid on
view([-55, 15])
xlabel('U1')
ylabel('U2')
zlabel('U3')
title('Built-in Gaussian Samples')

h2 = subplot(1,2,2);
plot3(uu_gen(:,1),uu_gen(:,2),uu_gen(:,3),'.')
grid on
view([-55, 15])
xlabel('U1')
ylabel('U2')
zlabel('U3')
title('Generated Gaussian Samples')

hlink = linkprop([h1,h2],{'CameraPosition','CameraUpVector'});
rotate3d on