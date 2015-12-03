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
