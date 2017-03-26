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

%% Test 2-D
clear;
clc;

K = 25;
copulaType = 'Frank';
alpha = 5;
Rho = [1 -0.3; -0.3, 1];

u = linspace(0.01,0.99,K);
[U1,U2] = ndgrid(u);
C = copulacdf(copulaType, [U1(:) U2(:)],alpha);
C = reshape(C, K,K);

U = copularnd(copulaType,alpha,1000);
C_est = empcopulacdf(U, K, 'deheuvels');

h1 = subplot(1,2,1);
surf(U1,U2,C);
xlabel('u1')
ylabel('u2')
title('C')

h2 = subplot(1,2,2);
surf(U1,U2,C_est);
xlabel('u1')
ylabel('u2')
title('C Est')

mse = mean((C_est(:)-C(:)).^2);
fprintf('MSE = %f\n', mse);

hlink = linkprop([h1,h2],{'CameraPosition','CameraUpVector'});
rotate3d on

%% Test 3-D
clear;
clc;

K = 25;
copulaType = 'Gaussian';
Rho = [1 -0.3 0.6; -0.3 1 0.2; 0.6 0.2 1];

u = linspace(0.01,0.99,K);
[U1,U2,U3] = ndgrid(u);
C = copulacdf(copulaType, [U1(:) U2(:) U3(:)],Rho);
C = reshape(C, [K,K,K]);

U = copularnd(copulaType,Rho,1000);
C_est = empcopulacdf(U, K, 'deheuvels');

% make sure we have our orientation properly by manually generating 2-D
% copula also
[UU1,UU2] = ndgrid(u);
C3_u1u2 = copulacdf('Gaussian', [UU1(:) UU2(:)], [1 -0.3; -0.3 1]); C3_u1u2 = reshape(C3_u1u2,[K,K]);
C3_u2u3 = copulacdf('Gaussian', [UU1(:) UU2(:)], [1 0.2; 0.2 1]); C3_u2u3 = reshape(C3_u2u3,[K,K]);
C3_u1u3 = copulacdf('Gaussian', [UU1(:) UU2(:)], [1 0.6; 0.6 1]); C3_u1u3 = reshape(C3_u1u3,[K,K]);

h1 = subplot(3,3,1);
surf(UU1,UU2,squeeze(sum(C_est,3))); xlabel('u_1'); ylabel('u_2')
h2 = subplot(3,3,2);
surf(UU1,UU2,squeeze(sum(C_est,2))); xlabel('u_1'); ylabel('u_3')
title('empcopulacdf')
h3 = subplot(3,3,3);
surf(UU1,UU2,squeeze(sum(C_est,1))); xlabel('u_2'); ylabel('u_3')

h4 = subplot(3,3,4);
surf(UU1,UU2,squeeze(sum(C,3))); xlabel('u_1'); ylabel('u_2')
h5 = subplot(3,3,5);
surf(UU1,UU2,squeeze(sum(C,2))); xlabel('u_1'); ylabel('u_3')
title('ACTUAL MARGINALIZED')
h6 = subplot(3,3,6);
surf(UU1,UU2,squeeze(sum(C,1))); xlabel('u_2'); ylabel('u_3')

h7 = subplot(3,3,7);
surf(UU1,UU2,C3_u1u2); xlabel('u_1'); ylabel('u_2')
h8 = subplot(3,3,8);
surf(UU1,UU2,C3_u2u3); xlabel('u_1'); ylabel('u_3')
title('ACTUAL MARGINAL')
h9 = subplot(3,3,9);
surf(UU1,UU2,C3_u1u3); xlabel('u_2'); ylabel('u_3')

rotate3d on

%% Generate Density from 2-D Frank Copula by first taking CDF and differentiating,
% and using empcopulapdf, to see hte difference

clear;
clc;

K = 25;
copulaType = 'Frank';
alpha = 5;
Rho = [1 -0.3; -0.3, 1];

u = linspace(0.01,0.99,K);
[U1,U2] = ndgrid(u);
c = copulapdf(copulaType, [U1(:) U2(:)],alpha);
c = reshape(c, K,K);

U = copularnd(copulaType,alpha,1000);
C_est = empcopulacdf(U, K, 'deheuvels');
[cc] = gradient(C_est);
[~,c_est_grad] = gradient(cc);

% estimate w/ empcopulapdf
h = 0.05;
K = 25;
c_direct = empcopulapdf(U, h, K, 'betak');

% normalize all of the c's for plotting effect
c = c./max(c(:));
c_est_grad = c_est_grad./max(c_est_grad(:));
c_direct = c_direct./max(c_direct(:));

h1 = subplot(1,3,1); surf(U1,U2,c); xlabel('u'); ylabel('v'); title('$$c_{actual} (\mathbf{u})$$', 'interpreter', 'latex'); grid on;
h2 = subplot(1,3,2); surf(U1,U2,c_est_grad); xlabel('u'); ylabel('v'); 
title('$$ \hat{c}(\mathbf{u}) = \frac{\partial \hat{C}(\mathbf{u})}{d \mathbf{u}}$$','interpreter','latex'); grid on;
h3 = subplot(1,3,3); surf(U1,U2,c_direct); xlabel('u'); ylabel('v'); 
title('$$ \hat{c}_\beta(\mathbf{u}) $$','interpreter','latex'); grid on;

%% A small test of functional dependnece

clear;
clc;
M = 5000;
K = 100; h = 0.1;
minVal = 0.0; maxVal = 1.0;
uu = linspace(minVal,maxVal,K); vv = linspace(minVal,maxVal,K);

x = rand(M,1); x = sort(x);
y = (x-0.5).^2;
u = pobs(x); v = pobs(y);
U = [u v];
C = empcopulacdf(U,K,'deheuvels');
c = empcopulapdf(U,h,K,'betak-matlab');
CC = C'; cc = c';       % CC and cc are v/u rather than u/v oriented.

% now construct the copula separately for the discordant section, and the
% concordant section
I = find(x<=0.5);
x1 = x(I); y1 = y(I);
I = find(x>0.5);
x2 = x(I); y2 = y(I);
u1 = pobs(x1); v1 = pobs(y1); U1 = [u1 v1];
u2 = pobs(x2); v2 = pobs(y2); U2 = [u2 v2];
C1 = empcopulacdf(U1,K,'deheuvels'); c1 = empcopulapdf(U1,h,K,'betak-matlab');
CC1 = C1'; cc1 = c1';
C2 = empcopulacdf(U2,K,'deheuvels'); c2 = empcopulapdf(U2,h,K,'betak-matlab');
CC2 = C2'; cc2 = c2';

% integrate from u=[0,uu], v=[0,1]
tauVec = zeros(1,K);
tauVec1 = zeros(1,K); tauVec2 = zeros(1,K);
umaxVec = uu;
for ii=1:K
    matSubset = CC(:,1:ii).*cc(:,1:ii);
    matSubset1 = CC1(:,1:ii).*cc1(:,1:ii);
    matSubset2 = CC2(:,1:ii).*cc2(:,1:ii);
    
    uIntegralRange = uu(1:ii);
    
    tauVec(ii) = 4*trapz(vv, trapz(uIntegralRange, matSubset, 2))-1;
    tauVec1(ii) = 4*trapz(vv, trapz(uIntegralRange, matSubset1, 2))-1;
    tauVec2(ii) = 4*trapz(vv, trapz(uIntegralRange, matSubset2, 2))-1;
end

subplot(4,2,1); surf(uu,vv,CC); xlabel('v'); ylabel('u'); zlabel('CC');
subplot(4,2,2); surf(uu,vv,cc); xlabel('v'); ylabel('u'); zlabel('cc');
subplot(4,2,3); surf(uu,vv,CC.*cc); xlabel('v'); ylabel('u'); zlabel('CC*cc');
subplot(4,2,4); plot(umaxVec, tauVec); grid on; xlabel('umax'); ylabel('Q(\tau) | y=x^2');
subplot(4,2,5); surf(uu,vv,CC1.*cc1); xlabel('v'); ylabel('u'); zlabel('CC1*cc1');
subplot(4,2,6); plot(umaxVec, tauVec1); grid on; xlabel('umax'); ylabel('Q(\tau) | y=x^2 -M');
subplot(4,2,7); surf(uu,vv,CC2.*cc2); xlabel('v'); ylabel('u'); zlabel('CC2*cc2');
subplot(4,2,8); plot(umaxVec, tauVec2); grid on; xlabel('umax'); ylabel('Q(\tau) | y=x^2 +M');
