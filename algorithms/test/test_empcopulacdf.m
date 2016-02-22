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