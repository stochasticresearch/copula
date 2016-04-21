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

%% cmpare different ways to generate empirical copula cdf (direct vs via pdf)
clear;
clc;

K = 100;

copulaType = 'Clayton';

u = linspace(0,1,K);
[U1,U2] = ndgrid(u);
c = copulapdf(copulaType, [U1(:) U2(:)],5);
C = copulacdf(copulaType, [U1(:) U2(:)],5);
c = reshape(c, K,K); C = reshape(C, K, K);

U = copularnd(copulaType,5,1000); h = 0.01;
c_est_betak = empcopulapdf(U, h, K, 'betak');
C_est_dehuvels = empcopulacdf(U, K, 'deheuvels');
C_est_betakIntegration = c_est_betak;
c_est_dehuvels = C_est_dehuvels;
for dim=1:2
    C_est_betakIntegration = cumtrapz(u,C_est_betakIntegration,dim);
    c_est_dehuvels = diff(c_est_dehuvels, 1, dim);
end
c_est_dehuvels_pad = zeros(size(c_est_dehuvels)+1);
c_est_dehuvels_pad(2:end,2:end) = c_est_dehuvels;
c_est_dehuvels = c_est_dehuvels_pad;

subplot(2,6,[1 2]);
contour(U1,U2,c);
xlabel('U_1'); ylabel('U_2'); title('Matlab PDF'); grid on;

subplot(2,6,[3 4]);
contour(U1,U2,c_est_dehuvels);
xlabel('U_1'); ylabel('U_2'); title('Dehuvels e-PDF'); grid on;

subplot(2,6, [5 6]);
contour(U1,U2,c_est_betak);
xlabel('U_1'); ylabel('U_2'); title('BetaK e-PDF'); grid on;


subplot(2,6,[7 8]);
surf(U1,U2,C); xlabel('U_1'); ylabel('U_2'); title('Matlab CDF'); grid on;

% compute the CDF via Dehuvels method
subplot(2,6,[9 10]);
surf(U1,U2,C_est_dehuvels); xlabel('U_1'); ylabel('U_2'); title('Dehuvels e-CDF'); grid on;

% compute the CDF by integrating the estimated PDF
subplot(2,6,[11 12]);
% C_est_betakIntegration = cumtrapz(u,c_est_betak,1);       % how to integrate selective dimensions
surf(U1,U2,C_est_betakIntegration); xlabel('U_1'); ylabel('U_2'); title('BetaK e-CDF'); grid on;

