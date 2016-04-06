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

% A script to help us understand the value of the edges for different
% copulas

%% 2-D tests

startPnt = 0.90;
interval = .001;
u1 = startPnt:interval:1.0;
u2 = u1;

[U1,U2] = ndgrid(u1,u2);
U = [U1(:) U2(:)];

% compute gaussian copula values
Rho = [1 0.4; 0.4 1];
gaussian_pdf = copulapdf('Gaussian', U, Rho);
gaussian_pdf = reshape(gaussian_pdf, length(u1), length(u1));
gaussian_pdf_fixed = copulapdf('Gaussian', fixU(U), Rho);
gaussian_pdf_fixed = reshape(gaussian_pdf_fixed, length(u1), length(u1));

alpha = 4;
% compute frank copula values
frank_pdf = copulapdf('Frank', U, alpha);
frank_pdf = reshape(frank_pdf, length(u1), length(u1));
frank_pdf_fixed = copulapdf('Frank', fixU(U), alpha);
frank_pdf_fixed = reshape(frank_pdf_fixed, length(u1), length(u1));

% compute gumbel copula values
gumbel_pdf = copulapdf('Gumbel', U, alpha);
gumbel_pdf = reshape(gumbel_pdf, length(u1), length(u1));
gumbel_pdf_fixed = copulapdf('Gumbel', fixU(U), alpha);
gumbel_pdf_fixed = reshape(gumbel_pdf_fixed, length(u1), length(u1));

% compute clayton copula values
clayton_pdf = copulapdf('Clayton', U, alpha);
clayton_pdf = reshape(clayton_pdf, length(u1), length(u1));
clayton_pdf_fixed = copulapdf('Clayton', fixU(U), alpha);
clayton_pdf_fixed = reshape(clayton_pdf_fixed, length(u1), length(u1));

figure(1);
h1=subplot(2,2,1); surf(U1,U2,gaussian_pdf); title('Gaussian Copula PDF'); grid on; xlabel('u_1'); ylabel('u_2');
h2=subplot(2,2,2); surf(U1,U2,frank_pdf); title('Frank Copula PDF'); grid on; xlabel('u_1'); ylabel('u_2');
h3=subplot(2,2,3); surf(U1,U2,gumbel_pdf); title('Gumbel Copula PDF'); grid on; xlabel('u_1'); ylabel('u_2');
h4=subplot(2,2,4); surf(U1,U2,clayton_pdf); title('Clayton Copula PDF'); grid on; xlabel('u_1'); ylabel('u_2');
linkprop([h1,h2,h3,h4],{'CameraPosition','CameraUpVector'}); rotate3d on;

figure(2);
h1=subplot(2,2,1); surf(U1,U2,gaussian_pdf_fixed); title('Gaussian Copula PDF - Fixed'); grid on; xlabel('u_1'); ylabel('u_2');
h2=subplot(2,2,2); surf(U1,U2,frank_pdf_fixed); title('Frank Copula PDF - Fixed'); grid on; xlabel('u_1'); ylabel('u_2');
h3=subplot(2,2,3); surf(U1,U2,gumbel_pdf_fixed); title('Gumbel Copula PDF - Fixed'); grid on; xlabel('u_1'); ylabel('u_2');
h4=subplot(2,2,4); surf(U1,U2,clayton_pdf_fixed); title('Clayton Copula PDF - Fixed'); grid on; xlabel('u_1'); ylabel('u_2');
linkprop([h1,h2,h3,h4],{'CameraPosition','CameraUpVector'}); rotate3d on;

find(clayton_pdf_fixed==0)
find(gaussian_pdf_fixed==0)
find(frank_pdf_fixed==0)
find(gumbel_pdf_fixed==0)