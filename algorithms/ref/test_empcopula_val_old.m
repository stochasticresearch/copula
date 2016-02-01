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

% tests the empcopula_val function

clear;
clc;

error_tol = 0.02;

n = 1000;
rho = .7;
Z = mvnrnd([0 0], [1 rho; rho 1], n);
U_real = normcdf(Z);
X_real = [gaminv(U_real(:,1),2,1) tinv(U_real(:,2),5)];

K = 200;
[C,U,c] = empcopula(X_real,K); U1 = U{1}; U2 = U{2};
cc = c{end};
% reconstruct C and c from empcopula_val function, and compare the plots
% and errors between the two to verify accuracy of empcopula_val function
for ii=1:K
    for jj=1:K
        U_1 = U1(ii,jj); U_2 = U2(ii,jj);
        U_query = [U_1 U_2]';
        [C_u, c_u] = empcopula_val(C,cc,U_query);
        if( abs(C(ii,jj)-C_u) > error_tol )  
            fprintf('***********************************************\n');
            fprintf('Mismatch: U_1=%f U_2=%f, C(U_1,U_2)=%f C_u=%f\n', U_1,U_2,C(ii,jj),C_u);
        end
        if( abs(cc(ii,jj)-c_u) > error_tol )
            fprintf('Mismatch: U_1=%f U_2=%f, c(U_1,U_2)=%f c_u=%f\n', ...
                U_1,U_2,cc(ii,jj),c_u);
            fprintf('***********************************************\n');
        end
    end
end
fprintf('Test complete!\n');