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

% test script for empcopdens_betak
clear;
clc;

K = 25;

copulaType = 'Clayton';

u = linspace(0,1,K);
[U1,U2] = ndgrid(u);
c2 = copulapdf(copulaType, [U1(:) U2(:)],5);
c2 = reshape(c2, K,K);
h3 = subplot(1,2,1);
surf(U1,U2,c2);
xlabel('u1')
ylabel('u2')
title('Matlab Implementation')

X = copularnd(copulaType,5,1000);
h = 0.05;

c1 = empcopulapdf(X, h, K, 'betak');
h2 = subplot(1,2,2);
surf(U1,U2,c1);
xlabel('u1')
ylabel('u2')
title('Estimated w/ Beta-Kernel')

mse = mean((c1(:)-c2(:)).^2);
fprintf('MSE = %f\n', mse);

hlink = linkprop([h1,h2],{'CameraPosition','CameraUpVector'});
rotate3d on
%%
clear;
clc;

% Test the 3-D
h = 0.05;
K = 25;
M = 1000;
D = 3;
Rho = [1 .4 .2; .4 1 -.8; .2 -.8 1];
Z = mvnrnd([0 0 0], Rho, M);
U = normcdf(Z,0,1);

c1 = empcopulapdf(U, h, K, 'betak');

u = linspace(0,1,K);
[U1,U2,U3] = ndgrid(u);
c2 = copulapdf('Gaussian', [U1(:) U2(:) U3(:)],Rho);
c2 = reshape(c2, [K,K,K]);

% make sure we have our orientation properly by manually generating 2-D
% copula also
[UU1,UU2] = ndgrid(u);
c3_u1u2 = copulapdf('Gaussian', [UU1(:) UU2(:)], [1 0.4; 0.4 1]); c3_u1u2 = reshape(c3_u1u2,[K,K]);
c3_u2u3 = copulapdf('Gaussian', [UU1(:) UU2(:)], [1 -0.8; -0.8 1]); c3_u2u3 = reshape(c3_u2u3,[K,K]);
c3_u1u3 = copulapdf('Gaussian', [UU1(:) UU2(:)], [1 0.2; 0.2 1]); c3_u1u3 = reshape(c3_u1u3,[K,K]);

h1 = subplot(3,3,1);
surf(UU1,UU2,squeeze(sum(c1,3))); xlabel('u_1'); ylabel('u_2')
h2 = subplot(3,3,2);
surf(UU1,UU2,squeeze(sum(c1,2))); xlabel('u_1'); ylabel('u_3')
title('empcopulapdf')
h3 = subplot(3,3,3);
surf(UU1,UU2,squeeze(sum(c1,1))); xlabel('u_2'); ylabel('u_3')

h4 = subplot(3,3,4);
surf(squeeze(sum(c2,3))); xlabel('u_1'); ylabel('u_2')
h5 = subplot(3,3,5);
surf(squeeze(sum(c2,2))); xlabel('u_1'); ylabel('u_3')
title('ACTUAL')
h6 = subplot(3,3,6);
surf(squeeze(sum(c2,1))); xlabel('u_2'); ylabel('u_3')

h7 = subplot(3,3,7);
surf(UU1,UU2,c3_u1u2); xlabel('u_1'); ylabel('u_2')
h8 = subplot(3,3,8);
surf(UU1,UU2,c3_u1u3); xlabel('u_1'); ylabel('u_3')
title('ACTUAL MARGINAL')
h9 = subplot(3,3,9);
surf(UU1,UU2,c3_u2u3); xlabel('u_2'); ylabel('u_3')


% mse = mean((c2(:)-c1(:)).^2);
% fprintf('MSE 3D = %f\n', mse);

% hlink = linkprop([h1,h2,h3,h7,h8,h9],{'CameraPosition','CameraUpVector'});
rotate3d on

%%
clear;
clc;

K = 25;
h = 0.05;

M = 1000;
D = 5;

% Generate samples from C1 (A,B,C) [Gaussian Copula]
Rho = [1 .4 .2; .4 1 -.8; .2 -.8 1];
Z = mvnrnd([0 0 0], Rho, M);
U_C1 = normcdf(Z,0,1);

% Generate samples from C2 (B,D) [Clayton Copula]
U_C2_1 = U_C1(:,2); c2_alpha = 1; p = rand(M,1);
U_C2_2 = U_C2_1.*(p.^(-c2_alpha./(1+c2_alpha)) - 1 + U_C2_1.^c2_alpha).^(-1./c2_alpha);
U_C2 = [U_C2_1 U_C2_2];

% Generate samples from C3 (C,E) [Clayton Copula]
U_C3_1 = U_C1(:,3); c3_alpha = 4; p = rand(M,1);
U_C3_2 = U_C3_1.*(p.^(-c3_alpha./(1+c3_alpha)) - 1 + U_C3_1.^c3_alpha).^(-1./c3_alpha);
U_C3 = [U_C3_1 U_C3_2];

U = [U_C1 U_C2(:,2) U_C3(:,2)];

X = [gaminv(U(:,1),2,1) ...
       betainv(U(:,2),2,2) ...
       unidinv(U(:,3),5) ...
       unidinv(U(:,4),3) ...
       norminv(U(:,5),0,1)];

X_vec = [continueRv(X(:,4)) X(:,2)];

% create pseudo-observations
U_in = pseudoobs(X_vec);

c1_ref = empcopulapdf(U_C2, h, K, 'betak');
c1_proper = empcopulapdf(U_in, h, K, 'betak');

u = linspace(0,1,K);
[U1,U2] = ndgrid(u);

h1 = subplot(1,2,1); surf(U1,U2,c1_ref); title('Reference')
h2 = subplot(1,2,2); surf(U1,U2,c1_proper); title('Estimated w/ KSDENSITY')
hlink = linkprop([h1,h2],{'CameraPosition','CameraUpVector'});
rotate3d on

%% Test empcopulapdf heuristically w/ discrete marginals & Gaussian Copula
clear;
clc;

M = 1000;
D = 2;

% Generate samples from C1 (A,B) [Gaussian Copula]
Rho = [1 -0.8; -0.8 1];
U = copularnd('Gaussian', Rho, M);

X = [unidinv(U(:,1),3) ...
     unidinv(U(:,2),4)];

% transform the data by dithering according to Michel, Denuit, Neslehova
X_xform = continueRv(X);

% generate pseudoobservations from X_xform
U_in = pseudoobs(X_xform);

% estimate the copula density
h = 0.05;
K = 25;
c = empcopulapdf(U_in, h, K, 'betak'); 

u = linspace(0,1,K);
[U1,U2] = ndgrid(u);
c_expect = copulapdf('Gaussian', [U1(:) U2(:)],Rho);
c_expect = reshape(c_expect,K,K);


h1 = subplot(1,2,1); surf(U1,U2,c_expect); grid on; title('Reference'); xlabel('U_1'); ylabel('U_2');
h2 = subplot(1,2,2); surf(U1,U2,c); grid on; title('Estimated'); xlabel('U_1'); ylabel('U_2');
linkprop([h1,h2],{'CameraPosition','CameraUpVector'}); rotate3d on;

%% Test empcopulapdf heuristically w/ discrete marginals & Gumbel Copula
clear;
clc;

M = 1000;
D = 2;

alpha = 9.5; copType = 'Clayton';
U = copularnd(copType, alpha, M);

X = [unidinv(U(:,1),3) ...
     unidinv(U(:,2),4)];

% transform the data by dithering according to Michel, Denuit, Neslehova
X_xform = continueRv(X);

% generate pseudoobservations from X_xform
U_in = pseudoobs(X_xform);

% estimate the copula density
h = 0.02;
K = 50;
c = empcopulapdf(U_in, h, K, 'betak');

[U1,U2] = ndgrid(linspace(0,1,K));
c_expect = copulapdf(copType, [U1(:) U2(:)],alpha);
c_expect = reshape(c_expect,K,K);

h1 = subplot(1,2,1); surf(U1,U2,c_expect); grid on; title('Reference'); xlabel('U_1'); ylabel('U_2');
h2 = subplot(1,2,2); surf(U1,U2,c); grid on; title('Estimated'); xlabel('U_1'); ylabel('U_2');
linkprop([h1,h2],{'CameraPosition','CameraUpVector'}); rotate3d on;

%% Test empcopulapdf heuristically w/ discrete marginals & Gumbel Copula
clear;
clc;

M = 1000;
D = 2;

Rho = [1 0.4; 0.4 1]; copType = 'Gaussian';
U = copularnd(copType, Rho, M);

X1 = unidinv(U(:,1),5);
X2 = unidinv(U(:,2),5);
X = [X1 X2];
 
[nn, ctrs] = hist3(X);
probMatrix = zeros(5,5);
probMatrix(1,1) = nn(1,1)/M; probMatrix(1,2) = nn(1,3)/M; probMatrix(1,3) = nn(1,6)/M; probMatrix(1,4) = nn(1,8)/M; probMatrix(1,5) = nn(1,10)/M;
probMatrix(2,1) = nn(3,1)/M; probMatrix(2,2) = nn(3,3)/M; probMatrix(2,3) = nn(3,6)/M; probMatrix(2,4) = nn(3,8)/M; probMatrix(2,5) = nn(3,10)/M;
probMatrix(3,1) = nn(6,1)/M; probMatrix(3,2) = nn(6,3)/M; probMatrix(3,3) = nn(6,6)/M; probMatrix(3,4) = nn(6,8)/M; probMatrix(3,5) = nn(6,10)/M;
probMatrix(4,1) = nn(8,1)/M; probMatrix(4,2) = nn(8,3)/M; probMatrix(4,3) = nn(8,6)/M; probMatrix(4,4) = nn(8,8)/M; probMatrix(4,5) = nn(8,10)/M;
probMatrix(5,1) = nn(10,1)/M; probMatrix(5,2) = nn(10,3)/M; probMatrix(5,3) = nn(10,6)/M; probMatrix(5,4) = nn(10,8)/M; probMatrix(5,5) = nn(10,10)/M;

X_xform = continueRv(X);
U_in = pseudoobs(X_xform);

% estimate the copula density
h = 0.02;
K = 25;
c_est = empcopulapdf(U_in, h, K, 'betak');

[U1,U2] = ndgrid(linspace(0,1,K));
c_actual = copulapdf(copType, [U1(:) U2(:)], Rho);
c_actual = reshape(c_actual,K,K);

h1 = subplot(1,3,1); hist3(X);                   title('Joint Distribution')
h2 = subplot(1,3,2); surf(U1,U2,c_est);          title('Estimated Copula')
h3 = subplot(1,3,3); surf(U1,U2,c_actual);       title('Actual Copula')
linkprop([h1,h2,h3],{'CameraPosition','CameraUpVector'}); rotate3d on;

% now, since we have the full PDF, we can compare how accurate hte
% estimated copula is to the true, with the definition: 
% c(u1,u2) = f(x1,x2)/[ f(x1) * f(x2) ]

[FX_in1, domain1] = ecdf(X1); FX_in1 = FX_in1(2:end); domain1 = domain1(2:end);
empInfoObj1 = rvEmpiricalInfo(domain1, [], FX_in1);

[FX_in2, domain2] = ecdf(X2); FX_in2 = FX_in2(2:end); domain2 = domain2(2:end);
empInfoObj2 = rvEmpiricalInfo(domain2, [], FX_in2);

for x1=1:5
    for x2=1:5
        f_val = probMatrix(x1,x2)/(0.33*0.33);  % f(x1,x2)/[f(x1)*f(x2)]
        u1 = empInfoObj1.queryDistribution(x1);
        u2 = empInfoObj2.queryDistribution(x2);
        if(abs(u1-1)<=.01)
            u1 = 0.99;
        elseif(abs(u1-.01)<=0.01)
            u1 = 0.01;
        end
        if(abs(u2-1)<=.01)
            u2 = 0.99;
        elseif(abs(u2-.01)<=0.01)
            u2 = 0.01;
        end
        u = [u1 u2];
        
        c_est_val = empcopula_val(c_est, u); 
        c_actual_val = copulapdf('Gaussian', u, Rho);
        
        f_val_mat(x1,x2) = f_val;
        c_est_mat(x1,x2) = c_est_val;
        c_actual_mat(x1,x2) = c_actual_val;
    end
end

f_val_mat
c_est_mat
c_actual_mat

est_err_mat = (f_val_mat-c_est_mat).^2
actual_err_mat = (f_val_mat-c_actual_mat).^2

%% Test empcopulapdf heuristically w/ discrete marginals & Gumbel Copula
clear;
clc;

M = 1000;
D = 2;

Rho = [1 0.4; 0.4 1]; copType = 'Gaussian';
U = copularnd(copType, Rho, M);

% Create a multivariate gaussian :) easy for testing
X = [norminv(U(:,1),0,1) ...
     norminv(U(:,2),0,1)];
U_in = pseudoobs(X);

% estimate the copula density
h = 0.01;
K = 50;
c_est = empcopulapdf(U_in, h, K, 'betak');

ii = 1;
for x1=-3:.1:3
    jj = 1;
    for x2=-3:.1:3
        f_val = mvnpdf([x1,x2], [0 0], Rho)/(normpdf(x1)*normpdf(x2)); u=[normcdf(x1) normcdf(x2)];
        c_est_val = empcopula_val(c_est, u); 
        c_actual_val = copulapdf('Gaussian', u, Rho);
        
        f_val_mat(ii,jj) = f_val;
        c_est_mat(ii,jj) = c_est_val;
        c_actual_mat(ii,jj) = c_actual_val;
        jj = jj + 1;
    end
    ii = ii + 1;
end

[U1,U2] = ndgrid(linspace(0,1,length(-3:.1:3)));

h1 = subplot(2,2,1); surf(U1,U2,f_val_mat); title('f calc'); xlabel('U_1'); ylabel('U_2')
h2 = subplot(2,2,2); surf(U1,U2,c_actual_mat); title('c actual'); xlabel('U_1'); ylabel('U_2')
h3 = subplot(2,2,3); surf(U1,U2,c_est_mat); title('c est'); xlabel('U_1'); ylabel('U_2')
h4 = subplot(2,2,4); surf(U1,U2,(c_est_mat-f_val_mat).^2); title('(c est - f calc)^2'); xlabel('U_1'); ylabel('U_2')
linkprop([h1,h2,h3,h4],{'CameraPosition','CameraUpVector'}); rotate3d on;
