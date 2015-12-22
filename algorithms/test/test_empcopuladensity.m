% test script for empcopdens_betak
clear;
clc;

K = 25;

copulaType = 'Clayton';

u = linspace(0,1,K);
[U1,U2] = ndgrid(u);
c3 = copulapdf(copulaType, [U1(:) U2(:)],5);
c3 = reshape(c3, K,K);
% c3 = c3./sum(c3(:));
h3 = subplot(1,3,3);
surf(U1,U2,c3);
xlabel('u1')
ylabel('u2')

X = copularnd(copulaType,5,1000);

h = .05;
[c1] = empcopdens_betak_v2(X(:,1), X(:,2), h, K);
% c1 = c1./sum(c1(:));
h1 = subplot(1,3,1);
surf(U1,U2,c1);
xlabel('u1')
ylabel('u2')

c2 = empcopuladensity(X, h, K, 'betak');
% c2 = c2./sum(c2(:));
h2 = subplot(1,3,2);
surf(U1,U2,c2);
xlabel('u1')
ylabel('u2')

mse = mean((c2(:)-c3(:)).^2);
fprintf('MSE = %f\n', mse);

hlink = linkprop([h1,h2,h3],{'CameraPosition','CameraUpVector'});
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
% c1 = empcopdens_betak_3d_v2(U(:,1), U(:,2), U(:,3), h, K);

c2 = empcopuladensity(U, h, K, 'betak');

u = linspace(0,1,K);
[U1,U2,U3] = ndgrid(u);
c3 = copulapdf('Gaussian', [U1(:) U2(:) U3(:)],Rho);
c3 = reshape(c3, [K,K,K]);

% make sure we have our orientation properly by manually generating 2-D
% copula also
[UU1,UU2] = ndgrid(u);
c3_u1u2 = copulapdf('Gaussian', [UU1(:) UU2(:)], [1 0.4; 0.4 1]); c3_u1u2 = reshape(c3_u1u2,[K,K]);
c3_u2u3 = copulapdf('Gaussian', [UU1(:) UU2(:)], [1 -0.8; -0.8 1]); c3_u2u3 = reshape(c3_u2u3,[K,K]);
c3_u1u3 = copulapdf('Gaussian', [UU1(:) UU2(:)], [1 0.2; 0.2 1]); c3_u1u3 = reshape(c3_u1u3,[K,K]);

h1 = subplot(3,3,1);
surf(UU1,UU2,squeeze(sum(c2,3))); xlabel('u_1'); ylabel('u_2')
h2 = subplot(3,3,2);
surf(UU1,UU2,squeeze(sum(c2,2))); xlabel('u_1'); ylabel('u_3')
title('EMPCOPULADENSITY')
h3 = subplot(3,3,3);
surf(UU1,UU2,squeeze(sum(c2,1))); xlabel('u_2'); ylabel('u_3')

h4 = subplot(3,3,4);
surf(squeeze(sum(c3,3))); xlabel('u_1'); ylabel('u_2')
h5 = subplot(3,3,5);
surf(squeeze(sum(c3,2))); xlabel('u_1'); ylabel('u_3')
title('ACTUAL')
h6 = subplot(3,3,6);
surf(squeeze(sum(c3,1))); xlabel('u_2'); ylabel('u_3')

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
Rho_C1 = [1 .4 .2; .4 1 -.8; .2 -.8 1];
Z = mvnrnd([0 0 0], Rho_C1, M);
U_C1 = normcdf(Z,0,1);

% Generate samples from C2 (B,D) [Clayton Copula]
U_C2_1 = U_C1(:,2); c2_alpha = 2; p = rand(M,1);
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

numPts = 100;
X_vec = [continueRv(X(:,4)) X(:,2)];
F1 = ksdensity(X_vec(:,1),linspace(min(X_vec(:,1)),max(X_vec(:,1)), numPts) ,'support','positive','function','cdf');
F2 = ksdensity(X_vec(:,2),linspace(min(X_vec(:,2)),max(X_vec(:,2)), numPts) ,'support','positive','function','cdf');
FX_vec = [F1' F2'];

c1_ref = empcopuladensity(U_C2, h, K, 'betak');
c1_proper = empcopuladensity(FX_vec, h, K, 'betak');

u = linspace(0,1,K);
[U1,U2] = ndgrid(u);

h1 = subplot(1,2,1); surf(U1,U2,c1_ref); title('Reference')
h2 = subplot(1,2,2); surf(U1,U2,c1_proper); title('Estimated w/ KSDENSITY')
hlink = linkprop([h1,h2],{'CameraPosition','CameraUpVector'});
rotate3d on