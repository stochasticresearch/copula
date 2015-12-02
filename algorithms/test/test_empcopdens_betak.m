% test script for empcopdens_betak
clear;
clc;

K = 50;

u = linspace(0,1,K);
[U1,U2] = meshgrid(u);
c3 = copulapdf('Gumbel', [U1(:) U2(:)],5);
c3 = reshape(c3, K,K);
h3 = subplot(1,3,3);
surf(U1,U2,c3);
xlabel('u1')
ylabel('u2')

X = copularnd('Gumbel',5,1000);

h = .1;
[c1] = empcopdens_betak_v2(X(:,1), X(:,2), h, K);
h1 = subplot(1,3,1);
surf(U1,U2,c1);
xlabel('u1')
ylabel('u2')

c2 = empcopuladensity(X, h, K, 'betak');
h2 = subplot(1,3,2);
surf(U1,U2,c2);
xlabel('u1')
ylabel('u2')

mse = mean((c2(:)-c1(:)).^2);
fprintf('MSE = %f\n', mse);

% hlink = linkprop([h1,h2,h3],{'CameraPosition','CameraUpVector'});
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
