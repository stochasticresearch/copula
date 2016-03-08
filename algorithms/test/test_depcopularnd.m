% Tests the depcopularnd script

%% Test Gaussian copula CI test, with BNT's Fisher test implementation
clear;
clc;

Rho1 = [1 0.5; 0.5 1];
Rho2 = [1 -0.3; -0.3 1];
M = 1000;

U_init = copularnd('Gaussian', Rho1, M);
U_dep = depcopularnd(U_init(:,1), 2, 'Gaussian', Rho2);

Z1 = U_init(:,1); Z2 = U_dep(:,1);
X = U_init(:,2);
Y = U_dep(:,2);

partialcorr([X Y Z1])
partialcorr([X Y Z2])

Rho3 = [1 0.5 0.3; 0.5 1 -0.2; 0.3 -0.2 1];
U = mvnrnd([0 0 0], Rho3, M);
partialcorr(U)



%% Test the Clayton Copula
alpha = 5;
M = 1000;
N = 2;
[U_init,X_i] = claytoncopularnd(M, N, alpha);
U_dep = depcopularnd(X_i(:,1), N, 'Clayton', alpha+3);

subplot(1,3,1); scatter(U_init(:,1),U_init(:,2)); grid on; title('Clayton Copula Samples');
subplot(1,3,2); scatter(U_dep(:,1), U_dep(:,2)); grid on; title('Clayton Copula Samples [DEP]');
subplot(1,3,3); scatter(U_init(:,1), U_dep(:,1)); grid on; title('DEP comparison');
