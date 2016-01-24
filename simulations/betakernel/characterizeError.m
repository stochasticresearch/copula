% Script which generates multivariate copulas of different types with discrete 
% marginal data, and then reconstructs the copula using beta kernels. 
% It then characterizes the error between the constructed copula and 
%  1.) the actual copula
%  2.) c(F(x1) ... f(Xn)) and f(x1 ... xn) / [ f(x1) * ... * f(xn) ]

clear;
clc;

K = 25;

copulaType = 'Clayton';

u = linspace(0,1,K);
[U1,U2] = ndgrid(u);
c2 = copulapdf(copulaType, [U1(:) U2(:)],5);
c2 = reshape(c2, K,K);
h3 = subplot(1,3,3);
surf(U1,U2,c2);
xlabel('u1')
ylabel('u2')

alpha = 5;
M = 1000;

X = copularnd(copulaType,alpha,M);

h = .05;
[c1] = empcopdens_betak_v2(X(:,1), X(:,2), h, K);
h1 = subplot(1,3,1);
surf(U1,U2,c1);
xlabel('u1')
ylabel('u2')

c1 = empcopuladensity(X, h, K, 'betak');
h2 = subplot(1,3,2);
surf(U1,U2,c1);
xlabel('u1')
ylabel('u2')

mse = mean((c1(:)-c2(:)).^2);
fprintf('MSE = %f\n', mse);

hlink = linkprop([h1,h2,h3],{'CameraPosition','CameraUpVector'});
rotate3d on