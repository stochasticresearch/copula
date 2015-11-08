% Test's the emp_copularnd function w/ various test patterns

clear;
clc;
close all;

% Generate a Joint Distribution of Gamma and T Random Variable's that are
% correlated with rho=0.7
n = 1000;
rho = .7;
Z = mvnrnd([0 0], [1 rho; rho 1], n);
U_real = normcdf(Z);
X_real = [gaminv(U_real(:,1),2,1) tinv(U_real(:,2),5)];

% calculate the empirical copula & generate samples from this copula
K = 200;
D = size(X_real,2);
[ U_gen, ~, U_emp ] = emp_copularnd( X_real, n, K );
X_gen = empdistrnd( U_gen, X_real );

% Setup for Plotting
U_real = U_real';
U_gen = U_gen';
U_emp = U_emp';

d = 1; ds = 2;

% plot original dependency structure
subplot(1,3,1)
scatter(U_real(d,1:n),U_real(ds,1:n),2.5,'*')
axis([-.1 1.1 -.1 1.1])
xlabel(sprintf('U_%d',d))
ylabel(sprintf('U_%d',ds))
title('Real Dependency Structure')
grid on

% plot empirical dependency structure
subplot(1,3,2)
scatter(U_emp(d,1:n),U_emp(ds,1:n),2.5,'*')
axis([-.1 1.1 -.1 1.1])
xlabel(sprintf('U_%d',d))
ylabel(sprintf('U_%d',ds))
title('Empirical Detected Dependency Structure')
grid on

% plot generated dependency structure
subplot(1,3,3)
scatter(U_gen(d,1:n),U_gen(ds,1:n),2.5,'*')
axis([-.1 1.1 -.1 1.1])
xlabel(sprintf('U_%d',d))
ylabel(sprintf('U_%d',ds))
title('Generated Dependency Structure')
grid on

% plot original samples of distribution
figure;
subplot(1,2,1)
scatter(X_real(:,d),X_real(:,ds),2.5,'*')
xlabel(sprintf('X_%d',d))
ylabel(sprintf('X_%d',ds))
title('Real Samples')
grid on

% plot generated samples of the distribution
subplot(1,2,2)
scatter(X_gen(:,d),X_gen(:,ds),2.5,'*')
xlabel(sprintf('X_%d',d))
ylabel(sprintf('X_%d',ds))
title('Generated Samples')
grid on