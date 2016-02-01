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

% Main script which compares the performance of the copula vs the
% conditional gaussian model for hybrid distributions

%% 2-D simulation
clear;
clc;

n = 1000;
[X_real, U_real] = gen_hybriddata_2D( n );
x1_domain = unique(X_real(:,1));
h = histogram(X_real(:,1),length(x1_domain));
x1_values = h.Values/n;
[f2,x2] = ksdensity(X_real(:,2));

% fit X_real to a conditional gaussian model
[x1_domain, x1_multinomial_est, x2_mle_params] = fit_clg_2D( X_real );
% generate samples from the derived generative model
X_clg_gen = gen_samples_clg_2D(x1_domain, x1_multinomial_est, x2_mle_params, n);

% fit X_real w/ empirical copula and empirical marginals
% first continue the random variable
X1_continued = X_real(:,1) + (rand(n,1)-1);
X_transformed = [X1_continued X_real(:,2)];
K = 100;
D = size(X_real,2);

% [ U_gen, ~, U_emp ] = emp_copularnd( X_transformed, n, K );
[C,U,c] = empcopula(X_transformed,K);    
U_gen = empcopularnd(c, n);

X_gen = empdistrnd(U_gen, X_real);

%%%%%%%%%%%%%% VISUALIZATION %%%%%%%%%%%%%%

% visualize the simulation data
subplot(2,2,1)
stem(x1_domain, x1_values)
xlabel('X_1')
ylabel('"Real" f(x_1)')
grid on

subplot(2,2,2)
plot(x2,f2/sum(f2))
xlabel('X_2')
ylabel('"Real" f(x_2)')
grid on

subplot(2,2,3)
scatter(U_real(:,1),U_real(:,2))
xlabel('U_1')
ylabel('U_2')
grid on
title('"Real" Dependency Structure')

subplot(2,2,4)
scatter(X_real(:,1),X_real(:,2))
xlabel('X_1')
ylabel('Y_1')
grid on
title('"Real" Joint Distribution')

figure;
subplot(1,3,1);
scatter(X_real(:,1),X_real(:,2));
grid on
xlabel('X_1')
ylabel('X_2')
title('Empirical Data')

subplot(1,3,2);
scatter(X_gen(:,1),X_gen(:,2));
grid on
xlabel('X_1')
ylabel('X_2')
title('Copula Generative Model')

subplot(1,3,3);
scatter(X_clg_gen(:,1),X_clg_gen(:,2));
grid on
xlabel('X_1')
ylabel('X_2')
title('CLG Generative Model')

fprintf('Copulas WIN\n')

%% 2-D all continuous simulation
clear;
clc;

n = 1000;
rho = .7;
Z = mvnrnd([0 0], [1 rho; rho 1], n);
U_real = normcdf(Z);
X_real = [gaminv(U_real(:,1),2,1) tinv(U_real(:,2),5)];

K = 100;
n_gen = 1000;
D = size(X_real,2);

% [ U_gen, ~, U_emp ] = emp_copularnd( X_real, n_gen, K );
[C,U,c] = empcopula(X_real,K);    
U_gen = empcopularnd(c, n_gen);

X_gen = empdistrnd(U_gen, X_real);

%%%%%%%%%%%%%% VISUALIZATION %%%%%%%%%%%%%%
subplot(2,2,1)
scatter(X_real(:,1),X_real(:,2))
xlabel('X_1')
ylabel('Y_1')
grid on
title('"Real" Joint Distribution')

subplot(2,2,2)
scatter(U_real(:,1),U_real(:,2))
xlabel('U_1')
ylabel('U_2')
grid on
title('"Real" Dependency Structure')

subplot(2,2,3);
scatter(X_gen(:,1),X_gen(:,2));
grid on
xlabel('X_1')
ylabel('X_2')
title('Copula Generative Model')

subplot(2,2,4);
scatter(U_gen(:,1),U_gen(:,2));
grid on
xlabel('X_1')
ylabel('X_2')
title('Generated Dependency Structure')

%% 3-D simulation
clear;
clc;
setupPath;

n = 1000;
D = 3;
[X_real, U_real] = gen_hybriddata_3D( n );

[ x1_domain, x1_multinomial_est, mle_params ] = fit_clg_3D( X_real );
[X1X2_clg_gen, X1X3_clg_gen, X2X3_clg_gen] = ...
    gen_samples_clg_3D(x1_domain, x1_multinomial_est, mle_params, n);

X1_continued = X_real(:,1) + (rand(n,1)-1);
X_transformed = [X1_continued X_real(:,2) X_real(:,3)];
K = 100;
D = size(X_real,2);

% [ U_gen, ~, U_emp ] = emp_copularnd_old( X_transformed, n, K );
[C,U,c] = empcopula(X_transformed,K);    
U_gen = empcopularnd(c, n);

X_gen = empdistrnd(U_gen, X_real);

% do an R-style "pairs" plot of the empirical data
subplot(3,3,1)
axis([0 1 0 1])
text(0.5,0.5,'X_1')

subplot(3,3,2)
scatter(X_real(:,1),X_real(:,2))
grid on
title('Empirical Data')

subplot(3,3,3)
scatter(X_real(:,1),X_real(:,3))
grid on

subplot(3,3,4)
scatter(X_real(:,2),X_real(:,1))
grid on

subplot(3,3,5)
axis([0 1 0 1])
text(0.5,0.5,'X_2')

subplot(3,3,6)
scatter(X_real(:,2),X_real(:,3))
grid on

subplot(3,3,7)
scatter(X_real(:,3),X_real(:,1))
grid on

subplot(3,3,8)
scatter(X_real(:,3),X_real(:,2))
grid on

subplot(3,3,9)
axis([0 1 0 1])
text(0.5,0.5,'X_3')

% do an R-style "pairs" plot of the CLG Generated data
figure;
subplot(3,3,1)
axis([0 1 0 1])
text(0.5,0.5,'X_1')

subplot(3,3,2)
scatter(X1X2_clg_gen(:,1),X1X2_clg_gen(:,2))
grid on
title('CLG Generative Model')

subplot(3,3,3)
scatter(X1X3_clg_gen(:,1),X1X3_clg_gen(:,2))
grid on

subplot(3,3,4)
scatter(X1X2_clg_gen(:,2),X1X2_clg_gen(:,1))
grid on

subplot(3,3,5)
axis([0 1 0 1])
text(0.5,0.5,'X_2')

subplot(3,3,6)
scatter(X2X3_clg_gen(:,1),X2X3_clg_gen(:,2))
grid on

subplot(3,3,7)
scatter(X1X3_clg_gen(:,2),X1X3_clg_gen(:,1))
grid on

subplot(3,3,8)
scatter(X2X3_clg_gen(:,2),X2X3_clg_gen(:,1))
grid on

subplot(3,3,9)
axis([0 1 0 1])
text(0.5,0.5,'X_3')

% do an R-style "pairs" plot of the Copula Generated data
figure;
subplot(3,3,1)
axis([0 1 0 1])
text(0.5,0.5,'X_1')

subplot(3,3,2)
scatter(X_gen(:,1),X_gen(:,2))
grid on
title('Copula Generative Data')

subplot(3,3,3)
scatter(X_gen(:,1),X_gen(:,3))
grid on

subplot(3,3,4)
scatter(X_gen(:,2),X_gen(:,1))
grid on

subplot(3,3,5)
axis([0 1 0 1])
text(0.5,0.5,'X_2')

subplot(3,3,6)
scatter(X_gen(:,2),X_gen(:,3))
grid on

subplot(3,3,7)
scatter(X_gen(:,3),X_gen(:,1))
grid on

subplot(3,3,8)
scatter(X_gen(:,3),X_gen(:,2))
grid on

subplot(3,3,9)
axis([0 1 0 1])
text(0.5,0.5,'X_3')

fprintf('Copulas WIN Again -- what did you expect?\n')