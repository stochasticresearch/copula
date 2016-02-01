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

clear;
clc;
close all;

% Generate the data
% range = 1.5;
% n = 1000;
% xx = -range + 2*range*rand(n,1);       % generate uniform RV's
% yy = 2*xx + 1;
% x = [xx yy];

n = 1000;
rho = .7;
Z = mvnrnd([0 0], [1 rho; rho 1], n);
U_real = normcdf(Z);
gamma_a = 2;
gamma_b = 1;
x = [gaminv(U_real(:,1),gamma_a,gamma_b) tinv(U_real(:,2),5)];
K = 200;
M = size(x,1);
D = size(x,2);


[C,U,c] = empcopula_old(x,200); U1 = U(:,:,1); U2 = U(:,:,2);
% subplot(1,3,1); surf(U1,U2,C); grid on;

K = 200;
D = size(x,2);

[C_new,U_new,c_new] = empcopula(x,K);     

% This is incorrect checks, b/c before we misunderstood how to understand
% U1 and U2
% % make sure U_new and U are the same
% U1_new = U_new{1};
% U2_new = U_new{2};
% u1_check = isequal(U1,U1_new);
% u2_check = isequal(U2,U2_new);
% if(u1_check && u2_check)
%     fprintf('U checks passed!\n');
% else
%     fprintf('U checks failed!\n');
% end

figure;
subplot(1,3,1); surf(U1,U2,C);
subplot(1,3,2); surf(U1,U2,C_new);
subplot(1,3,3); surf(U1,U2,C-C_new);

ff = c_new{2};
figure;
subplot(1,3,1); surf(U1,U2,c)
subplot(1,3,2); surf(U1,U2,ff);
subplot(1,3,3); surf(U1,U2,c-ff);

% generate samples the old way
n_gen = 1000;
U_gen_old = emp_copularnd_old( x, n_gen, K );

% generate samples the new way
U_gen_new = empcopularnd(c_new, n_gen);

figure;
subplot(1,2,1); scatter(U_gen_old(:,1),U_gen_old(:,2))
subplot(1,2,2); scatter(U_gen_new(:,1),U_gen_new(:,2))

%% Test 3-D
clear;
clc;

n = 1000;
[X_real, U_real] = gen_hybriddata_3D( n );
x = X_real;
K = 200;
D = size(x,2);
n_gen = 1000;

% [C,U,c] = empcopula_old(x,K); U1 = U(:,:,1); U2 = U(:,:,2);

[C,U,c] = empcopula(x,K);     

% generate samples the old way
[U_gen_old,~,~,f] = emp_copularnd_old2( x, n_gen, K );

% generate samples the new way
U_gen_new = empcopularnd(c, n_gen);

figure;
subplot(1,3,1); scatter(U_real(:,1),U_real(:,2)); title('U_REAL');
subplot(1,3,2); scatter(U_gen_old(:,1),U_gen_old(:,2)); title('U_GEN_OLD');
subplot(1,3,3); scatter(U_gen_new(:,1),U_gen_new(:,2)); title('U_GEN_NEW');


%% 

[X1,X2] = meshgrid(linspace(-3,3,25)',linspace(-3,3,25)');
X = [X1(:) X2(:)];
p = mvncdf(X);
pp = reshape(p,25,25);
% subplot(1,2,1);
% surf(X1,X2,pp);

aa = diff(pp,1,1);
aa_tmp = diff(pp,1,2);
bb = diff(aa,1,2);
% subplot(1,2,2);
% surf(bb);

[cc] = gradient(pp);
[~,dd] = gradient(cc);

figure; 
subplot(1,3,1);
surf(bb);
title('Diff version')
subplot(1,3,2);
surf(dd)
title('Gradient version #1')

dd2 = gradient(cc');
subplot(1,3,3);
surf(dd2)
title('Gradient version #2')
