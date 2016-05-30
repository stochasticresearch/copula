%**************************************************************************
%*                                                                        *
%* Copyright (C) 2016  Kiran Karra <kiran.karra@gmail.com>                *
%*                                                                        *
%* This program is free software: you can redistribute it and/or modify   *
%* it under the terms of the GNU General Public License as published by   *
%* the Free Software Foundation, either version 3 of the License, or      *
%* (at your option) any later version.                                    *
%*                                                                        *
%* This program is distributed in the hope that it will be useful,        *
%* but WITHOUT ANY WARRANTY; without even the implied warranty of         *
%* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          *
%* GNU General Public License for more details.                           *
%*                                                                        *
%* You should have received a copy of the GNU General Public License      *
%* along with this program.  If not, see <http://www.gnu.org/licenses/>.  *
%*                                                                        *
%**************************************************************************

% Tests the depcopularnd function

%% Test Gaussian Copula w/ depcopularnd
clear;
clc;

Rho1 = [1 0.5; 0.5 1];
Rho2 = [1 -0.3; -0.3 1];
M = 1000;

% generate the PDF
K = 100;
uu = linspace(0,1,K);
[U1,U2] = meshgrid(uu,uu);
c = copulapdf('Gaussian', [U1(:) U2(:)], Rho2);
c = reshape(c,K,K);

numMCsims = 100;
rho_XY_Z_hat_vec  = zeros(1,numMCsims);
rho_X_Y_hat_vec   = zeros(1,numMCsims);
rho_X_Z_hat_vec   = zeros(1,numMCsims);
rho_Y_Z_hat_vec   = zeros(1,numMCsims);
rho_Z1_Z2_hat_vec = zeros(1,numMCsims);

lastMsgLen = 0;
dispstat('','init'); % One time only initialization
dispstat(sprintf('Begining the simulation...'),'keepthis','timestamp');
for simnum=1:numMCsims
    msg = sprintf('Percentage Complete=%0.2f', simnum/numMCsims*100);
    dispstat(msg,'timestamp');
    
    
    U_init = copularnd('Gaussian', Rho1, M);        % generates [Z X]
    Z1 = U_init(:,1); X = U_init(:,2);
    U_dep = depcopularnd(c, Z1);
    Z2 = U_dep(:,1); Y = U_dep(:,2);
    
    rho_XY_Z_hat = partialcorr(X,Y,Z1); 
    rho_XY_Z_hat_vec(simnum) = rho_XY_Z_hat;
    
    rho_X_Y_hat  = corr(X,Y);           
    rho_X_Y_hat_vec(simnum) = rho_X_Y_hat;
    
    rho_X_Z_hat  = corr(X,Z1);          
    rho_X_Z_hat_vec(simnum) = rho_X_Z_hat;
    
    rho_Y_Z_hat  = corr(Y,Z1);          
    rho_Y_Z_hat_vec(simnum) = rho_Y_Z_hat;
    
    rho_Z1_Z2_hat = corr(Z1,Z2);        
    rho_Z1_Z2_hat_vec(simnum) = rho_Z1_Z2_hat;
end

fprintf(repmat('\b',1,lastMsgLen));
fprintf('********** Gaussian RV Testing - X indep Y | Z **********\n');
fprintf('rho(X,Y|Z) = %f\n', mean(rho_XY_Z_hat_vec) );
fprintf('rho(X,Y) = %f\n', mean(rho_X_Y_hat_vec) );
fprintf('rho(X,Z) = %f || \t\t ==> %f\n', mean(rho_X_Z_hat_vec), 0.5 );
fprintf('rho(Y,Z) = %f || \t ==> %f\n', mean(rho_Y_Z_hat_vec), -0.3 );
fprintf('rho(Z1,Z2)=%f\n', mean(rho_Z1_Z2_hat_vec) );

%% Test Clayton Copula w/ depcopularnd
clear;
clc;

alpha = 5;
M = 1000;
N = 2;

% generate the PDF
K = 100;
uu = linspace(0,1,K);
[U1,U2] = meshgrid(uu,uu);
c = claytoncopulapdf([U1(:) U2(:)], alpha);
c = reshape(c,K,K);

numMCsims = 100;
srho_XY_Z_hat_vec  = zeros(1,numMCsims);
srho_X_Y_hat_vec   = zeros(1,numMCsims);
srho_X_Z_hat_vec   = zeros(1,numMCsims);
srho_Y_Z_hat_vec   = zeros(1,numMCsims);
srho_Z1_Z2_hat_vec = zeros(1,numMCsims);

lastMsgLen = 0;
dispstat('','init'); % One time only initialization
dispstat(sprintf('Begining the simulation...'),'keepthis','timestamp');
for simnum=1:numMCsims
    msg = sprintf('Percentage Complete=%0.2f\n', simnum/numMCsims*100);
    dispstat(msg,'timestamp');

    
    [U_init] = claytoncopularnd(M, N, alpha);
    Z1 = U_init(:,1); X = U_init(:,2);
    U_dep = depcopularnd(c, Z1);
    Z2 = U_dep(:,1); Y = U_dep(:,2);    
    
    srho_XY_Z_hat = partialcorr(X,Y,Z1,'type','Spearman'); 
    srho_XY_Z_hat_vec(simnum) = srho_XY_Z_hat;
    
    srho_X_Y_hat  = corr(X,Y,'type','Spearman');           
    srho_X_Y_hat_vec(simnum) = srho_X_Y_hat;
    
    srho_X_Z_hat  = corr(X,Z1,'type','Spearman');          
    srho_X_Z_hat_vec(simnum) = srho_X_Z_hat;
    
    srho_Y_Z_hat  = corr(Y,Z1,'type','Spearman');          
    srho_Y_Z_hat_vec(simnum) = srho_Y_Z_hat;
    
    srho_Z1_Z2_hat = corr(Z1,Z2,'type','Spearman');        
    srho_Z1_Z2_hat_vec(simnum) = srho_Z1_Z2_hat;
end

fprintf(repmat('\b',1,lastMsgLen));
fprintf('********** Clayton Dependency Testing - X indep Y | Z **********\n');
fprintf('rho_s(X,Y|Z) = %f\n', mean(srho_XY_Z_hat_vec) );
fprintf('rho_s(X,Y) = %f\n', mean(srho_X_Y_hat_vec) );
fprintf('rho_s(X,Z) = %f || \t\t ==> %f\n', mean(srho_X_Z_hat_vec), copulastat('Clayton', alpha, 'type', 'Spearman') );
fprintf('rho_s(Y,Z) = %f || \t\t ==> %f\n', mean(srho_Y_Z_hat_vec), copulastat('Clayton', alpha, 'type', 'Spearman') );
fprintf('rho_s(Z1,Z2)=%f\n', mean(srho_Z1_Z2_hat_vec) );

%% Test Frank Copula w/ depcopularnd
clear;
clc;

alpha = 5;
M = 1000;
N = 2;

% generate the PDF
K = 100;
uu = linspace(0,1,K);
[U1,U2] = meshgrid(uu,uu);
c = frankcopulapdf([U1(:) U2(:)], alpha);
c = reshape(c,K,K);

numMCsims = 100;
srho_XY_Z_hat_vec  = zeros(1,numMCsims);
srho_X_Y_hat_vec   = zeros(1,numMCsims);
srho_X_Z_hat_vec   = zeros(1,numMCsims);
srho_Y_Z_hat_vec   = zeros(1,numMCsims);
srho_Z1_Z2_hat_vec = zeros(1,numMCsims);

lastMsgLen = 0;
dispstat('','init'); % One time only initialization
dispstat(sprintf('Begining the simulation...'),'keepthis','timestamp');
for simnum=1:numMCsims
    msg = sprintf('Percentage Complete=%0.02f\n', simnum/numMCsims*100);
    dispstat(msg,'timestamp');
    
    [U_init] = frankcopularnd(M, N, alpha);
    Z1 = U_init(:,1); X = U_init(:,2);
    U_dep = depcopularnd(c, Z1);
    Z2 = U_dep(:,1); Y = U_dep(:,2);

    srho_XY_Z_hat = partialcorr(X,Y,Z1,'type','Spearman'); 
    srho_XY_Z_hat_vec(simnum) = srho_XY_Z_hat;
    
    srho_X_Y_hat  = corr(X,Y,'type','Spearman');           
    srho_X_Y_hat_vec(simnum) = srho_X_Y_hat;
    
    srho_X_Z_hat  = corr(X,Z1,'type','Spearman');          
    srho_X_Z_hat_vec(simnum) = srho_X_Z_hat;
    
    srho_Y_Z_hat  = corr(Y,Z1,'type','Spearman');          
    srho_Y_Z_hat_vec(simnum) = srho_Y_Z_hat;
    
    srho_Z1_Z2_hat = corr(Z1,Z2,'type','Spearman');        
    srho_Z1_Z2_hat_vec(simnum) = srho_Z1_Z2_hat;
end

fprintf(repmat('\b',1,lastMsgLen));
fprintf('********** Frank Dependency Testing - X indep Y | Z **********\n');
fprintf('rho_s(X,Y|Z) = %f\n', mean(srho_XY_Z_hat_vec) );
fprintf('rho_s(X,Y) = %f\n', mean(srho_X_Y_hat_vec) );
fprintf('rho_s(X,Z) = %f || \t\t ==> %f\n', mean(srho_X_Z_hat_vec), copulastat('Frank', alpha, 'type', 'Spearman') );
fprintf('rho_s(Y,Z) = %f || \t\t ==> %f\n', mean(srho_Y_Z_hat_vec), copulastat('Frank', alpha, 'type', 'Spearman') );
fprintf('rho_s(Z1,Z2)=%f\n', mean(srho_Z1_Z2_hat_vec) );

%% Test Gumbel Copula w/ depcopularnd
clear;
clc;

alpha = 2;
M = 1000;
N = 2;

% generate the PDF
K = 100;
uu = linspace(0,1,K);
[U1,U2] = meshgrid(uu,uu);
c = gumbelcopulapdf([U1(:) U2(:)], alpha);
c = reshape(c,K,K);

numMCsims = 100;
srho_XY_Z_hat_vec  = zeros(1,numMCsims);
srho_X_Y_hat_vec   = zeros(1,numMCsims);
srho_X_Z_hat_vec   = zeros(1,numMCsims);
srho_Y_Z_hat_vec   = zeros(1,numMCsims);
srho_Z1_Z2_hat_vec = zeros(1,numMCsims);

lastMsgLen = 0;
dispstat('','init'); % One time only initialization
dispstat(sprintf('Begining the simulation...'),'keepthis','timestamp');
for simnum=1:numMCsims
    msg = sprintf('Percentage Complete=%0.02f\n', simnum/numMCsims*100);
    dispstat(msg,'timestamp');

    [U_init] = gumbelcopularnd(M, N, alpha);
    Z1 = U_init(:,1); X = U_init(:,2);
    U_dep = depcopularnd(c, Z1);
    Z2 = U_dep(:,1); Y = U_dep(:,2);
    
    srho_XY_Z_hat = partialcorr(X,Y,Z1,'type','Spearman'); 
    srho_XY_Z_hat_vec(simnum) = srho_XY_Z_hat;
    
    srho_X_Y_hat  = corr(X,Y,'type','Spearman');           
    srho_X_Y_hat_vec(simnum) = srho_X_Y_hat;
    
    srho_X_Z_hat  = corr(X,Z1,'type','Spearman');          
    srho_X_Z_hat_vec(simnum) = srho_X_Z_hat;
    
    srho_Y_Z_hat  = corr(Y,Z1,'type','Spearman');          
    srho_Y_Z_hat_vec(simnum) = srho_Y_Z_hat;
    
    srho_Z1_Z2_hat = corr(Z1,Z2,'type','Spearman');        
    srho_Z1_Z2_hat_vec(simnum) = srho_Z1_Z2_hat;
end

fprintf(repmat('\b',1,lastMsgLen));
fprintf('********** Gumbel Dependency Testing - X indep Y | Z **********\n');
fprintf('rho_s(X,Y|Z) = %f\n', mean(srho_XY_Z_hat_vec) );
fprintf('rho_s(X,Y) = %f\n', mean(srho_X_Y_hat_vec) );
fprintf('rho_s(X,Z) = %f || \t\t ==> %f\n', mean(srho_X_Z_hat_vec), copulastat('Gumbel', alpha, 'type', 'Spearman') );
fprintf('rho_s(Y,Z) = %f || \t\t ==> %f\n', mean(srho_Y_Z_hat_vec), copulastat('Gumbel', alpha, 'type', 'Spearman') );
fprintf('rho_s(Z1,Z2)=%f\n', mean(srho_Z1_Z2_hat_vec) );

%% Test 3-D Gaussian depcopularnd for the following model
%    Z1      Z2
%     |\    /| 
%     | \  / | 
%     |  \/  |
%     |  /\  |
%     | /  \ |
% X<---      --->Y
%  Arrows point FROM Z1 --> X, Y, and FROM Z2 --> X, Y

clear;
clc;

% define correlations
rho_Z1_Z2 = 0;
rho_Z1_X  = 0.2;
rho_Z1_Y  = 0.2;
rho_Z2_X  = -0.8;
rho_Z2_Y  = -0.6;

Rho1 = [1 rho_Z1_Z2 rho_Z1_X ; rho_Z1_Z2 1 rho_Z2_X; rho_Z1_X rho_Z2_X 1];
Rho2 = [1 rho_Z1_Z2 rho_Z1_Y ; rho_Z1_Z2 1 rho_Z2_Y; rho_Z1_Y rho_Z2_Y 1];

M = 1000;

numMCsims = 100;
rho_Z1_Z2_hat_vec = zeros(1,numMCsims);
rho_X_Z1_hat_vec = zeros(1,numMCsims);
rho_Y_Z1_hat_vec = zeros(1,numMCsims);
rho_X_Z2_hat_vec = zeros(1,numMCsims);
rho_Y_Z2_hat_vec = zeros(1,numMCsims);
rho_X_Y_hat_vec = zeros(1,numMCsims);
rho_X_Y_given_Z1Z2_hat_vec = zeros(1,numMCsims);
rho_X_Y_given_Z1_hat_vec = zeros(1,numMCsims);
rho_X_Y_given_Z2_hat_vec = zeros(1,numMCsims);

% generate the PDF
K = 100;
uu = linspace(0,1,K);
[U1,U2,U3] = ndgrid(uu);
c = copulapdf('Gaussian', [U1(:) U2(:) U3(:)], Rho2);
c = reshape(c,K,K,K);

lastMsgLen = 0;
dispstat('','init'); % One time only initialization
dispstat(sprintf('Begining the simulation...'),'keepthis','timestamp');
for simnum=1:numMCsims
    msg = sprintf('Percentage Complete=%0.02f\n', simnum/numMCsims*100);
    dispstat(msg,'timestamp');
    
    U_init = copularnd('Gaussian', Rho1, M);        % generates [Z1 Z2 X]
    
    Z1Z2_1 = U_init(:,1:2); X = U_init(:,3);
    U_dep = depcopularnd(c, Z1Z2_1);
    Z1Z2_2 = U_dep(:,1:2); Y = U_dep(:,3);
    
    % compute and store statistics
    rho_Z1_Z2_hat = corr(Z1Z2_1(:,1),Z1Z2_1(:,2));          rho_Z1_Z2_hat_vec(simnum) = rho_Z1_Z2_hat;
    rho_X_Z1_hat = corr(X,Z1Z2_1(:,1));                     rho_X_Z1_hat_vec(simnum) = rho_X_Z1_hat;
    rho_Y_Z1_hat = corr(Y,Z1Z2_1(:,1));                     rho_Y_Z1_hat_vec(simnum) = rho_Y_Z1_hat;
    rho_X_Z2_hat = corr(X,Z1Z2_1(:,2));                     rho_X_Z2_hat_vec(simnum) = rho_X_Z2_hat;
    rho_Y_Z2_hat = corr(Y,Z1Z2_1(:,2));                     rho_Y_Z2_hat_vec(simnum) = rho_Y_Z2_hat;
    rho_X_Y_hat  = corr(X,Y);                               rho_X_Y_hat_vec(simnum)  = rho_X_Y_hat;
    rho_X_Y_given_Z1Z2_hat = partialcorr(X,Y,Z1Z2_1);       rho_X_Y_given_Z1Z2_hat_vec(simnum) = rho_X_Y_given_Z1Z2_hat;
    rho_X_Y_given_Z1_hat   = partialcorr(X,Y,Z1Z2_1(:,1));  rho_X_Y_given_Z1_hat_vec(simnum) = rho_X_Y_given_Z1_hat;
    rho_X_Y_given_Z2_hat   = partialcorr(X,Y,Z1Z2_1(:,2));  rho_X_Y_given_Z2_hat_vec(simnum) = rho_X_Y_given_Z2_hat;

end

fprintf(repmat('\b',1,lastMsgLen));
fprintf('********** Gaussian D=3 Dependency Testing - X indep Y | {Z1,Z2} **********\n');
fprintf('rho(Z1,Z2)=%f \t\t\t||\t\t ==> 0\n', mean(rho_Z1_Z2_hat) );
fprintf('rho(X,Z1) = %f \t\t||\t\t ==> 0.2\n', mean(rho_X_Z1_hat) );
fprintf('rho(Y,Z1) = %f \t\t||\t\t ==> 0.2\n', mean(rho_Y_Z1_hat) );
fprintf('rho(X,Z2) = %f \t\t||\t\t ==> -0.8\n', mean(rho_X_Z2_hat) );
fprintf('rho(Y,Z2) = %f \t\t||\t\t ==> -0.6\n', mean(rho_Y_Z2_hat) );
fprintf('rho(X,Y) = %f\n', mean(rho_X_Y_hat) );
fprintf('rho(X,Y|Z1,Z2) = %f \t\t||\t\t ==> 0\n', mean(rho_X_Y_given_Z1Z2_hat_vec) );
fprintf('rho(X,Y|Z1) = %f\n', mean(rho_X_Y_given_Z1_hat_vec) );
fprintf('rho(X,Y|Z2) = %f\n', mean(rho_X_Y_given_Z2_hat_vec) );

%% Test 3-D Gaussian depcopularnd for the following model
%    Z1----->Z2
%     |\    /| 
%     | \  / | 
%     |  \/  |
%     |  /\  |
%     | /  \ |
% X<---      --->Y
%  Arrows point FROM Z1 --> X, Y, and FROM Z2 --> X, Y

clear;
clc;

% define correlations
rho_Z1_Z2 = 0.4;
rho_Z1_X  = 0.2;
rho_Z1_Y  = 0.2;
rho_Z2_X  = -0.8;
rho_Z2_Y  = -0.6;

Rho1 = [1 rho_Z1_Z2 rho_Z1_X ; rho_Z1_Z2 1 rho_Z2_X; rho_Z1_X rho_Z2_X 1];
Rho2 = [1 rho_Z1_Z2 rho_Z1_Y ; rho_Z1_Z2 1 rho_Z2_Y; rho_Z1_Y rho_Z2_Y 1];

M = 1000;

numMCsims = 100;
rho_Z1_Z2_hat_vec = zeros(1,numMCsims);
rho_X_Z1_hat_vec = zeros(1,numMCsims);
rho_Y_Z1_hat_vec = zeros(1,numMCsims);
rho_X_Z2_hat_vec = zeros(1,numMCsims);
rho_Y_Z2_hat_vec = zeros(1,numMCsims);
rho_X_Y_hat_vec = zeros(1,numMCsims);
rho_X_Y_given_Z1Z2_hat_vec = zeros(1,numMCsims);
rho_X_Y_given_Z1_hat_vec = zeros(1,numMCsims);
rho_X_Y_given_Z2_hat_vec = zeros(1,numMCsims);

% generate the PDF
K = 100;
uu = linspace(0,1,K);
[U1,U2,U3] = ndgrid(uu);
c = copulapdf('Gaussian', [U1(:) U2(:) U3(:)], Rho2);
c = reshape(c,K,K,K);

lastMsgLen = 0;
dispstat('','init'); % One time only initialization
dispstat(sprintf('Begining the simulation...'),'keepthis','timestamp');
for simnum=1:numMCsims
    msg = sprintf('Percentage Complete=%0.02f\n', simnum/numMCsims*100);
    dispstat(msg,'timestamp');
    
    U_init = copularnd('Gaussian', Rho1, M);        % generates [Z1 Z2 X]
    
    Z1Z2_1 = U_init(:,1:2); X = U_init(:,3);
    U_dep = depcopularnd(c, Z1Z2_1);
    Z1Z2_2 = U_dep(:,1:2); Y = U_dep(:,3);
    
    % compute and store statistics
    rho_Z1_Z2_hat = corr(Z1Z2_1(:,1),Z1Z2_1(:,2));          rho_Z1_Z2_hat_vec(simnum) = rho_Z1_Z2_hat;
    rho_X_Z1_hat = corr(X,Z1Z2_1(:,1));                     rho_X_Z1_hat_vec(simnum) = rho_X_Z1_hat;
    rho_Y_Z1_hat = corr(Y,Z1Z2_1(:,1));                     rho_Y_Z1_hat_vec(simnum) = rho_Y_Z1_hat;
    rho_X_Z2_hat = corr(X,Z1Z2_1(:,2));                     rho_X_Z2_hat_vec(simnum) = rho_X_Z2_hat;
    rho_Y_Z2_hat = corr(Y,Z1Z2_1(:,2));                     rho_Y_Z2_hat_vec(simnum) = rho_Y_Z2_hat;
    rho_X_Y_hat  = corr(X,Y);                               rho_X_Y_hat_vec(simnum)  = rho_X_Y_hat;
    rho_X_Y_given_Z1Z2_hat = partialcorr(X,Y,Z1Z2_1);       rho_X_Y_given_Z1Z2_hat_vec(simnum) = rho_X_Y_given_Z1Z2_hat;
    rho_X_Y_given_Z1_hat   = partialcorr(X,Y,Z1Z2_1(:,1));  rho_X_Y_given_Z1_hat_vec(simnum) = rho_X_Y_given_Z1_hat;
    rho_X_Y_given_Z2_hat   = partialcorr(X,Y,Z1Z2_1(:,2));  rho_X_Y_given_Z2_hat_vec(simnum) = rho_X_Y_given_Z2_hat;

end

fprintf(repmat('\b',1,lastMsgLen));
fprintf('********** Gaussian D=3 Dependency Testing - X indep Y | {Z1,Z2} **********\n');
fprintf('rho(Z1,Z2)=%f \t\t\t||\t\t ==> 0.4\n', mean(rho_Z1_Z2_hat) );
fprintf('rho(X,Z1) = %f \t\t||\t\t ==> 0.2\n', mean(rho_X_Z1_hat) );
fprintf('rho(Y,Z1) = %f \t\t||\t\t ==> 0.2\n', mean(rho_Y_Z1_hat) );
fprintf('rho(X,Z2) = %f \t\t||\t\t ==> -0.8\n', mean(rho_X_Z2_hat) );
fprintf('rho(Y,Z2) = %f \t\t||\t\t ==> -0.6\n', mean(rho_Y_Z2_hat) );
fprintf('rho(X,Y) = %f\n', mean(rho_X_Y_hat) );
fprintf('rho(X,Y|Z1,Z2) = %f \t\t||\t\t ==> 0\n', mean(rho_X_Y_given_Z1Z2_hat_vec) );
fprintf('rho(X,Y|Z1) = %f\n', mean(rho_X_Y_given_Z1_hat_vec) );
fprintf('rho(X,Y|Z2) = %f\n', mean(rho_X_Y_given_Z2_hat_vec) );

%% Test the 3-D Clayton by modeling data w/ the affine transform and compare
%    Z1----->Z2
%     |\    /| 
%     | \  / | 
%     |  \/  |
%     |  /\  |
%     | /  \ |
% X<---      --->Y
%  Arrows point FROM Z1 --> X, Y, and FROM Z2 --> X, Y

clear;
clc;

alpha = 5;  % remember, we can only have one alpha w/out nested archimedean copulas :(
M = 1000;
N = 3;

% generate the PDF
K = 100;
uu = linspace(0,1,K);
[U1,U2,U3] = ndgrid(uu);
c = claytoncopulapdf([U1(:) U2(:) U3(:)], alpha);
c = reshape(c,K,K,K);

numMCsims = 100;
srho_Z1_Z2_hat_vec = zeros(1,numMCsims);
srho_X_Z1_hat_vec = zeros(1,numMCsims);
srho_Y_Z1_hat_vec = zeros(1,numMCsims);
srho_X_Z2_hat_vec = zeros(1,numMCsims);
srho_Y_Z2_hat_vec = zeros(1,numMCsims);
srho_X_Y_hat_vec = zeros(1,numMCsims);
srho_X_Y_given_Z1Z2_hat_vec = zeros(1,numMCsims);
srho_X_Y_given_Z1_hat_vec = zeros(1,numMCsims);
srho_X_Y_given_Z2_hat_vec = zeros(1,numMCsims);

lastMsgLen = 0;
dispstat('','init'); % One time only initialization
dispstat(sprintf('Begining the simulation...'),'keepthis','timestamp');
for simnum=1:numMCsims
    msg = sprintf('Percentage Complete=%0.2f\n', simnum/numMCsims*100);
    dispstat(msg,'timestamp');
    
    [U_init] = claytoncopularnd(M, N, alpha);
    Z1Z2_1 = U_init(:,1:2); X = U_init(:,3);
    U_dep = depcopularnd(c, Z1Z2_1);
    Z1Z2_2 = U_dep(:,1:2); Y = U_dep(:,3);
    
    % compute and store statistics
    srho_Z1_Z2_hat = corr(Z1Z2_1(:,1),Z1Z2_1(:,2),'type','Spearman');          
    srho_Z1_Z2_hat_vec(simnum) = srho_Z1_Z2_hat;
    
    srho_X_Z1_hat = corr(X,Z1Z2_1(:,1),'type','Spearman');                     
    srho_X_Z1_hat_vec(simnum) = srho_X_Z1_hat;
    
    srho_Y_Z1_hat = corr(Y,Z1Z2_1(:,1),'type','Spearman');                     
    srho_Y_Z1_hat_vec(simnum) = srho_Y_Z1_hat;
    
    srho_X_Z2_hat = corr(X,Z1Z2_1(:,2),'type','Spearman');                     
    srho_X_Z2_hat_vec(simnum) = srho_X_Z2_hat;
    
    srho_Y_Z2_hat = corr(Y,Z1Z2_1(:,2),'type','Spearman');                     
    srho_Y_Z2_hat_vec(simnum) = srho_Y_Z2_hat;
    
    srho_X_Y_hat  = corr(X,Y,'type','Spearman');                               
    srho_X_Y_hat_vec(simnum)  = srho_X_Y_hat;
    
    srho_X_Y_given_Z1Z2_hat = partialcorr(X,Y,Z1Z2_1,'type','Spearman');       
    srho_X_Y_given_Z1Z2_hat_vec(simnum) = srho_X_Y_given_Z1Z2_hat;
    
    srho_X_Y_given_Z1_hat   = partialcorr(X,Y,Z1Z2_1(:,1),'type','Spearman');  
    srho_X_Y_given_Z1_hat_vec(simnum) = srho_X_Y_given_Z1_hat;
    
    srho_X_Y_given_Z2_hat   = partialcorr(X,Y,Z1Z2_1(:,2),'type','Spearman');  
    srho_X_Y_given_Z2_hat_vec(simnum) = srho_X_Y_given_Z2_hat;
end

fprintf(repmat('\b',1,lastMsgLen));
fprintf('********** Clayton D=3 Dependency Testing - X indep Y | {Z1,Z2} **********\n');
fprintf('rho(Z1,Z2)=%f \t\t||\t\t ==> %f\n', mean(srho_Z1_Z2_hat), copulastat('Clayton',alpha,'type','Spearman'));
fprintf('rho(X,Z1) = %f \t\t||\t\t ==> %f\n', mean(srho_X_Z1_hat), copulastat('Clayton',alpha,'type','Spearman'));
fprintf('rho(Y,Z1) = %f \t\t||\t\t ==> %f\n', mean(srho_Y_Z1_hat), copulastat('Clayton',alpha,'type','Spearman'));
fprintf('rho(X,Z2) = %f \t\t||\t\t ==> %f\n', mean(srho_X_Z2_hat), copulastat('Clayton',alpha,'type','Spearman'));
fprintf('rho(Y,Z2) = %f \t\t||\t\t ==> %f\n', mean(srho_Y_Z2_hat), copulastat('Clayton',alpha,'type','Spearman'));
fprintf('rho(X,Y) = %f \t\t||\t\t ==> %f\n', mean(srho_X_Y_hat), copulastat('Clayton',alpha,'type','Spearman'));
fprintf('rho(X,Y|Z1,Z2) = %f \t\t||\t\t ==> 0\n', mean(srho_X_Y_given_Z1Z2_hat_vec) );
fprintf('rho(X,Y|Z1) = %f\n', mean(srho_X_Y_given_Z1_hat_vec) );
fprintf('rho(X,Y|Z2) = %f\n', mean(srho_X_Y_given_Z2_hat_vec) );

%% Test the 3-D Frank by modeling data w/ the affine transform and compare
%    Z1----->Z2
%     |\    /| 
%     | \  / | 
%     |  \/  |
%     |  /\  |
%     | /  \ |
% X<---      --->Y
%  Arrows point FROM Z1 --> X, Y, and FROM Z2 --> X, Y

clear;
clc;

alpha = 5;  % remember, we can only have one alpha w/out nested archimedean copulas :(
M = 1000;
N = 3;

% generate the PDF
K = 100;
uu = linspace(0,1,K);
[U1,U2,U3] = ndgrid(uu);
c = frankcopulapdf([U1(:) U2(:) U3(:)], alpha);
c = reshape(c,K,K,K);

numMCsims = 100;
srho_Z1_Z2_hat_vec = zeros(1,numMCsims);
srho_X_Z1_hat_vec = zeros(1,numMCsims);
srho_Y_Z1_hat_vec = zeros(1,numMCsims);
srho_X_Z2_hat_vec = zeros(1,numMCsims);
srho_Y_Z2_hat_vec = zeros(1,numMCsims);
srho_X_Y_hat_vec = zeros(1,numMCsims);
srho_X_Y_given_Z1Z2_hat_vec = zeros(1,numMCsims);
srho_X_Y_given_Z1_hat_vec = zeros(1,numMCsims);
srho_X_Y_given_Z2_hat_vec = zeros(1,numMCsims);

lastMsgLen = 0;
dispstat('','init'); % One time only initialization
dispstat(sprintf('Begining the simulation...'),'keepthis','timestamp');
for simnum=1:numMCsims
    msg = sprintf('Percentage Complete=%0.2f\n', simnum/numMCsims*100);
    dispstat(msg,'timestamp');
    
    [U_init] = frankcopularnd(M, N, alpha);
    Z1Z2_1 = U_init(:,1:2); X = U_init(:,3);
    U_dep = depcopularnd(c, Z1Z2_1);
    Z1Z2_2 = U_dep(:,1:2); Y = U_dep(:,3);
    
    % compute and store statistics
    srho_Z1_Z2_hat = corr(Z1Z2_1(:,1),Z1Z2_1(:,2),'type','Spearman');          
    srho_Z1_Z2_hat_vec(simnum) = srho_Z1_Z2_hat;
    
    srho_X_Z1_hat = corr(X,Z1Z2_1(:,1),'type','Spearman');                     
    srho_X_Z1_hat_vec(simnum) = srho_X_Z1_hat;
    
    srho_Y_Z1_hat = corr(Y,Z1Z2_1(:,1),'type','Spearman');                     
    srho_Y_Z1_hat_vec(simnum) = srho_Y_Z1_hat;
    
    srho_X_Z2_hat = corr(X,Z1Z2_1(:,2),'type','Spearman');                     
    srho_X_Z2_hat_vec(simnum) = srho_X_Z2_hat;
    
    srho_Y_Z2_hat = corr(Y,Z1Z2_1(:,2),'type','Spearman');                     
    srho_Y_Z2_hat_vec(simnum) = srho_Y_Z2_hat;
    
    srho_X_Y_hat  = corr(X,Y,'type','Spearman');                               
    srho_X_Y_hat_vec(simnum)  = srho_X_Y_hat;
    
    srho_X_Y_given_Z1Z2_hat = partialcorr(X,Y,Z1Z2_1,'type','Spearman');       
    srho_X_Y_given_Z1Z2_hat_vec(simnum) = srho_X_Y_given_Z1Z2_hat;
    
    srho_X_Y_given_Z1_hat   = partialcorr(X,Y,Z1Z2_1(:,1),'type','Spearman');  
    srho_X_Y_given_Z1_hat_vec(simnum) = srho_X_Y_given_Z1_hat;
    
    srho_X_Y_given_Z2_hat   = partialcorr(X,Y,Z1Z2_1(:,2),'type','Spearman');  
    srho_X_Y_given_Z2_hat_vec(simnum) = srho_X_Y_given_Z2_hat;
end

fprintf(repmat('\b',1,lastMsgLen));
fprintf('********** Frank D=3 Dependency Testing - X indep Y | {Z1,Z2} **********\n');
fprintf('rho(Z1,Z2)=%f \t\t||\t\t ==> %f\n', mean(srho_Z1_Z2_hat), copulastat('Frank',alpha,'type','Spearman'));
fprintf('rho(X,Z1) = %f \t\t||\t\t ==> %f\n', mean(srho_X_Z1_hat), copulastat('Frank',alpha,'type','Spearman'));
fprintf('rho(Y,Z1) = %f \t\t||\t\t ==> %f\n', mean(srho_Y_Z1_hat), copulastat('Frank',alpha,'type','Spearman'));
fprintf('rho(X,Z2) = %f \t\t||\t\t ==> %f\n', mean(srho_X_Z2_hat), copulastat('Frank',alpha,'type','Spearman'));
fprintf('rho(Y,Z2) = %f \t\t||\t\t ==> %f\n', mean(srho_Y_Z2_hat), copulastat('Frank',alpha,'type','Spearman'));
fprintf('rho(X,Y) = %f \t\t||\t\t ==> %f\n', mean(srho_X_Y_hat), copulastat('Frank',alpha,'type','Spearman'));
fprintf('rho(X,Y|Z1,Z2) = %f \t\t||\t\t ==> 0\n', mean(srho_X_Y_given_Z1Z2_hat_vec) );
fprintf('rho(X,Y|Z1) = %f\n', mean(srho_X_Y_given_Z1_hat_vec) );
fprintf('rho(X,Y|Z2) = %f\n', mean(srho_X_Y_given_Z2_hat_vec) );

%% Test the 3-D Gumbel by modeling data w/ the affine transform and compare
%    Z1----->Z2
%     |\    /| 
%     | \  / | 
%     |  \/  |
%     |  /\  |
%     | /  \ |
% X<---      --->Y
%  Arrows point FROM Z1 --> X, Y, and FROM Z2 --> X, Y

clear;
clc;

alpha = 2;  % remember, we can only have one alpha w/out nested archimedean copulas :(
M = 1000;
N = 3;

% generate the PDF
K = 100;
uu = linspace(0,1,K);
[U1,U2,U3] = ndgrid(uu);
c = gumbelcopulapdf([U1(:) U2(:) U3(:)], alpha);
c = reshape(c,K,K,K);

numMCsims = 100;
srho_Z1_Z2_hat_vec = zeros(1,numMCsims);
srho_X_Z1_hat_vec = zeros(1,numMCsims);
srho_Y_Z1_hat_vec = zeros(1,numMCsims);
srho_X_Z2_hat_vec = zeros(1,numMCsims);
srho_Y_Z2_hat_vec = zeros(1,numMCsims);
srho_X_Y_hat_vec = zeros(1,numMCsims);
srho_X_Y_given_Z1Z2_hat_vec = zeros(1,numMCsims);
srho_X_Y_given_Z1_hat_vec = zeros(1,numMCsims);
srho_X_Y_given_Z2_hat_vec = zeros(1,numMCsims);

lastMsgLen = 0;
dispstat('','init'); % One time only initialization
dispstat(sprintf('Begining the simulation...'),'keepthis','timestamp');
for simnum=1:numMCsims
    msg = sprintf('Percentage Complete=%0.2f\n', simnum/numMCsims*100);
    dispstat(msg,'timestamp');
    
    [U_init] = gumbelcopularnd(M, N, alpha);
    Z1Z2_1 = U_init(:,1:2); X = U_init(:,3);
    U_dep = depcopularnd(c, Z1Z2_1);
    Z1Z2_2 = U_dep(:,1:2); Y = U_dep(:,3);
    
    % compute and store statistics
    srho_Z1_Z2_hat = corr(Z1Z2_1(:,1),Z1Z2_1(:,2),'type','Spearman');          
    srho_Z1_Z2_hat_vec(simnum) = srho_Z1_Z2_hat;
    
    srho_X_Z1_hat = corr(X,Z1Z2_1(:,1),'type','Spearman');                     
    srho_X_Z1_hat_vec(simnum) = srho_X_Z1_hat;
    
    srho_Y_Z1_hat = corr(Y,Z1Z2_1(:,1),'type','Spearman');                     
    srho_Y_Z1_hat_vec(simnum) = srho_Y_Z1_hat;
    
    srho_X_Z2_hat = corr(X,Z1Z2_1(:,2),'type','Spearman');                     
    srho_X_Z2_hat_vec(simnum) = srho_X_Z2_hat;
    
    srho_Y_Z2_hat = corr(Y,Z1Z2_1(:,2),'type','Spearman');                     
    srho_Y_Z2_hat_vec(simnum) = srho_Y_Z2_hat;
    
    srho_X_Y_hat  = corr(X,Y,'type','Spearman');                               
    srho_X_Y_hat_vec(simnum)  = srho_X_Y_hat;
    
    srho_X_Y_given_Z1Z2_hat = partialcorr(X,Y,Z1Z2_1,'type','Spearman');       
    srho_X_Y_given_Z1Z2_hat_vec(simnum) = srho_X_Y_given_Z1Z2_hat;
    
    srho_X_Y_given_Z1_hat   = partialcorr(X,Y,Z1Z2_1(:,1),'type','Spearman');  
    srho_X_Y_given_Z1_hat_vec(simnum) = srho_X_Y_given_Z1_hat;
    
    srho_X_Y_given_Z2_hat   = partialcorr(X,Y,Z1Z2_1(:,2),'type','Spearman');  
    srho_X_Y_given_Z2_hat_vec(simnum) = srho_X_Y_given_Z2_hat;
end

fprintf(repmat('\b',1,lastMsgLen));
fprintf('********** Gumbel D=3 Dependency Testing - X indep Y | {Z1,Z2} **********\n');
fprintf('rho(Z1,Z2)=%f \t\t||\t\t ==> %f\n', mean(srho_Z1_Z2_hat), copulastat('Gumbel',alpha,'type','Spearman'));
fprintf('rho(X,Z1) = %f \t\t||\t\t ==> %f\n', mean(srho_X_Z1_hat), copulastat('Gumbel',alpha,'type','Spearman'));
fprintf('rho(Y,Z1) = %f \t\t||\t\t ==> %f\n', mean(srho_Y_Z1_hat), copulastat('Gumbel',alpha,'type','Spearman'));
fprintf('rho(X,Z2) = %f \t\t||\t\t ==> %f\n', mean(srho_X_Z2_hat), copulastat('Gumbel',alpha,'type','Spearman'));
fprintf('rho(Y,Z2) = %f \t\t||\t\t ==> %f\n', mean(srho_Y_Z2_hat), copulastat('Gumbel',alpha,'type','Spearman'));
fprintf('rho(X,Y) = %f \t\t||\t\t ==> %f\n', mean(srho_X_Y_hat), copulastat('Gumbel',alpha,'type','Spearman'));
fprintf('rho(X,Y|Z1,Z2) = %f \t\t||\t\t ==> 0\n', mean(srho_X_Y_given_Z1Z2_hat_vec) );
fprintf('rho(X,Y|Z1) = %f\n', mean(srho_X_Y_given_Z1_hat_vec) );
fprintf('rho(X,Y|Z2) = %f\n', mean(srho_X_Y_given_Z2_hat_vec) );

%%
%% Tests w/ depcopularnd_old
%% 
%% Test Gaussian copula CI w/ depcopularnd_old
clear;
clc;

rng(12345);
Rho1 = [1 0.5; 0.5 1];
Rho2 = [1 -0.3; -0.3 1];
M = 1000;

U_init = copularnd('Gaussian', Rho1, M);        % generates [Z X]
U_dep = depcopularnd_old(U_init(:,1), 2, 'Gaussian', Rho2); % input [Z] to generate [Z Y]

Z1 = U_init(:,1); Z2 = U_dep(:,1);
X = U_init(:,2);
Y = U_dep(:,2);

fprintf('********** Gaussian RV Testing - X indep Y | Z **********\n');
fprintf('rho(X,Y|Z1) = %f\n', partialcorr(X,Y,Z1));
fprintf('rho(X,Y) = %f\n', corr(X,Y));
fprintf('rho(X,Z1) = %f\n', corr(X,Z1));
fprintf('rho(Y,Z1) = %f\n', corr(Y,Z1));

fprintf('rho(X,Y|Z2) = %f\n', partialcorr(X,Y,Z2));
fprintf('rho(X,Z2) = %f\n', corr(X,Z2));
fprintf('rho(Y,Z2) = %f\n', corr(Y,Z2));

fprintf('\nrho(Z1,Z2)=%f\n', corr(Z1,Z2));

Rho3 = [1 0.5 0.3; 0.5 1 -0.2; 0.3 -0.2 1];
U = mvnrnd([0 0 0], Rho3, M);
X = U(:,1); Y = U(:,2); Z1 = U(:,3);

fprintf('********** Gaussian RV Testing - X !indep Y | Z **********\n');
fprintf('rho(X,Y|Z) = %f\n', partialcorr(X,Y,Z1));
fprintf('rho(X,Y) = %f\n', corr(X,Y));
fprintf('rho(X,Z) = %f\n', corr(X,Z1));
fprintf('rho(Y,Z) = %f\n', corr(Y,Z1));

%% Test the Clayton Copula w/ depcopularnd_old
clear;
clc;

alpha = 5;
M = 1000;
N = 2;
[U_init] = claytoncopularnd(M, N, alpha);
U_dep = depcopularnd_old(U_init(:,1), N, 'Clayton', alpha+3);

Z1 = U_dep(:,1); Z2 = U_init(:,1);
X = U_init(:,2);
Y = U_dep(:,2);

fprintf('********** CLAYTON Dependency Testing - X indep Y | Z ******\n');
fprintf('rho_s(Z1,Z2)=%f\n', corr(Z1,Z2,'type','Spearman'));

fprintf('rho_s(X,Y|Z1)=%f\n', partialcorr(X, Y, Z1, 'Type', 'Spearman'));
fprintf('rho_s(X,Y)=%f\n', corr(X, Y, 'Type', 'Spearman'));
fprintf('rho_s(X,Z1)=%f\n', corr(X, Z1, 'Type', 'Spearman'));
fprintf('rho_s(Y,Z1)=%f\n', corr(Y, Z1, 'Type', 'Spearman'));

fprintf('rho_s(X,Y|Z2)=%f\n', partialcorr(X, Y, Z2, 'Type', 'Spearman'));
fprintf('rho_s(X,Z2)=%f\n', corr(X, Z2, 'Type', 'Spearman'));
fprintf('rho_s(X,Z2)=%f\n', corr(Y, Z2, 'Type', 'Spearman'));

%% Test the Frank Copula w/ depcopularnd_old
clear;
clc;

alpha = 5;
M = 1000;
N = 2;
[U_init] = frankcopularnd(M, N, alpha);
U_dep = depcopularnd_old(U_init(:,1), N, 'Frank', alpha-3);

Z1 = U_dep(:,1); Z2 = U_init(:,1);
X = U_init(:,2);
Y = U_dep(:,2);

fprintf('********** FRANK Dependency Testing - X indep Y | Z ******\n');
fprintf('rho_s(Z1,Z2)=%f\n', corr(Z1,Z2,'type','Spearman'));

fprintf('rho_s(X,Y|Z1)=%f\n', partialcorr(X, Y, Z1, 'Type', 'Spearman'));
fprintf('rho_s(X,Y)=%f\n', corr(X, Y, 'Type', 'Spearman'));
fprintf('rho_s(X,Z1)=%f\n', corr(X, Z1, 'Type', 'Spearman'));
fprintf('rho_s(Y,Z1)=%f\n', corr(Y, Z1, 'Type', 'Spearman'));

fprintf('rho_s(X,Y|Z2)=%f\n', partialcorr(X, Y, Z2, 'Type', 'Spearman'));
fprintf('rho_s(X,Z2)=%f\n', corr(X, Z2, 'Type', 'Spearman'));
fprintf('rho_s(X,Z2)=%f\n', corr(Y, Z2, 'Type', 'Spearman'));
fprintf('**************************************************************\n');