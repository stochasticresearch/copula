%% Generate CI data w/ Z, X, Y w/ some random dependencies
clear;
clc;
close all;

num_samps = 1000;
model = 'CI';
n = 1;

% The following configuration generates dependencies as follows:
%  Z = Normal(0,1);
%  X = 0.5*Normal(0,1) + Z;
%  Y = 0.1*Normal(0,1) + Z^2;

distX = makedist('Normal', 0, 1);
distY = makedist('Normal', 0, 1);
distZ = cell(n,1); distZ{1,1} = makedist('Normal', 0, 1);

fpre_XY = cell(1,2);  fpre_XY{1} = @opNo; fpre_XY{2} = @opNo;
fcmb_XY = cell(1,2);  fcmb_XY{1} = @(x)opWeightedSum(x, repmat([0.5 1], num_samps, 1));
                      fcmb_XY{2} = @(y)opWeightedSum(y, repmat([0.1 1], num_samps, 1));
fpost_XY = cell(1,2); fpost_XY{1} = @opNo; fpost_XY{2} = @opNo;

fpre_Z = cell(n, 2); fpre_Z{1,1} = @opNo; fpre_Z{1,2} = @(z)power(z,2);
fcmb_Z = cell(1, 2); fcmb_Z{1} = @opNo; fcmb_Z{2} = @opNo;
fpost_Z = cell(1,2); fpost_Z{1} = @opNo; fpost_Z{2} = @opNo;

[ Xt, Yt, Z, X, Y ] = gen_conditionaldata( num_samps, model, n, ...
    distX, distY, distZ, fpre_XY, fcmb_XY, fpost_XY, ...
    fpre_Z, fcmb_Z, fpost_Z );

subplot(3,3,1); ksdensity(X); title('f_X(x)'); grid on
subplot(3,3,2); ksdensity(Y); title('f_Y(y)'); grid on
subplot(3,3,3); ksdensity(Z); title('f_Z(z)'); grid on
subplot(3,3,4); ksdensity(Xt); title('f_{X_t}(x_t)'); grid on;
subplot(3,3,5); ksdensity(Yt); title('f_{Y_t}(y_t)'); grid on;
subplot(3,3,6); ksdensity(Z); title('f_Z(z)'); grid on;
subplot(3,3,7); scatter(Xt,Z); title('X_t = 0.5X + Z'); grid on;
subplot(3,3,8); scatter(Yt,Z); title('Y_t = 0.1Y + Z^2'); grid on;
subplot(3,3,9); scatter(Xt,Yt); title('X_t-Y_t'); grid on;

%% Generate CI data w/ Z1, Z2, X, Y w/ some random dependencies
clear;
clc;
close all;

num_samps = 1000;
model = 'CI';
n = 2;

% The following configuration generates dependencies as follows:
%  Z1 = Normal(0,1); Z2 = Exp(3);
%  X = log(0.4*sqrt[Normal(0,1)] + 0.6*exp(0.2*Z1 + 0.3*Z2^3) );
%  Y = sqrt(expm1[Normal(0,1)] * log(0.4*Z1^2 * 0.5*Z2) );
distX = makedist('Normal', 0, 1);
distY = makedist('Normal', 0, 1);
distZ = cell(n,1); distZ{1,1} = makedist('Normal', 0, 1); distZ{2,1} = makedist('Exponential', 3);

fpre_XY = cell(1,2);  fpre_XY{1} = @sqrt; fpre_XY{2} = @expm1;
fcmb_XY = cell(1,2);  fcmb_XY{1} = @(x)opWeightedSum(x, repmat([0.4 0.6], num_samps, 1));
                      fcmb_XY{2} = @(y)prod(y,2);       % element wise multiplication along a row
fpost_XY = cell(1,2); fpost_XY{1} = @log; fpost_XY{2} = @sqrt;

fpre_Z = cell(n, 2); fpre_Z{1,1} = @opNo; fpre_Z{2,1} = @(z2)power(z2,3);
                     fpre_Z{1,2} = @(z1)power(z1,2);  fpre_Z{2,2} = @opNo;
fcmb_Z = cell(1, 2); fcmb_Z{1} = @(x)opWeightedSum(x, repmat([0.2 0.3], num_samps, 1)); 
                     fcmb_Z{2} = @(y)opWeightedProduct(y, repmat([0.4 0.5], num_samps, 1));
fpost_Z = cell(1,2); fpost_Z{1} = @exp; fpost_Z{2} = @log;

[ Xt, Yt, Z, X, Y ] = gen_conditionaldata( num_samps, model, n, ...
    distX, distY, distZ, fpre_XY, fcmb_XY, fpost_XY, ...
    fpre_Z, fcmb_Z, fpost_Z );

subplot(4,4,1); ksdensity(X); title('f_X(x)'); grid on
subplot(4,4,2); ksdensity(Y); title('f_Y(y)'); grid on
subplot(4,4,3); ksdensity(Z(:,1)); title('f_Z_1(z)'); grid on
subplot(4,4,4); ksdensity(Z(:,2)); title('f_Z_2(z)'); grid on;

subplot(4,4,5); ksdensity(Xt); title('f_{X_t}(x_t)'); grid on;
subplot(4,4,6); ksdensity(Yt); title('f_{Y_t}(y_t)'); grid on;
subplot(4,4,7); ksdensity(Z(:,1)); title('f_Z_1(z)'); grid on;
subplot(4,4,8); ksdensity(Z(:,2)); title('f_Z_2(z)'); grid on;

subplot(4,4,9); scatter(Xt,Z(:,1)); title('X_t-Z_1'); grid on;
subplot(4,4,10); scatter(Yt,Z(:,1)); title('X_t-Z_1'); grid on;
subplot(4,4,[11 12 15 16]); scatter(Xt,Yt); title('X_t-Y_t'); grid on;

subplot(4,4,13); scatter(Xt,Z(:,2)); title('X_t-Z_2'); grid on;
subplot(4,4,14); scatter(Yt,Z(:,2)); title('X_t-Z_2'); grid on;

%% Generate CI data w/ Z1,Z2,Z3, X, Y w/ some random dependencies
clear;
clc;
close all;

num_samps = 1000;
model = 'CI';
n = 3;

% The following configuration generates dependencies as follows:
%  Z1 = Normal(0,1); Z2 = Normal(0,1); Z3 = Normal(0,1);
%  X = 0.5*Normal(0,1) + 2*(Z1 + Z2^2);
%  Y = 0.1*Normal(0,1) + 3*(Z1^2 - sqrt(Z2) + Z3^3);

distX = makedist('Normal', 0, 1);
distY = makedist('Normal', 0, 1);
distZ = cell(n,1); distZ{1,1} = makedist('Normal', 0, 1);
                   distZ{2,1} = makedist('Normal', 0, 1);
                   distZ{3,1} = makedist('Normal', 0, 1);

fpre_XY = cell(1,2);  fpre_XY{1} = @opNo; fpre_XY{2} = @opNo;
fcmb_XY = cell(1,2);  fcmb_XY{1} = @(x)opWeightedSum(x, repmat([0.5 2], num_samps, 1));
                      fcmb_XY{2} = @(y)opWeightedSum(y, repmat([0.1 3], num_samps, 1));
fpost_XY = cell(1,2); fpost_XY{1} = @opNo; fpost_XY{2} = @opNo;

fpre_Z = cell(n, 2); fpre_Z{1,1} = @opNo; fpre_Z{2,1} = @(z2x)power(z2x,2); fpre_Z{3,1} = @opNo; 
                     fpre_Z{1,2} = @(z1y)power(z1y,2); fpre_Z{2,2} = @sqrt; fpre_Z{3,2} = @(z3y)power(z3y,3);
fcmb_Z = cell(1, 2); fcmb_Z{1} = @(zx)opWeightedSum(zx, repmat([1 1 0], num_samps, 1));  % notice the 0 to disclude Z3
                     fcmb_Z{2} = @(zy)opWeightedSum(zy, repmat([1 -1 1], num_samps, 1));  % notice the 1 to include Z3
fpost_Z = cell(1,2); fpost_Z{1} = @opNo; fpost_Z{2} = @opNo;

[ Xt, Yt, Z, X, Y ] = gen_conditionaldata( num_samps, model, n, ...
    distX, distY, distZ, fpre_XY, fcmb_XY, fpost_XY, ...
    fpre_Z, fcmb_Z, fpost_Z );

subplot(2,4,1); scatter(Xt,Yt); title('X_t-Y_t'); grid on;
subplot(2,4,2); scatter(Xt,Z(:,1)); title('X_t-Z_1'); grid on;
subplot(2,4,3); scatter(Xt,Z(:,2)); title('X_t-Z_2'); grid on;
subplot(2,4,4); scatter(Xt,Z(:,3)); title('X_t-Z_3'); grid on;

subplot(2,4,5); scatter(Xt,Yt); title('X_t-Y_t'); grid on;
subplot(2,4,6); scatter(Yt,Z(:,1)); title('Y_t-Z_1'); grid on;
subplot(2,4,7); scatter(Yt,Z(:,2)); title('Y_t-Z_2'); grid on;
subplot(2,4,8); scatter(Yt,Z(:,3)); title('Y_t-Z_3'); grid on;
