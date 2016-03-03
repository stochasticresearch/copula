%%
clear;
clc;
close all;

x = randn(1000,1);
y = randn(1000,1);
z = randn(1000,1);
x = 0.5*x + z;
y = 0.1*y + z.^2;

subplot(1,3,1); scatter(x,z);
subplot(1,3,2); scatter(y,z);
subplot(1,3,3); scatter(x,y);

%% Generate CI data w/ Z1, X, Y w/ some random dependencies
clear;
clc;
close all;

num_samps = 1000;
model = 'CI';
n = 1;

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

[ X, Y, Z ] = gen_conditionaldata( num_samps, model, n, ...
    distX, distY, distZ, fpre_XY, fcmb_XY, fpost_XY, ...
    fpre_Z, fcmb_Z, fpost_Z );

subplot(1,3,1); scatter(X,Z);
subplot(1,3,2); scatter(Y,Z);
subplot(1,3,3); scatter(X,Y);

%% Generate CI data w/ Z1, X, Y w/ some random dependencies
clear;
clc;
close all;

num_samps = 1000;
model = 'CI';
n = 1;

distX = makedist('Normal', 0, 1);
distY = makedist('Normal', 0, 1);
distZ = cell(n,1); distZ{1,1} = makedist('Normal', 0, 1);

fpre_XY = cell(1,2);  fpre_XY{1} = @tanh; fpre_XY{2} = @(x)power(x,2);
fcmb_XY = cell(1,2);  fcmb_XY{1} = @(x)opWeightedSum(x, repmat([0.5 0.5],num_samps,1));
                      fcmb_XY{2} = @(y)opWeightedSum(y, repmat([0.5 0.5],num_samps,1));
fpost_XY = cell(1,2); fpost_XY{1} = @opNo; fpost_XY{2} = @opNo;

fpre_Z = cell(n, 1); fpre_Z{1} = @(x)opAffine(x,2,4);
fcmb_Z = @opNo;
fpost_Z = @opNo;

[ X, Y, Z ] = gen_conditionaldata( num_samps, model, n, ...
    distX, distY, distZ, fpre_XY, fcmb_XY, fpost_XY, ...
    fpre_Z, fcmb_Z, fpost_Z );

%% Generate nCI data w/ Z1, X, Y