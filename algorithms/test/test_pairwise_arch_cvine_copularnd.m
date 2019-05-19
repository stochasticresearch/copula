% test the pairwise_arch_cvine_copularnd.m functionality
clear;
clc;

copula_type='Clayton';
d = 3;
% corrvec = linspace(0.7,0.95,d);
corrvec = [0.7, 0.8, 0.9];
num_samps = 500;
U = pairwise_arch_cvine_copularnd(copula_type,corrvec,num_samps);
plotmatrix(U)

%% test how we generate the gaussian version of the data
clear;
clc;
close all;

d = 4;
corrvec = linspace(0.15, 0.4, d);
RR = eye(d+1);
RR(d+1,1:d) = corrvec;
RR(1:d,d+1) = corrvec;
R = corrcov(nearcorr(RR));

U1 = copularnd('Gaussian', R, 500);
figure;
plotmatrix(U1)

df = 2;
U2 = copularnd('t', R, df, 500);
figure;
plotmatrix(U2)