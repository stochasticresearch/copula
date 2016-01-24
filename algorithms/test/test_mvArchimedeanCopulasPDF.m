% Script which characterizes the error between the multivariate Archimedean
% copula PDF calculation and the results from R

clear;
clc;

alpha = 4;

% read input test file
u = csvread('../r_playground/testfiles/mvArchimedeanCopula_input.csv');
y_claytonR = csvread('../r_playground/testfiles/claytonPdf3D_output.csv');

y_claytonMlab = claytoncopulapdf(u, alpha);
mse_clayton_3D = sum((y_claytonR-y_claytonMlab).^2)