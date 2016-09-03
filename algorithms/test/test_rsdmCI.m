% Tests Conditional Independence w/ RSDM and the residuals processing
% algorithm

clear;
clc;

M = 500;
x = rand(M,1);
y = x.^2;
z = 2*y;

subplot(3,3,1);
scatter(x,y); grid on; xlabel('x'); ylabel('y');

subplot(3,3,2);
scatter(y,z); grid on; xlabel('y'); ylabel('z');

subplot(3,3,3);
scatter(x,z); grid on; xlabel('x'); ylabel('z');

subplot(3,3,4);
scatter(pobs(x),pobs(y)); grid on; xlabel('u'); ylabel('v');

subplot(3,3,5);
scatter(pobs(y),pobs(z)); grid on; xlabel('v'); ylabel('w');

subplot(3,3,6);
scatter(pobs(x),pobs(z)); grid on; xlabel('u'); ylabel('w');

rsdm10_rectWidths = [0.5 0.25 0.125];
rsdm10_overlapIncr = 0.005;
rsdm10_alpha = 0.01;
rsdm10_MA_size = 5;
rsdm10_peakFinderSel = 8;
rsdm10_boostFactor = 2;

[rsdm1, Rx, rect1] = copuladeptest10(x, y, 'kendall', rsdm10_rectWidths, ...
    rsdm10_overlapIncr, rsdm10_alpha, rsdm10_MA_size, rsdm10_peakFinderSel, ...
    rsdm10_boostFactor);
[rsdm2, Ry, rect2] = copuladeptest10(z, y, 'kendall', rsdm10_rectWidths, ...
    rsdm10_overlapIncr, rsdm10_alpha, rsdm10_MA_size, rsdm10_peakFinderSel, ...
    rsdm10_boostFactor);

subplot(3,3,7);
scatter(Rx,Ry); grid on; xlabel('xx'); ylabel('yy');

subplot(3,3,8);
scatter(pobs(Rx),pobs(Ry)); grid on; xlabel('uu'); ylabel('vv');