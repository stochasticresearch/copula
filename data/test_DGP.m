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

%% A script to plot the linear DGP and see what the dependencies look like
clear;
clc;

rng(1234);

M = 500;
W1_prev = [0 0];
W4_prev = 0;
W5_prev = [0 0];

[W1, W2, W3, W4, W5] = linearDGP(M, W1_prev, W4_prev, W5_prev);

figure;
subplot(4,2,1);
x = W1(1:end-1); y = W1(2:end);
scatter(x,y); 
grid on; xlabel('W_1(t)'); ylabel('W_1(t-1)');
title(sprintf('RSDM=%0.02f\n', rsdm(x,y)));

subplot(4,2,2);
x = W1(1:end-2); y = W1(3:end);
scatter(x, y); 
grid on; xlabel('W_1(t)'); ylabel('W_1(t-2)');
title(sprintf('RSDM=%0.02f\n', rsdm(x,y)));

subplot(4,2,3);
x = W2(1:end-2); y = W1(3:end);
scatter(x,y); 
grid on; xlabel('W_2(t)'); ylabel('W_1(t-2)');
title(sprintf('RSDM=%0.02f\n', rsdm(x,y)));

subplot(4,2,4);
x = W3(1:end-2); y = W1(3:end);
scatter(x,y); 
grid on; xlabel('W_3(t)'); ylabel('W_1(t-2)');
title(sprintf('RSDM=%0.02f\n', rsdm(x,y)));

subplot(4,2,5);
x = W4(1:end-1); y = W1(2:end);
scatter(x,y); 
grid on; xlabel('W_4(t)'); ylabel('W_1(t-1)');
title(sprintf('RSDM=%0.02f\n', rsdm(x,y)));

subplot(4,2,6);
x = W4(1:end-1); y = W4(2:end);
scatter(x,y); 
grid on; xlabel('W_4(t)'); ylabel('W_4(t-1)');
title(sprintf('RSDM=%0.02f\n', rsdm(x,y)));

subplot(4,2,7);
x = W4(1:end-1); y = W5(2:end);
scatter(x,y); 
grid on; xlabel('W_4(t)'); ylabel('W_5(t-1)');
title(sprintf('RSDM=%0.02f\n', rsdm(x,y)));

subplot(4,2,8);
x = W5(1:end-1); y = W4(2:end);
scatter(x,y); 
grid on; xlabel('W_5(t)'); ylabel('W_4(t-1)');
title(sprintf('RSDM=%0.02f\n', rsdm(x,y)));

% Explore RSCDM
% See if W1 has an effect on W2
figure;
x = W1(1:end-1); y = W1(2:end); z = W2(2:end);
[rscdmVal, RxStacked, RyStacked] = rscdm(x,y,z);

subplot(2,2,1); scatter(x,y); 
grid on; xlabel('W_1(t)'); ylabel('W_1(t-1)');
title(sprintf('RSDM=%0.02f\n', rsdm(x,y)));

subplot(2,2,4); scatter(x,z); 
grid on; xlabel('W_1(t)'); ylabel('W_2(t-1)');
title(sprintf('RSDM=%0.02f\n', rsdm(x,z)));

subplot(2,2,3); scatter(y,z); 
grid on; xlabel('W_1(t-1)'); ylabel('W_2(t-1)');
title(sprintf('RSDM=%0.02f\n', rsdm(y,z)));

subplot(2,2,2); scatter(RxStacked, RyStacked); 
grid on; xlabel('W_1(t)|W_2(t-1)'); ylabel('W_1(t-1)|W_2(t-1)');
title(sprintf('RSCDM=%0.02f\n', rscdmVal));

% now see if W1 has an effect on W2 (it should, looking at the DGP)

%% A script to plot the nonlinear DGP and see what the dependencies look like
clear;
clc;

rng(1234);

M = 500;
W1_prev = 0;
W2_prev = 0;

[W1, W2, W3, W4] = nonlinearDGP(M, W1_prev, W2_prev);

subplot(4,2,1);
x = W1(1:end-1); y = W1(2:end);
scatter(x,y);
grid on; xlabel('W_1(t)'); ylabel('W_1(t-1)');
title(sprintf('RSDM=%0.02f\n', rsdm(x,y)));

subplot(4,2,2);
x = W1(1:end-1); y = W2(2:end);
scatter(x,y);
grid on; xlabel('W_1(t)'); ylabel('W_2(t-1)');
title(sprintf('RSDM=%0.02f\n', rsdm(x,y)));

subplot(4,2,3);
x = W2(1:end-1); y = W1(2:end);
scatter(x,y);
grid on; xlabel('W_2(t)'); ylabel('W_1(t-1)');
title(sprintf('RSDM=%0.02f\n', rsdm(x,y)));

subplot(4,2,4);
x = W2(1:end-1); y = W2(2:end);
scatter(x,y);
grid on; xlabel('W_2(t)'); ylabel('W_2(t-1)');
title(sprintf('RSDM=%0.02f\n', rsdm(x,y)));

subplot(4,2,5);
x = W2(1:end-1); y = W4(2:end);
scatter(x,y);
grid on; xlabel('W_2(t)'); ylabel('W_4(t-1)');
title(sprintf('RSDM=%0.02f\n', rsdm(x,y)));

subplot(4,2,6);
x = W3(1:end-1); y = W1(2:end);
scatter(x,y);
grid on; xlabel('W_3(t)'); ylabel('W_1(t-1)');
title(sprintf('RSDM=%0.02f\n', rsdm(x,y)));

subplot(4,2,7);
x = W3(1:end-1); y = W2(2:end);
scatter(x,y);
grid on; xlabel('W_3(t)'); ylabel('W_2(t-1)');
title(sprintf('RSDM=%0.02f\n', rsdm(x,y)));
