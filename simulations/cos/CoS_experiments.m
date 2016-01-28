%******************************************************************************
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

% A script which explores functional dependence, and how the empirical
% copula changes w/ different functional dependencies

clear;
clc;

% Generate the data
range = 1.5;
n = 1000;
x = -range + 2*range*rand(n,1);       % generate uniform RV's
% x = rand(1000,1);       % generate uniform RV's
% transform to desired distribution

% x = norminv(x, 0, 1);

% x = tinv(x, 3);

% x = gaminv(x,2,1);

% alpha = 5;
% beta = 1;
% betaDist = makedist('Beta', alpha, beta);
% x = icdf(betaDist, x);

figtitlestr = sprintf('X ~ Uniform(-%d,%d)',range,range);
% figtitlestr = 'X ~ Beta(5,1)';
neg = 0;

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


N = 500;
K = 100;

subplot_idx = 1;
% affine dependence
fprintf('Processing Affine Dependence\n');
y = 2*x + 1;
titStr = 'Y=2X+1';
if(neg)
    y = y*-1;
    titStr = 'Y=-1*(2X+1)';
end
subplot(4,4,subplot_idx); subplot_idx = subplot_idx + 1;
scatter(x, y, 'r')
title(titStr);
xlabel('X')
ylabel('Y')
grid on

subplot(4,4,subplot_idx); subplot_idx = subplot_idx + 1;
[ecop,U,c] = empcopula([x y],K); U1 = U{1}; U2 = U{2};
contour(U1,U2,ecop)
xlabel('U_1')
ylabel('U_2')
title('$$\hat{C}(F_X,G_Y)$$', 'Interpreter', 'Latex');

U_gen = empcopularnd(c, N);
subplot(4,4,subplot_idx); subplot_idx = subplot_idx + 1;
scatter(U_gen(:,1), U_gen(:,2))
title('Genearted Samples from $$\hat{C}(F_X,G_Y)$$', 'Interpreter', 'Latex');
xlabel('U_1')
ylabel('U_2')
grid on

subplot(4,4,subplot_idx); subplot_idx = subplot_idx + 1;
hist3([x y], [100 100]);
set(get(gca,'child'),'FaceColor','interp','CDataMode','auto');
title('f_X_Y(X,Y)');
xlabel('U_1')
ylabel('U_2')
grid on


% quadratic dependence
fprintf('Processing Quadratic Dependence\n');
y = x.^2;
titStr = 'Y=X^2';
if(neg)
    y = y*-1;
    titStr = 'Y=-X^2';
end
subplot(4,4,subplot_idx); subplot_idx = subplot_idx + 1;
scatter(x, y, 'r')
title(titStr);
xlabel('X')
ylabel('Y')
grid on

subplot(4,4,subplot_idx); subplot_idx = subplot_idx + 1;
[ecop,U,c] = empcopula([x y],K); U1 = U{1}; U2 = U{2};
contour(U1,U2,ecop)
xlabel('U_1')
ylabel('U_2')
title('$$\hat{C}(F_X,G_Y)$$', 'Interpreter', 'Latex');

U_gen = empcopularnd(c, N);
subplot(4,4,subplot_idx); subplot_idx = subplot_idx + 1;
scatter(U_gen(:,1), U_gen(:,2))
title('Samples of $$\hat{C}(F_X,G_Y)$$', 'Interpreter', 'Latex');
xlabel('U_1')
ylabel('U_2')
grid on

subplot(4,4,subplot_idx); subplot_idx = subplot_idx + 1;
hist3([x y], [100 100]);
set(get(gca,'child'),'FaceColor','interp','CDataMode','auto');
title('f_X_Y(X,Y)');
xlabel('U_1')
ylabel('U_2')
grid on

% sinusoidal dependence
fprintf('Processing Sinusoidal Dependence\n');
y = sin(x);
titStr = 'Y=sin(X)';
if(neg)
    y = y*-1;
    titStr = 'Y= -sin(X)';
end
subplot(4,4,subplot_idx); subplot_idx = subplot_idx + 1;
scatter(x, y, 'r')
title(titStr);
xlabel('X')
ylabel('Y')
grid on

subplot(4,4,subplot_idx); subplot_idx = subplot_idx + 1;
[ecop,U,c] = empcopula([x y],K); U1 = U{1}; U2 = U{2};
contour(U1,U2,ecop)
xlabel('U_1')
ylabel('U_2')
title('$$\hat{C}(F_X,G_Y)$$', 'Interpreter', 'Latex');

U_gen = empcopularnd(c, N);
subplot(4,4,subplot_idx); subplot_idx = subplot_idx + 1;
scatter(U_gen(:,1), U_gen(:,2))
title('Samples of $$\hat{C}(F_X,G_Y)$$', 'Interpreter', 'Latex');
xlabel('U_1')
ylabel('U_2')
grid on

subplot(4,4,subplot_idx); subplot_idx = subplot_idx + 1;
hist3([x y], [100 100]);
set(get(gca,'child'),'FaceColor','interp','CDataMode','auto');
title('f_X_Y(X,Y)');
xlabel('U_1')
ylabel('U_2')
grid on

% 4th order Poly dependence
fprintf('Processing 4th order polynomial dependence\n');
y = x.^4-2*x.^2;
titStr = 'Y=X^4-2X^2';
if(neg)
    y = y*-1;
    titStr = 'Y= -(X^4-2X^2)';
end
subplot(4,4,subplot_idx); subplot_idx = subplot_idx + 1;
scatter(x, y, 'r')
title(titStr);
xlabel('X')
ylabel('Y')
grid on

subplot(4,4,subplot_idx); subplot_idx = subplot_idx + 1;
[ecop,U,c] = empcopula([x y],K); U1 = U{1}; U2 = U{2};
contour(U1,U2,ecop)
xlabel('U_1')
ylabel('U_2')
title('$$\hat{C}(F_X,G_Y)$$', 'Interpreter', 'Latex');

U_gen = empcopularnd(c, N);
subplot(4,4,subplot_idx); subplot_idx = subplot_idx + 1;
scatter(U_gen(:,1), U_gen(:,2))
title('Samples of $$\hat{C}(F_X,G_Y)$$', 'Interpreter', 'Latex');
xlabel('U_1')
ylabel('U_2')
grid on

subplot(4,4,subplot_idx); subplot_idx = subplot_idx + 1;
hist3([x y], [100 100]);
set(get(gca,'child'),'FaceColor','interp','CDataMode','auto');
title('f_X_Y(X,Y)');
xlabel('U_1')
ylabel('U_2')
grid on

h = figtitle(figtitlestr);

%%  Experiment w/ Polynomial functions only to understand effect of slope
% quadratic dependence
y = x.^2;
titStr = 'Y=X^2';
if(neg)
    y = y*-1;
    titStr = 'Y=-X^2';
end
subplot_idx = 1;
subplot(2,4,subplot_idx); subplot_idx = subplot_idx + 1;
scatter(x, y, 'r')
title(titStr);
xlabel('X')
ylabel('Y')
grid on

subplot(2,4,subplot_idx); subplot_idx = subplot_idx + 1;
[ecop,U,c] = empcopula([x y],K); U1 = U{1}; U2 = U{2};
contour(U1,U2,ecop)
xlabel('U_1')
ylabel('U_2')
title('$$\hat{C}(F_X,G_Y)$$', 'Interpreter', 'Latex');

% quadratic dependence
y = x.^3;
titStr = 'Y=X^3';
if(neg)
    y = y*-1;
    titStr = 'Y=-X^3';
end
subplot(2,4,subplot_idx); subplot_idx = subplot_idx + 1;
scatter(x, y, 'r')
title(titStr);
xlabel('X')
ylabel('Y')
grid on

subplot(2,4,subplot_idx); subplot_idx = subplot_idx + 1;
[ecop,U,c] = empcopula([x y],K); U1 = U{1}; U2 = U{2};
contour(U1,U2,ecop)
xlabel('U_1')
ylabel('U_2')
title('$$\hat{C}(F_X,G_Y)$$', 'Interpreter', 'Latex');

% quadratic dependence
y = x.^4;
titStr = 'Y=X^4';
if(neg)
    y = y*-1;
    titStr = 'Y=-X^4';
end
subplot(2,4,subplot_idx); subplot_idx = subplot_idx + 1;
scatter(x, y, 'r')
title(titStr);
xlabel('X')
ylabel('Y')
grid on

subplot(2,4,subplot_idx); subplot_idx = subplot_idx + 1;
[ecop,U,c] = empcopula([x y],K); U1 = U{1}; U2 = U{2};
contour(U1,U2,ecop)
xlabel('U_1')
ylabel('U_2')
title('$$\hat{C}(F_X,G_Y)$$', 'Interpreter', 'Latex');

% quadratic dependence
y = x.^4-2*x.^2;
titStr = 'Y=X^4-2X^2';
if(neg)
    y = y*-1;
    titStr = 'Y=-(X^4-2X^2)';
end
subplot(2,4,subplot_idx); subplot_idx = subplot_idx + 1;
scatter(x, y, 'r')
title(titStr);
xlabel('X')
ylabel('Y')
grid on

subplot(2,4,subplot_idx); subplot_idx = subplot_idx + 1;
[ecop, U1, U2] = ecopula([x y]);
contour(U1,U2,ecop)
xlabel('U_1')
ylabel('U_2')
title('$$\hat{C}(F_X,G_Y)$$', 'Interpreter', 'Latex');


h = figtitle(figtitlestr);

%% Polynomial Experiments
figure;

y = -10*x.^2;
subplot_idx = 1;
subplot(2,3,subplot_idx); subplot_idx = subplot_idx + 1;
[ecop,U,c] = empcopula([x y],K); U1 = U{1}; U2 = U{2};
contour(U1,U2,ecop)
xlabel('U_1')
ylabel('U_2')
title('$$\hat{C}(F_X,G_Y) | Y = -10X^2$$', 'Interpreter', 'Latex');

y = -5*x.^2;
subplot(2,3,subplot_idx); subplot_idx = subplot_idx + 1;
[ecop,U,c] = empcopula([x y],K); U1 = U{1}; U2 = U{2};
contour(U1,U2,ecop)
xlabel('U_1')
ylabel('U_2')
title('$$\hat{C}(F_X,G_Y) | Y = -5X^2$$', 'Interpreter', 'Latex');

y = -x.^2;
subplot(2,3,subplot_idx); subplot_idx = subplot_idx + 1;
[ecop,U,c] = empcopula([x y],K); U1 = U{1}; U2 = U{2};
contour(U1,U2,ecop)
xlabel('U_1')
ylabel('U_2')
title('$$\hat{C}(F_X,G_Y) | Y = -X^2$$', 'Interpreter', 'Latex');

y = x.^2;
subplot(2,3,subplot_idx); subplot_idx = subplot_idx + 1;
[ecop,U,c] = empcopula([x y],K); U1 = U{1}; U2 = U{2};
contour(U1,U2,ecop)
xlabel('U_1')
ylabel('U_2')
title('$$\hat{C}(F_X,G_Y) | Y = X^2$$', 'Interpreter', 'Latex');

y = 5*x.^2;
subplot(2,3,subplot_idx); subplot_idx = subplot_idx + 1;
[ecop,U,c] = empcopula([x y],K); U1 = U{1}; U2 = U{2};
contour(U1,U2,ecop)
xlabel('U_1')
ylabel('U_2')
title('$$\hat{C}(F_X,G_Y) | Y = 5X^2$$', 'Interpreter', 'Latex');

y = 10*x.^2;
subplot(2,3,subplot_idx); subplot_idx = subplot_idx + 1;
[ecop,U,c] = empcopula([x y],K); U1 = U{1}; U2 = U{2};
contour(U1,U2,ecop)
xlabel('U_1')
ylabel('U_2')
title('$$\hat{C}(F_X,G_Y) | Y = 10X^2$$', 'Interpreter', 'Latex');

h = figtitle(figtitlestr);

figure;
y = -(x-5).^2;
subplot_idx = 1;
subplot(2,3,subplot_idx); subplot_idx = subplot_idx + 1;
[ecop,U,c] = empcopula([x y],K); U1 = U{1}; U2 = U{2};
contour(U1,U2,ecop)
xlabel('U_1')
ylabel('U_2')
title('$$\hat{C}(F_X,G_Y) | Y = -(X-5)^2$$', 'Interpreter', 'Latex');

y = -x.^2;
subplot(2,3,subplot_idx); subplot_idx = subplot_idx + 1;
[ecop,U,c] = empcopula([x y],K); U1 = U{1}; U2 = U{2};
contour(U1,U2,ecop)
xlabel('U_1')
ylabel('U_2')
title('$$\hat{C}(F_X,G_Y) | Y = -X^2$$', 'Interpreter', 'Latex');

y = -(x+5).^2;
subplot(2,3,subplot_idx); subplot_idx = subplot_idx + 1;
[ecop,U,c] = empcopula([x y],K); U1 = U{1}; U2 = U{2};
contour(U1,U2,ecop)
xlabel('U_1')
ylabel('U_2')
title('$$\hat{C}(F_X,G_Y) | Y = -(X+5)^2$$', 'Interpreter', 'Latex');

y = (x-5).^2;
subplot(2,3,subplot_idx); subplot_idx = subplot_idx + 1;
[ecop,U,c] = empcopula([x y],K); U1 = U{1}; U2 = U{2};
contour(U1,U2,ecop)
xlabel('U_1')
ylabel('U_2')
title('$$\hat{C}(F_X,G_Y) | Y = (X-5)^2$$', 'Interpreter', 'Latex');

y = x.^2;
subplot(2,3,subplot_idx); subplot_idx = subplot_idx + 1;
[ecop,U,c] = empcopula([x y],K); U1 = U{1}; U2 = U{2};
contour(U1,U2,ecop)
xlabel('U_1')
ylabel('U_2')
title('$$\hat{C}(F_X,G_Y) | Y = X^2$$', 'Interpreter', 'Latex');

y = (x+5).^2;
subplot(2,3,subplot_idx); subplot_idx = subplot_idx + 1;
[ecop,U,c] = empcopula([x y],K); U1 = U{1}; U2 = U{2};
contour(U1,U2,ecop)
xlabel('U_1')
ylabel('U_2')
title('$$\hat{C}(F_X,G_Y) | Y = (X+5)^2$$', 'Interpreter', 'Latex');

h = figtitle(figtitlestr);

%% sinusoidal experiments
figure;

y = sin(-10*x);
subplot_idx = 1;
subplot(2,3,subplot_idx); subplot_idx = subplot_idx + 1;
[ecop,U,c] = empcopula([x y],K); U1 = U{1}; U2 = U{2};
contour(U1,U2,ecop)
xlabel('U_1')
ylabel('U_2')
title('$$\hat{C}(F_X,G_Y) | Y = sin(-10X)$$', 'Interpreter', 'Latex');

y = sin(-5*x);
subplot(2,3,subplot_idx); subplot_idx = subplot_idx + 1;
[ecop,U,c] = empcopula([x y],K); U1 = U{1}; U2 = U{2};
contour(U1,U2,ecop)
xlabel('U_1')
ylabel('U_2')
title('$$\hat{C}(F_X,G_Y) | Y = sin(-5X)$$', 'Interpreter', 'Latex');

y = sin(-x);
subplot(2,3,subplot_idx); subplot_idx = subplot_idx + 1;
[ecop,U,c] = empcopula([x y],K); U1 = U{1}; U2 = U{2};
contour(U1,U2,ecop)
xlabel('U_1')
ylabel('U_2')
title('$$\hat{C}(F_X,G_Y) | Y = sin(-X)$$', 'Interpreter', 'Latex');

y = sin(x);
subplot(2,3,subplot_idx); subplot_idx = subplot_idx + 1;
[ecop,U,c] = empcopula([x y],K); U1 = U{1}; U2 = U{2};
contour(U1,U2,ecop)
xlabel('U_1')
ylabel('U_2')
title('$$\hat{C}(F_X,G_Y) | Y = sin(X)$$', 'Interpreter', 'Latex');

y = sin(5*x);
subplot(2,3,subplot_idx); subplot_idx = subplot_idx + 1;
[ecop,U,c] = empcopula([x y],K); U1 = U{1}; U2 = U{2};
contour(U1,U2,ecop)
xlabel('U_1')
ylabel('U_2')
title('$$\hat{C}(F_X,G_Y) | Y = sin(5X)$$', 'Interpreter', 'Latex');

y = sin(10*x);
subplot(2,3,subplot_idx); subplot_idx = subplot_idx + 1;
[ecop,U,c] = empcopula([x y],K); U1 = U{1}; U2 = U{2};
contour(U1,U2,ecop)
xlabel('U_1')
ylabel('U_2')
title('$$\hat{C}(F_X,G_Y) | Y = sin(10X)$$', 'Interpreter', 'Latex');

h = figtitle(figtitlestr);

figure;
y = sin(x-10);
subplot_idx = 1;
subplot(2,2,subplot_idx); subplot_idx = subplot_idx + 1;
[ecop,U,c] = empcopula([x y],K); U1 = U{1}; U2 = U{2};
contour(U1,U2,ecop)
xlabel('U_1')
ylabel('U_2')
title('$$\hat{C}(F_X,G_Y) | Y = sin(X-10)$$', 'Interpreter', 'Latex');

y = sin(x-5);
subplot(2,2,subplot_idx); subplot_idx = subplot_idx + 1;
[ecop,U,c] = empcopula([x y],K); U1 = U{1}; U2 = U{2};
contour(U1,U2,ecop)
xlabel('U_1')
ylabel('U_2')
title('$$\hat{C}(F_X,G_Y) | Y = sin(X-5)$$', 'Interpreter', 'Latex');

y = sin(x+5);
subplot(2,2,subplot_idx); subplot_idx = subplot_idx + 1;
[ecop,U,c] = empcopula([x y],K); U1 = U{1}; U2 = U{2};
contour(U1,U2,ecop)
xlabel('U_1')
ylabel('U_2')
title('$$\hat{C}(F_X,G_Y) | Y = sin(X+5)$$', 'Interpreter', 'Latex');

y = sin(x+10);
subplot(2,2,subplot_idx); subplot_idx = subplot_idx + 1;
[ecop,U,c] = empcopula([x y],K); U1 = U{1}; U2 = U{2};
contour(U1,U2,ecop)
xlabel('U_1')
ylabel('U_2')
title('$$\hat{C}(F_X,G_Y) | Y = sin(X+10)$$', 'Interpreter', 'Latex');

h = figtitle(figtitlestr);