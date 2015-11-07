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
neg = 0;

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


N = 500;
K = 250;

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
[ecop, U1, U2] = ecopula([x y]);
contour(U1,U2,ecop)
xlabel('U_1')
ylabel('U_2')
title('$$\hat{C}(F_X,G_Y)$$', 'Interpreter', 'Latex');

[U_gen, Z_sorted, U_emp] = emp_copularnd([x y], N, K);
subplot(4,4,subplot_idx); subplot_idx = subplot_idx + 1;
scatter(U_emp(:,1), U_emp(:,2))
title('Empirical Samples of $$\hat{C}(F_X,G_Y)$$', 'Interpreter', 'Latex');
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
[ecop, U1, U2] = ecopula([x y]);
contour(U1,U2,ecop)
xlabel('U_1')
ylabel('U_2')
title('$$\hat{C}(F_X,G_Y)$$', 'Interpreter', 'Latex');

[U_gen, Z_sorted, U_emp] = emp_copularnd([x y], N, K);
subplot(4,4,subplot_idx); subplot_idx = subplot_idx + 1;
scatter(U_emp(:,1), U_emp(:,2))
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
[ecop, U1, U2] = ecopula([x y]);
contour(U1,U2,ecop)
xlabel('U_1')
ylabel('U_2')
title('$$\hat{C}(F_X,G_Y)$$', 'Interpreter', 'Latex');

[U_gen, Z_sorted, U_emp] = emp_copularnd([x y], N, K);
subplot(4,4,subplot_idx); subplot_idx = subplot_idx + 1;
scatter(U_emp(:,1), U_emp(:,2))
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
[ecop, U1, U2] = ecopula([x y]);
contour(U1,U2,ecop)
xlabel('U_1')
ylabel('U_2')
title('$$\hat{C}(F_X,G_Y)$$', 'Interpreter', 'Latex');

[U_gen, Z_sorted, U_emp] = emp_copularnd([x y], N, K);
subplot(4,4,subplot_idx); subplot_idx = subplot_idx + 1;
scatter(U_emp(:,1), U_emp(:,2))
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

% % polynomial order 4 dependence
% figure;
% subplot_idx = 1;
% 
% y = x.^4;
% titStr = 'Y=X^4';
% if(neg)
%     y = y*-1;
%     titStr = 'Y=-X^4';
% end
% [U_gen, Z_sorted, U_emp] = emp_copularnd([x y], N, K);
% subplot(3,3,subplot_idx); subplot_idx = subplot_idx + 1;
% scatter(U_emp(:,1), U_emp(:,2))
% title(sprintf('%s %s', 'Samples of C(F_X,G_Y) | ', titStr));
% xlabel('U_1')
% ylabel('U_2')
% grid on
% 
% subplot(3,3,subplot_idx); subplot_idx = subplot_idx + 1;
% scatter(x, y, 'r')
% title(titStr);
% xlabel('X')
% ylabel('Y')
% grid on
% 
% subplot(3,3,subplot_idx); subplot_idx = subplot_idx + 1;
% [ecop, U1, U2] = ecopula([x y]);
% contour(U1,U2,ecop)
% xlabel('U_1')
% ylabel('U_2')
% title(sprintf('%s %s','C(F_X,G_Y) | ', titStr));
% 
% % polynomial order 4 dependence
% y = (x.^2-0.25).*(x.^2-1);
% titStr = 'Y=(X^2-0.25)*(X^2-1)';
% if(neg)
%     y = y*-1;
%     titStr = 'Y=-(X^2-0.25)*(X^2-1)';
% end
% [U_gen, Z_sorted, U_emp] = emp_copularnd([x y], N, K);
% subplot(3,3,subplot_idx); subplot_idx = subplot_idx + 1;
% scatter(U_emp(:,1), U_emp(:,2))
% title(sprintf('%s %s', 'Samples of C(F_X,G_Y) | ', titStr));
% xlabel('U_1')
% ylabel('U_2')
% grid on
% 
% subplot(3,3,subplot_idx); subplot_idx = subplot_idx + 1;
% scatter(x, y, 'r')
% title(titStr);
% xlabel('X')
% ylabel('Y')
% grid on
% 
% subplot(3,3,subplot_idx); subplot_idx = subplot_idx + 1;
% [ecop, U1, U2] = ecopula([x y]);
% contour(U1,U2,ecop)
% xlabel('U_1')
% ylabel('U_2')
% title(sprintf('%s %s','C(F_X,G_Y) | ', titStr));
% 
% % sinusoidal dependence
% titStr = 'Y=sin(X)';
% if(neg)
%     y = y*-1;
%     titStr = 'Y= -sin(X)';
% end
% y = sin(x);
% [U_gen, Z_sorted, U_emp] = emp_copularnd([x y], N, K);
% subplot(3,3,subplot_idx); subplot_idx = subplot_idx + 1;
% scatter(U_emp(:,1), U_emp(:,2))
% title(sprintf('%s %s', 'Samples of C(F_X,G_Y) | ', titStr));
% xlabel('U_1')
% ylabel('U_2')
% grid on
% 
% subplot(3,3,subplot_idx); subplot_idx = subplot_idx + 1;
% scatter(x, y, 'r')
% title(titStr);
% xlabel('X')
% ylabel('Y')
% grid on
% 
% subplot(3,3,subplot_idx); subplot_idx = subplot_idx + 1;
% [ecop, U1, U2] = ecopula([x y]);
% contour(U1,U2,ecop)
% xlabel('U_1')
% ylabel('U_2')
% title(sprintf('%s %s','C(F_X,G_Y) | ', titStr));
% 
% h = figtitle(figtitlestr);
% 
% % exponential dependence
% figure;
% subplot_idx = 1;
% 
% y = exp(x);
% titStr = 'Y=e^X';
% if(neg)
%     y = y*-1;
%     titStr = 'Y=-e^X';
% end
% [U_gen, Z_sorted, U_emp] = emp_copularnd([x y], N, K);
% subplot(3,3,subplot_idx); subplot_idx = subplot_idx + 1;
% scatter(U_emp(:,1), U_emp(:,2))
% title(sprintf('%s %s', 'Samples of C(F_X,G_Y) | ', titStr));
% xlabel('U_1')
% ylabel('U_2')
% grid on
% 
% subplot(3,3,subplot_idx); subplot_idx = subplot_idx + 1;
% scatter(x, y, 'r')
% title(titStr);
% xlabel('X')
% ylabel('Y')
% grid on
% 
% subplot(3,3,subplot_idx); subplot_idx = subplot_idx + 1;
% [ecop, U1, U2] = ecopula([x y]);
% contour(U1,U2,ecop)
% xlabel('U_1')
% ylabel('U_2')
% title(sprintf('%s %s','C(F_X,G_Y) | ', titStr));
% 
% % circular dependence
% y = sqrt(x.^2+1);
% titStr = 'Y=sqrt(X^2+1)';
% if(neg)
%     y = y*-1;
%     titStr = 'Y=-sqrt(X^2+1)';
% end
% [U_gen, Z_sorted, U_emp] = emp_copularnd([x y], N, K);
% subplot(3,3,subplot_idx); subplot_idx = subplot_idx + 1;
% scatter(U_emp(:,1), U_emp(:,2))
% title(sprintf('%s %s', 'Samples of C(F_X,G_Y) | ', titStr));
% xlabel('U_1')
% ylabel('U_2')
% grid on
% 
% subplot(3,3,subplot_idx); subplot_idx = subplot_idx + 1;
% scatter(x, y, 'r')
% title(titStr);
% xlabel('X')
% ylabel('Y')
% grid on
% 
% subplot(3,3,subplot_idx); subplot_idx = subplot_idx + 1;
% [ecop, U1, U2] = ecopula([x y]);
% contour(U1,U2,ecop)
% xlabel('U_1')
% ylabel('U_2')
% title(sprintf('%s %s','C(F_X,G_Y) | ', titStr));
% 
% % circular dependence
% y = real(log(x));
% titStr = 'Y=Re[log(X)]';
% if(neg)
%     y = y*-1;
%     titStr = 'Y=-Re[log(X)]';
% end
% [U_gen, Z_sorted, U_emp] = emp_copularnd([x y], N, K);
% subplot(3,3,subplot_idx); subplot_idx = subplot_idx + 1;
% scatter(U_emp(:,1), U_emp(:,2))
% title(sprintf('%s %s', 'Samples of C(F_X,G_Y) | ', titStr));
% xlabel('U_1')
% ylabel('U_2')
% grid on
% 
% subplot(3,3,subplot_idx); subplot_idx = subplot_idx + 1;
% scatter(x, y, 'r')
% title(titStr);
% xlabel('X')
% ylabel('Y')
% grid on
% 
% subplot(3,3,subplot_idx); subplot_idx = subplot_idx + 1;
% [ecop, U1, U2] = ecopula([x y]);
% contour(U1,U2,ecop)
% xlabel('U_1')
% ylabel('U_2')
% title(sprintf('%s %s','C(F_X,G_Y) | ', titStr));
% 
% h = figtitle(figtitlestr);

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
[ecop, U1, U2] = ecopula([x y]);
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
[ecop, U1, U2] = ecopula([x y]);
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
[ecop, U1, U2] = ecopula([x y]);
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

