function [] = cdf_explore()
% shows how the cdf is affected by the function g(x)

M = 5000;
X = rand(M,1)*2-1;
Y = 4*(X).^2;

X1_idx = find(X<=0);
X2_idx = find(X>0);

X1 = X(X1_idx);
X2 = X(X2_idx);

% g1 = 4*(X1).^2;
% g1_inv = -sqrt(g1./4);
% g2 = 4*(X2).^2;
% g2_inv = sqrt(g2./4);

g1 = fun_g1(X1);
g2 = fun_g2(X2);
g1_inv = fun_g1_inv(g1);
g2_inv = fun_g2_inv(g2);

numEcdfPts = 250;

% here, we force X, X1, and X2 to be all uniform, rather than using
% rvEmpiricalInfo object, b/c we know that they are actually uniform ...,
% to get rid of the approximation (and thus curvature in the graph)
% domainX = linspace(min(X),max(X),numEcdfPts);
% FX = ksdensity(X, domainX, 'function', 'cdf');
% empInfoX = rvEmpiricalInfo(domainX, [], FX, 0);
domainX = linspace(min(X), max(X), numEcdfPts);
empInfoX = makedist('Uniform', 'lower', min(X), 'upper', max(X));
FX = empInfoX.cdf(domainX);

domainY = linspace(min(Y), max(Y), numEcdfPts);
FY = ksdensity(Y, domainY, 'function', 'cdf');
empInfoY = rvEmpiricalInfo(domainY, [], FY, 0);

% domainX1 = linspace(min(X1),max(X1), numEcdfPts);
% FX1 = ksdensity(X1, domainX1, 'function', 'cdf');
% empInfoX1 = rvEmpiricalInfo(domainX1, [], FX1, 0);
domainX1 = linspace(min(X1),max(X1), numEcdfPts);
empInfoX1 = makedist('Uniform', 'lower', min(X1), 'upper', max(X1));
FX1 = empInfoX1.cdf(domainX1);

% domainX2 = linspace(min(X2),max(X2), numEcdfPts);
% FX2 = ksdensity(X2, domainX2, 'function', 'cdf');
% empInfoX2 = rvEmpiricalInfo(domainX2, [], FX2, 0);
domainX2 = linspace(min(X2),max(X2), numEcdfPts);
empInfoX2 = makedist('Uniform', 'lower', min(X2), 'upper', max(X2));
FX2 = empInfoX2.cdf(domainX2);

u = zeros(1,length(X));
v = zeros(1,length(Y));
u1 = zeros(1,length(X1));
u2 = zeros(1,length(X2));
v1 = zeros(1,length(u1));
v2 = zeros(1,length(u2));

for ii=1:length(u)
    u(ii) = empInfoX.cdf(X(ii));
end
for ii=1:length(v)
    v(ii) = empInfoY.cdf(Y(ii));
end

for ii=1:length(u1)
    u1(ii) = empInfoX.cdf(X1(ii));
end
for ii=1:length(u2)
    u2(ii) = empInfoX.cdf(X2(ii));
end
for ii=1:length(v1)
    v1(ii) = empInfoX.cdf(g1_inv(ii));
end
for ii=1:length(v2)
    v2(ii) = empInfoX.cdf(g2_inv(ii));
end

y_vec = [min(Y):0.25:max(Y)];
uv = zeros(1,length(y_vec));
for y_idx=1:length(y_vec)
    y = y_vec(y_idx);
    uv(y_idx) = (1-empInfoX1.cdf(fun_g1_inv(y)))*0.25 + (empInfoX2.cdf(fun_g2_inv(y)))*0.75;
end


figure;

subplot(2,5,1);
scatter(X1, g1); hold on;
scatter(X2, g2);
legend('g_1', 'g_2');
grid on;

subplot(2,5,2);
scatter(g1, g1_inv); hold on;
scatter(g2, g2_inv);
legend('g_1^{-1}', 'g_2^{-1}');
grid on;

subplot(2,5,3);
scatter(u,v);
grid on;
xlabel('u'); ylabel('v');

subplot(2,5,4);
scatter(pobs(X), pobs(Y)); grid on;
grid on;
xlabel('pobs(x)'); ylabel('pobs(y)');

subplot(2,5,6);
scatter(u1,v1);
grid on;
xlabel('u'); ylabel('v_1');

subplot(2,5,7);
scatter(u2,v2);
grid on;
xlabel('u'); ylabel('v_2');

subplot(2,5,8);
scatter([u1 u2],[empInfoY.cdf(4*(1)^2)-2*v1 2*v2-1]);
grid on;
xlabel('u'); ylabel('1-2v_1 + 2v_2-1');

% subplot(2,5,9);
% scatter(y_vec,uv/max(uv));
% grid on;



subplot(2,5,5);
plot(domainX, FX); grid on;
xlabel('x'); ylabel('F_X');

subplot(2,5,[9 10]);
plot(domainY, FY); grid on;
xlabel('y'); ylabel('F_Y');
hold on;
scatter(y_vec,uv, 'r');

figure;
subplot(1,3,1);
plot(domainX1, FX1); grid on;
xlabel('x_1'); ylabel('F_{X_1}');

subplot(1,3,2);
plot(domainX2, FX2); grid on;
xlabel('x_2'); ylabel('F_{X_2}');

subplot(1,3,3);
x = linspace(min(X), max(X), 100);
FX_x = empInfoX.cdf(x);

FY_y = 0.25-0.25*(empInfoX1.cdf(fun_g1_inv(fun_g1(x)))) + ...
            0.75*(empInfoX2.cdf(fun_g2_inv(fun_g2(x))));


plot(FX_x,FY_y);

end

function [g1_val] = fun_g1(x)

g1_val = 4*(x).^2;

end

function [g2_val] = fun_g2(x)

g2_val = 4*(x).^2;

end

function [g1_inv_val] = fun_g1_inv(y)

g1_inv_val = -sqrt(y./4);

end

function [g2_inv_val] = fun_g2_inv(y)

g2_inv_val = sqrt(y./4);

end