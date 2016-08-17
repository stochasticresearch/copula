% shows how the cdf is affected by the function g(x)
clear;
clc;

M = 2000;
x = rand(M,1)*2-1;
y = 4*(x).^2;

x_less0_idx = find(x<=0);
x_greater0_idx = find(x>0);

x_less0 = x(x_less0_idx);
x_greater0 = x(x_greater0_idx);

g1 = 4*(x_less0).^2;
g1_inv = -sqrt(g1./4);
g2 = 4*(x_greater0).^2;
g2_inv = sqrt(g2./4);

domainX = linspace(min(x),max(x),100);  % 100 ecdf points
FX = ksdensity(x, domainX, 'function', 'cdf');
domainY = linspace(min(y), max(y), 100);
FY = ksdensity(y, domainY, 'function', 'cdf');

empInfoX = rvEmpiricalInfo(domainX, [], FX, 0);
empInfoY = rvEmpiricalInfo(domainY, [], FY, 0);

u = pobs(x);
v = zeros(1,length(y));
u1 = zeros(1,length(x_less0));
u2 = zeros(1,length(x_greater0));
v1 = zeros(1,length(u1));
v2 = zeros(1,length(u2));
v3 = zeros(1,length(y));

for ii=1:length(u1)
    u1(ii) = empInfoX.cdf(x_less0(ii));
end
for ii=1:length(u2)
    u2(ii) = empInfoX.cdf(x_greater0(ii));
end
for ii=1:length(v1)
    v1(ii) = empInfoX.cdf(g1_inv(ii));
end
for ii=1:length(v2)
    v2(ii) = empInfoX.cdf(g2_inv(ii));
end
for ii=1:length(v)
    v(ii) = empInfoY.cdf(y(ii));
end

% v1 = pobs(g1_inv); v1 = v1';
% v2 = pobs(g2_inv); v2 = v2';

subplot(2,5,1);
scatter(x_less0, g1); hold on;
scatter(x_greater0, g2);
legend('g_1', 'g_2');
grid on;

subplot(2,5,2);
scatter(x_less0, g1_inv); hold on;
scatter(x_greater0, g2_inv);
legend('g_1^{-1}', 'g_2^{-1}');
grid on;

subplot(2,5,3);
scatter(u,v);
grid on;
xlabel('u'); ylabel('v');

subplot(2,5,4);
scatter(pobs(x), pobs(y)); grid on;
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

subplot(2,5,5);
plot(domainX, FX); grid on;
grid on;
xlabel('x'); ylabel('F_X');

subplot(2,5,10);
plot(domainY, FY); grid on;
grid on;
xlabel('y'); ylabel('F_Y');
