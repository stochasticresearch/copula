function [] = cdf_explore_2()

M = 5000;
X = rand(M,1)*3-1;
Y = (4*(X).^2).*(X<=1) + 4*(X>1);

X1_idx = find(X<=0);
X2_idx = find(X>0 & X<=1);
X3_idx = find(X>1);

X1 = X(X1_idx);
X2 = X(X2_idx);
X3 = X(X3_idx);

numEcdfPts = 250;
domainY = linspace(min(Y), max(Y), numEcdfPts);
FY = ksdensity(Y, domainY, 'function', 'cdf');
empInfoY = rvEmpiricalInfo(domainY, [], FY, 0);

subplot(1,3,1); scatter(X,Y); xlabel('X'); ylabel('Y'); grid on;
subplot(1,3,2); scatter(pobs(X), pobs(Y)); xlabel('u'); ylabel('v'); grid on;
subplot(1,3,3); plot(domainY, FY); xlabel('y'); ylabel('F_Y(y)'); grid on;

end

function [g1_val] = fun_g1(x)

g1_val = 4*(x).^2;

end

function [g2_val] = fun_g2(x)

g2_val = 4*(x).^2;

end

function [g3_va] = fun_g3(x)

g3_val = 4;

end

function [g1_inv_val] = fun_g1_inv(y)

g1_inv_val = -sqrt(y./4);

end

function [g2_inv_val] = fun_g2_inv(y)

g2_inv_val = sqrt(y./4);

end