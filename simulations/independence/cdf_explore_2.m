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
domainX = linspace(min(X), max(X), numEcdfPts);
FX = ksdensity(X, domainX, 'function', 'cdf');
empInfoX = rvEmpiricalInfo(domainX, [], FX, 0);

domainY = linspace(min(Y), max(Y), numEcdfPts);
FY = ksdensity(Y, domainY, 'function', 'cdf');
% [FY, domainY] = ecdf(domainY); FY = FY'; domainY = domainY';
empInfoY = rvEmpiricalInfo(domainY, [], FY, 0);

domainX1 = linspace(min(X1),max(X1), numEcdfPts);
empInfoX1 = makedist('Uniform', 'lower', min(X1), 'upper', max(X1));
FX1 = empInfoX1.cdf(domainX1);

domainX2 = linspace(min(X2),max(X2), numEcdfPts);
empInfoX2 = makedist('Uniform', 'lower', min(X2), 'upper', max(X2));
FX2 = empInfoX2.cdf(domainX2);

y_vec = min(Y):0.1:max(Y)+0.1;
FY_calculated = zeros(1,length(y_vec));
for y_idx=1:length(y_vec)
    y = y_vec(y_idx);
    FY_calculated(y_idx) = 1/3*(1-empInfoX1.cdf(fun_g1_inv(y))) + ...
                           1/3*(empInfoX2.cdf(fun_g2_inv(y))) + ...
                           1/3*(empInfoX3(fun_g3_inv(y)));
end

x = linspace(min(X),max(X),100);

for x_idx=1:length(x)
    xx = x(x_idx);
    FX_x(x_idx) = empInfoX.cdf(xx);
    FY_y(x_idx) = 1/3*(1-empInfoX1.cdf(fun_g1_inv(fun_g1(xx)))) + ...
           1/3*(empInfoX2.cdf(fun_g2_inv(fun_g2(xx)))) + ...
           1/3*(empInfoX3(fun_g3_inv(fun_g3(xx)))) - 1/3;

end


subplot(3,2,1); scatter(X,Y); xlabel('X'); ylabel('Y'); grid on;

% instead of using pobs directly, we do the sort for Y b/c it has duplicate
% values ...
pobsX = pobs(X);
data_sorted = sort(Y);
[~, pobsY] = ismember(Y,data_sorted);
pobsY = pobsY/length(Y);

subplot(3,2,2); scatter(pobsX, pobsY); xlabel('u'); ylabel('v'); grid on;
subplot(3,2,3); plot(domainX, FX); xlabel('x'); ylabel('F_X(x)'); grid on;

subplot(3,2,4); 
plot(domainY, FY); xlabel('y'); ylabel('F_Y(y)'); grid on;
hold on;
plot(y_vec, FY_calculated, 'o');

subplot(3,2,[5 6])
plot(FX_x,FY_y); grid on; xlabel('F_X(x)'); ylabel('F_Y(y)');

end

function [g3_inv_val] = fun_g3_inv(y)

if(y>=4)
    g3_inv_val = 1;
else
    g3_inv_val = 0;
end

end

function [FX3_val] = empInfoX3(y)
FX3_val = (y>=1);
end

function [g1_val] = fun_g1(x)

g1_val = 4*(x).^2;

end

function [g2_val] = fun_g2(x)

g2_val = 4*(x).^2;

end

function [g3_val] = fun_g3(x)

g3_val = 4;

end

function [g1_inv_val] = fun_g1_inv(y)

g1_inv_val = -sqrt(y./4);

end

function [g2_inv_val] = fun_g2_inv(y)

g2_inv_val = sqrt(y./4);

end