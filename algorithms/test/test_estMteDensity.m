clear;
clc;

% Test script for estMteDensity
x = normrnd(0,1,1000,1);
[f,xi] = ksdensity(x);
mte_estimate = estMteDensity(x);

plot(xi,f,mte_estimate.domain,mte_estimate.density);
grid on;
xlabel('x')
legend('KDE', 'MTE')

refLikelihood = -1*normlike([0,1],x)

mteLLVal = 0;
for ii=1:1000
    mteLLVal = mteLLVal + log(mte_estimate.queryDensity(x(ii)));
end
mteLLVal