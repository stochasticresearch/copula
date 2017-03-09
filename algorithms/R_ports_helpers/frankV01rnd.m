function [y] = frankV01rnd(V0, theta0, theta1, rej, approx)
%FRANKV01 - LS^{-1}[exp(-V_0psi_0^{-1}(psi_1(t)))]

if(nargin<4)
    rej = 1;
end
if(nargin<5)
    approx = 10000;
end

y = F01Frankrnd(V0, theta0, theta1, rej, approx);

end