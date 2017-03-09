function [y] = claytonV01rnd(V0, theta0, theta1)
%CLAYTONV01 - LS^{-1}[exp(-V_0psi_0^{-1}(psi_1(t)))]

h = 1;
alpha = theta0/theta1;
y = etstablernd(V0, h, alpha );

end