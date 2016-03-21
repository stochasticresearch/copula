function [y] = claytonV01rnd(V0, theta0, theta1)
%CLAYTONV01 - LS^{-1}[exp(-V_0psi_0^{-1}(psi_1(t)))]

y = etstablernd(V0, 1, theta0/theta1 );

end