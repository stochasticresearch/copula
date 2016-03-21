function [y] = gumbelV01rnd(V0, theta0, theta1)
%GUMBELV01 - LS^{-1}[exp(-V_0psi_0^{-1}(psi_1(t)))]

% Sample from S(1,1,0,V0;1)
% with Laplace-Stieltjes transform exp(-V0*t)
  
alpha = theta0/theta1;
if(alpha==1)
    y = V0;
else
    y = stable1rnd(length(V0), alpha, (cos(pi/2.0*alpha)*V0).^(1/alpha), 0);
end

end