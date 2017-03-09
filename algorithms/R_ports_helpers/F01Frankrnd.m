%**************************************************************************
%*                                                                        *
%* Copyright (C) 2016  Kiran Karra <kiran.karra@gmail.com>                *
%*                                                                        *
%* This program is free software: you can redistribute it and/or modify   *
%* it under the terms of the GNU General Public License as published by   *
%* the Free Software Foundation, either version 3 of the License, or      *
%* (at your option) any later version.                                    *
%*                                                                        *
%* This program is distributed in the hope that it will be useful,        *
%* but WITHOUT ANY WARRANTY; without even the implied warranty of         *
%* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          *
%* GNU General Public License for more details.                           *
%*                                                                        *
%* You should have received a copy of the GNU General Public License      *
%* along with this program.  If not, see <http://www.gnu.org/licenses/>.  *
%*                                                                        *
%**************************************************************************

function [ V ] = F01Frankrnd( V0, theta0, theta1, rej, approx )
%F01JOERND Generate a vector of variates V01 ~ F01 with Laplace-Stieltjes 
% transform ((1-(1-exp(-t)*(1-e^(-theta1)))^alpha)/(1-e^(-theta0)))^V0.
% Inputs:
%  V0 - vector of random variates from F0
%  theta0 - parameter in (0,infinity)
%  theta1 - parameter in [theta0, infinity)
%  rej - if V0*theta_0*p0^(V0-1) > rej a rejection
%        from F01 of Joe is applied (otherwise, the sum is
%        sampled via a logarithmic envelope for the summands)
%  approx - approx largest number of summands before asymptotics is used
%
% Outputs:
%  y - vector of random variates V01
%
% Acknowledgements - R implementation of rF01Frank by:
%   Marius Hofert, Martin Maechler
p0 = -expm1(-theta0);
p1 = -expm1(-theta1);
iAlpha = (theta1 - theta0) / theta1;
gamma_1_a = gamma(iAlpha);

V = zeros(size(V0));
for idx=1:length(V0)
    V(idx) = F01Frankrnd_single(V0(idx), theta0, theta1, p0, p1, gamma_1_a, rej, approx);
end

end

function [ V ] = F01Frankrnd_single( V0, theta0, theta1, p0, p1, gamma_1_a, rej, approx )

alpha = theta0 / theta1;
iAlpha = (theta1-theta0)/theta1;

if(V0*theta0*p0^(V0-1.0) > rej)
    condition = 1;
    while(condition)
        U = rand();
        V = F01Joernd(V0, alpha,gamma_1_a, approx);
        condition = U > (p1^V);
    end
else
    Ip = exp(-theta1);
    V = 0;
    for j=1:fix(V0)
        condition = 1;
        while(condition)
            U = rand();
            X = logrnd(1,p1,Ip);
            condition = ( U*(X-alpha) ) > ( 1.0/beta(X, iAlpha) );
        end
        V = V + X;
    end
end

end