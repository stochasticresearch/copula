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

function [ y ] = F01Joernd( V0, alpha, gamma_1_a, approx )
%F01JOERND Sample V01 ~ F01 with Laplace-Stieltjes transform 
%((1-(1-exp(-t))^alpha))^V0
% Inputs:
%  V0 - a vector of V0 parameters
%  alpha - theta0/theta1 in (0,1]
%  approx - largest number of summands before asymptotics is used  
%
% Outputs:
%  y - a random variate from F
%
% Acknowledgements - R implementation of rF01Joe by:
%   Marius Hofert, Martin Maechler

y = zeros(size(V0));
for idx=1:length(V0)
    y(idx) = F01Joernd_single(V0(idx), alpha, gamma_1_a, approx);
end

end

function [ y ] = F01Joernd_single( V0, alpha, gamma_1_a, approx )

if(V0 > approx)
    y = V0^(1.0/alpha)*rstable0(alpha);
else
    y = sum(sibuyarnd(fix(V0),alpha, gamma_1_a));
end

end