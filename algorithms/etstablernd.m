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

function [ St ] = etstablernd( V0, h, alpha )
%RETSTABLE Samples from the exponentially tilted stable distribution
% Sample St ~ \tilde{S}(alpha, 1, (cos(alpha*pi/2)*V_0)^{1/alpha},
%			 V_0*I_{alpha = 1}, h*I_{alpha != 1}; 1)
% with Laplace-Stieltjes transform exp(-V_0((h+t)^alpha-h^alpha)),
% see Nolan's book for the parametrization, via the "fast rejection algorithm",
% see Hofert (2012).
%
% Inputs:
%  V0 - vector of random variates V0
%  h - parameter in [0,infinity)
%  alpha - parameter in (0,1].  WARNING, alpha IS ONLY SUPPORTED TO BE A
%                               SCALAR AND THUS APPLIES TO EVERY V0!!
% Outputs:
%  St - vector of random variates
%
% Acknowledgements - R implementation of retstable by:
%   Marius Hofert, Martin Maechler

V0 = V0(:);
n = length(V0);

% alpha == 1 => St corresponds to a point mass at V0 with Laplace-Stieltjes
% transform exp(-V0*t)
if(alpha == 1.)
    St = V0;
else
    St = zeros(length(V0),1);
    for ii=1:n
        m = max(1, round( V0(ii)* (h^alpha) ) );
        c = (V0(ii)/m)^(1/alpha);
        
        St(ii) = 0;     % will be result after summation
        for kk=1:m
            condition = 1;
            while condition
                St_k = c*stable0rnd(alpha);
                U = rand();
                condition = (U > exp( -h*St_k));
            end
            St(ii) = St(ii) + St_k;
        end
        
    end
end

end