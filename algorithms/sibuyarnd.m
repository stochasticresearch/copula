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

function [ y ] = sibuyarnd( n, alpha, gamma_1_a )

y = zeros(n,1);
for ii=1:n
    y(ii) = sibuyarnd_single(alpha, gamma_1_a);
end

end

function [ y ] = sibuyarnd_single( alpha, gamma_1_a )
%SIBUYARND generates a random variate from the Sibuya (alpha) distribution
%with CDF F(n) = 1-1/(n*B(n,1-alpha)), n in IN, with Laplace-Stieltjes 
% transform 1-(1-exp(-t))^alpha via the algorithm of Hofert (2011).
%
% Inputs:
%  alpha - theta0/theta1 in (0,1]
% Outputs:
%  y - a random variate from F
%
% Acknowledgements - R implementation of rSibuya by:
%   Marius Hofert, Martin Maechler

DBL_EPSILON = 2.2204460492503131e-16;

U = rand();
if(U <= alpha)
	y = 1.;    
else  % < alpha < U < 1 */
	xMax = 1.0/DBL_EPSILON; % ==> floor(x) == ceil(x)  for x >= xMax
	Ginv = ((1-U)*gamma_1_a)^(-1/alpha);
    fGinv = floor(Ginv);
    
    if(Ginv > xMax)
        y = fGinv;
    elseif( (1-U) < (1/(fGinv*beta(fGinv, 1.-alpha))) )
        y = ceil(Ginv);
    else
        y = fGinv;
    end
    
end

end