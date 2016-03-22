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

function [ y ] = logrnd( n, p, Ip )
%RLOG Sample a Log(p) distribution with the algorithm "LK" of Kemp (1981).
%Generating random variates from a Log(p) distribution with PMF
%          p_k = p^k/(-log(1-p)k), k in IN,

% Inputs:
%  n - the # of variables to generate
%  p - value between (0,1)
%
% Outputs:
%  y - random variate from this distribution
%
% Acknowledgements - R implementation of rLog by:
%   Marius Hofert, Martin Maechler

if(p <= 0. ||  p > 1.)
	error('rLog(): p must be inside (0,1)');
else
    y = zeros(n,1);
    for ii=1:n
        y(ii) = logrnd_single(p, Ip);
    end
end

end

function [ y ] = logrnd_single( p, Ip )

U = rand();
if(U>p)
    y = 1;
else
    if(p<=0.5)
        Q = -expm1(log1p(-p) * rand()); % = 1-(1-p)^unif
        logQ = log(Q);
    else
        iQ = Ip^rand(); % = (1-p)^unif
        Q = 1 - iQ;
        logQ = log1p(-iQ);
    end
    if(U<(Q*Q))
        y = floor(1. + log(U)/logQ);
    else
        if(U>Q)
            y = 1;
        else
            y = 2;
        end
    end
end

end