%******************************************************************************
%* 
%* Copyright (C) 2016  Kiran Karra <kiran.karra@gmail.com>
%*
%* This program is free software: you can redistribute it and/or modify
%* it under the terms of the GNU General Public License as published by
%* the Free Software Foundation, either version 3 of the License, or
%* (at your option) any later version.
%*
%* This program is distributed in the hope that it will be useful,
%* but WITHOUT ANY WARRANTY; without even the implied warranty of
%* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%* GNU General Public License for more details.
%*
%* You should have received a copy of the GNU General Public License
%* along with this program.  If not, see <http://www.gnu.org/licenses/>.

function [ y ] = rstable1(n, alpha, beta, gamma, delta, pm)
% Generates a value from the stable distribution.
% Algorithms copied directly from R source code of the copula package
%    - rstable1.R
%    - retstable.c

    y = rstable_c(n, alpha) * gamma + delta;

end


function [ y ] = rstable_c(n, alpha)
    y = cos(pi/2.0*alpha).^(-1.0/alpha) * rstable0(alpha);
end

function [ y ] = rstable0(alpha)
    if(alpha == 1) 
        y = 1.0;
    else
        U = rand();
        while 1
            % generate non-zero exponential random variable
            W = exprnd(1);
            if(W~=0)
                break;
            end
        end
        y = (A(pi*U,alpha)/(W.^(1.0-alpha)) ).^(1.0/alpha);
    end
end

function [ y ] = A(x, alpha)
    Ialpha = 1.0-alpha;
    y = A_3(x, alpha, Ialpha);
end

function [ y ] = A_3(x, alpha, Ialpha)
    y = ( (Ialpha* sinc(Ialpha*x/pi)).^Ialpha ) * ...
        ( ( (alpha * sinc(alpha *x/pi)).^alpha ) ./ sinc(x./pi) );
end