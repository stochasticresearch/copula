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

function [ y ] = stable0rnd(alpha)
    if(alpha == 1) 
        y = 1.0;
    else
        U = rand();
        condition = 1;
        while condition
            % generate non-zero exponential random variable
            W = exprnd(1);
            condition = (W==0);
        end
        y = ( A(pi*U,alpha)/(W^(1.0-alpha)) )^(1.0/alpha);
    end
end

function [ y ] = A(x, alpha)
    Ialpha = 1.0-alpha;
    y = A_3(x, alpha, Ialpha);
end

function [ y ] = A_3(x, alpha, Ialpha)
    y = sin(alpha*x)^alpha * sin((Ialpha)*x)^(Ialpha) / sin(x);
end