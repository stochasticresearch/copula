%**************************************************************************
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
%**************************************************************************

function [ y ] = stable1rnd(n, alpha, gamma, delta)
% Generates n random variates from the stable distribution.
% Algorithms copied directly from R source code of the copula package
%    - rstable1.R
%    - retstable.c
% Only valid when beta = 1, pm = 1

    y = stable_c_rnd(n, alpha) .* gamma + delta;

end


function [ y ] = stable_c_rnd(n, alpha)
    y = zeros(n,1);
    for ii=1:n
        y(ii) = cos(pi/2.0*alpha).^(-1.0/alpha) * stable0rnd(alpha);
    end
end

