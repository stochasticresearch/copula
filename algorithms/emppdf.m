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
%* 
%**************************************************************************

function [ f, xi ] = emppdf( x, isdiscrete )

M = size(x,1);

if(isdiscrete)
    [~,xi] = empcdf( x, isdiscrete );
    f = zeros(1,length(xi));
    idx = 1;
    for jj=1:length(xi)
        f(idx) = sum(x==xi(jj))/M;
        idx = idx + 1;
    end
else
    [f,xi] = ksdensity(x);
end

f = f(:);
xi = xi(:);
f = f';
xi = xi';

end