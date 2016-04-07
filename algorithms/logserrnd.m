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


function [ y ] = logserrnd(n, p)
% Performs the same as LOGSERRND above, but this is the way it was
% implemented in the copula package in R.

y = zeros(n,1);
thresh = 0.95;

if(length(p)==1)
    p = repmat(p, n, 1);
end

for nn=1:n
    if(p(nn)<thresh)
        y(nn) = logserrnd_LS ( p(nn) );
    else
        y(nn) = logserrnd_LK ( p(nn) );
    end
    
end

end

function [ y ] = logserrnd_LS( alpha )
    t = - alpha / log(1 - alpha);
    u = rand();
    p = t;
    x = 1;
    while (u > p) 
        u = u - p;
        x = x + 1;
        p = p * alpha * (x - 1) / x;
    end
    y = x;
end

function [ y ] = logserrnd_LK( alpha )

    h = log(1 - alpha);
    u2 = rand();
    x = 1.0;
    if ( u2 > alpha ) 
        y = fix(x);     % used to be int(x)
    else 
        u1 = rand;
        q = 1 - exp(u1 * h);
        if ( u2 < q * q ) 
            y = fix(1 + log(u2) / log(q));
        elseif ( u2 > q ) 
            y = 1;
        else
            y = 2;
        end
    end

end