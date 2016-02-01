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

% A script which contains some test simulations for the Paper "Discrete
% Copula's," published in IEEE Transactions on Fuzzy Systems 2006

function [C] = irreducibleCopula(n)

P = perms(1:n);
C = zeros(size(P,1), n*n);

for zz=1:size(P,1)
    fprintf('Calculating permutation = %s\n', sprintf('%d ', P(zz,:)))
    
    ti = 1;
    for ii=1:n
        for jj=1:n
            C(zz, ti) = func(P(zz,:), ii, jj);
            ti = ti + 1;
        end
    end
end

end

function [y] = func(perm, i, j)

n = length(perm);

y = 0;
for k=1:j
    indic = 1:i;
    if( any(perm(k)==indic) )
        y = y + 1;
    end
end

y = y / n;
end