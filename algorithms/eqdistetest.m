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

function [ p ] = eqdistetest( X, Y, nperm )
%MVDISTGOF Test if samples from two multivariate distributions of the same
%dimensionality come from the same distribution, based on the paper:
%TESTING FOR EQUAL DISTRIBUTIONS IN HIGH DIMENSION by 
%Gabor J. Szkely and Maria L. Rizzo.  Uses permutation testing to derive a
%p-value
% 
% Inputs:
%  X - samples from the first distribution, dimension [n1 x D]
%  Y - samples from the second distribution, dimension - [n2 x D]
%  alpha - significance level for p-value test
%  nperm - [optional] the number of permutations to test in the
%          computation of the p-value, defaults to 10000
%
% Outputs:
%  h - 1 if X and Y are from the same distribution according to this test
%  p - the computed p-value

% default the number of permutations to compare to 10000 (if arg is not
% provided)
if(nargin<3)
    nperm = 10000;
end

if(size(X,2)~=size(Y,2))
    error('X and Y must be of the same dimension!');
end

n1 = size(X,1);
n2 = size(Y,1);
n = n1 + n2;
XY_concat = [X;Y];
D = squareform(pdist(XY_concat));

perm = 1:n;
T = multisampleE(D, 2, [n1 n2], perm)
% T = twosampleE(D, n1, n2, perm)

ek = 0;
h = 0;
for ii=1:nperm
    perm = randperm(n);
    T_ii = multisampleE(D, 2, [n1 n2], perm);
    ek = ek + (T < T_ii);
    fprintf('T_ii=%f T=%f\n', T_ii, T);
end
p = (ek+1) / (nperm+1);

end

function [e] = multisampleE( D, nsamps, sizes, perm)

startIdxs = ones(1,nsamps);
for ii=2:nsamps
    startIdxs(ii) = startIdxs(ii-1) + sizes(ii-1);
end

e = 0;
for ii=1:nsamps
    m = sizes(ii);
    for jj=ii+1:nsamps
        n = sizes(jj);
        e = e + twosampleE(D, m, n, perm);
    end
end

end

function [e] = twosampleE( D, m, n, perm )
% computes the test statistic defined by Equation 5 in the paper referenced
% above

sumxx = 0;
for ii=1:m
    for jj=ii+1:m
        sumxx = sumxx + D(perm(ii),perm(jj));
    end
end
sumxx = sumxx* 2/(m*m);

sumyy = 0;
for ii=1:n
    for jj=ii+1:n
        sumyy = sumyy + D(perm(ii+m),perm(jj+m));
    end
end
sumyy = sumyy * 2/(n*n);

sumxy = 0;
for ii=1:m
    for jj=1:n
        sumxy = sumxy + D(perm(ii),perm(jj+m));
    end
end
sumxy = sumxy/(m*n);

e = m*n/(m+n) * (2*sumxy - sumxx - sumyy);

end