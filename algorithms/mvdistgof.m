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

function [ h, p ] = mvdistgof( X, Y, alpha, nperm )
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
if(nargin<4)
    nperm = 10000;
end

if(size(X,2)~=size(Y,2))
    error('X and Y must be of the same dimension!');
end

T = energyteststatistic(X, Y);

n1 = size(X,1);
n2 = size(Y,1);
n = n1 + n2;
XY_concat = [X;Y];

p = 0;
for ii=1:nperm
    permidx = randperm(n)';
    XY_perm = XY_concat(permidx,:);
    X_permute = XY_perm(1:n1,:);
    Y_permute = XY_perm(n1+1,:);
    T_ii = energyteststatistic(X_permute, Y_permute);
    p = p + (T_ii < T);
    fprintf('T_ii=%f T=%f\n', T_ii, T);
end
p = p / nperm;
h = ~(p < alpha);

end


function [e] = energyteststatistic( X, Y )
% computes the test statistic defined by Equation 5 in the paper referenced
% above

n1 = size(X,1);
n2 = size(Y,1);

v1 = 0;
for ii=1:n1
    for mm=1:n2
        v1 = v1 + norm(X(ii,:) - Y(mm,:));
    end
end

v2 = 0;
for ii=1:n1
    for jj=1:n1
        v2 = v2 + norm(X(ii,:) - X(jj,:));
    end
end

v3 = 0;
for ll=1:n2
    for mm=1:n2
        v3 = v3 + norm(Y(ll,:) - Y(mm,:));
    end
end

e = n1*n2/(n1 + n2) * ( (2/(n1*n2))*v1 + (1/n1^2)*v2 + (1/n2^2)*v3 );

end