function [ p, T ] = mvdistgof( X, sizes, nperm )
%MVDISTGOF Tests if samples from two multivariate distributions of the same
%dimensionality come from the same distribution, based on the paper:
%TESTING FOR EQUAL DISTRIBUTIONS IN HIGH DIMENSION by 
%Gabor J. Szkely and Maria L. Rizzo.  Uses permutation testing to derive a
%p-value.  Ported from the energy package in R.
% 
% Inputs:
%  X - pooled samples from the p distributions to compare
%       dimension [(n_1+n_2+...+n_p) x D]
%  sizes - a vector w/ the info [n_1 n_2 ... n_p]
%  nperm - [optional] the number of permutations to test in the
%          computation of the p-value, defaults to 500
%
% Outputs:
%  p - the computed p-value
%  T - the computed test statistic
%
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

% default the number of permutations to compare to 500 (if arg is not
% provided)
if(nargin<3)
    nperm = 500;
end

n = sum(sizes);
% nsamps = length(sizes);
% D = squareform(pdist(X));
% T = multivariateE(D, nsamps, sizes);
T = multivariateEnergyDistance(X, sizes);

ek = 0;
for ii=1:nperm
    perm = randperm(n);
    X_permute = X(perm,:);
%     D = squareform(pdist(X_permute));
%     T_ii = multivariateE(D, nsamps, sizes);
    T_ii = multivariateEnergyDistance(X_permute, sizes);
    ek = ek + (T < T_ii);
%     fprintf('T_ii=%f T=%f\n', T_ii, T);
end
p = (ek+1) / (nperm+1);

end