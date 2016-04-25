function [ e ] = multivariateEnergyDistance( X, sizes )
%MULTIVARIATEENERGYDISTANCE Computes the energy between p d-dimensional
%probability distributions, given samples from these distributions, based 
%on the paper: TESTING FOR EQUAL DISTRIBUTIONS IN HIGH DIMENSION by 
%Gabor J. Szkely and Maria L. Rizzo.  Uses permutation testing to derive a
%p-value.  Ported from the energy package in R.
% 
% Inputs:
%  X - pooled samples from the p distributions to compare
%       dimension [(n_1+n_2+...+n_p) x D]
%  sizes - a vector w/ the info [n_1 n_2 ... n_p], where n_i is the number
%  of samples from the i^th d-dimensional distribution
%
% Outputs:
%  e - the computed energy
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

nsamps = length(sizes);
D = squareform(pdist(X));
e = multivariateE(D, nsamps, sizes);

end

function [e] = multivariateE( D, nsamps, sizes)

startIdxs = ones(1,nsamps);
for ii=2:nsamps
    startIdxs(ii) = startIdxs(ii-1) + sizes(ii-1);
end

e = 0;
for ii=1:nsamps
    m = sizes(ii);
    for jj=ii+1:nsamps
        n = sizes(jj);
        e = e + bivariateE(D, m, n);
    end
end

end

function [e] = bivariateE( D, m, n )
% computes the test statistic defined by Equation 5 in the paper referenced
% above

xx = D(1:m,1:m);
sumxx = sum(xx(:))/2;
sumxx = sumxx* 2/(m*m);

yy = D(m+1:end,m+1:end);
sumyy = sum(yy(:))/2;
sumyy = sumyy * 2/(n*n);

xy = D(1:m,m+1:end);
sumxy = sum(xy(:));
sumxy = sumxy/(m*n);

e = m*n/(m+n) * (2*sumxy - sumxx - sumyy);

end