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

function [ H ] = empdistrnd( U, X )
%EMPDISTRND Generates samples from a specified joint distribution
% Inputs:
%  U - a [N x D] matrix of dependent uniform random variables (with the
%                desired dependency structure)
%  X - samples of the desired marginal distribution for output joint
%      distribution.  Must be a matrix of dimension [alpha*N x D], where
%      alpha is >= 1.  Each column d (of D), represents the samples of the
%      desired marginal distribution for that dimension d.  
% Outputs:
%  H - samples of the desired joint distribution, with the desired
%      dependency structure given by U and the marginal distributions,
%      given by X

N = size(X,1);
D = size(X,2);

n_gen = size(U,1);

Z_sorted = zeros(N,D);
for d=1:D
    [Z_sorted(:,d), ~] = sort(X(:,d));
end
if(N/n_gen > 1 )
    Z_sorted = downsample(Z_sorted,floor(N/n_gen));
end

% apply inverse transform to generate samples from this joint distribution
H = zeros(n_gen,D);
for d=1:D
    for j=1:n_gen
        un = U(j,d)*n_gen;
        i = ceil(un);
        H(j,d) = Z_sorted(i,d);
    end
end

end

