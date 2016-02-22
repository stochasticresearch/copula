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

function [ C ] = empcopulacdf( U, K, method )
%EMPCOPULACDF Calculates the empirical copula function over the grid
%points with K points in each dimension, using the specified method
% 
% Inputs
%  U - the pseudosamples of the copula (i.e.) F(X) = [F(X_1) ... F(X_D)],
%      should be a M x D input, where M is the number of samples and D is
%      the dimensionality
%  K - the spacing of the gridpoints over which to calculate the empirical
%      copula density
%  method - options are:
%           deheuvels - Use the Deheuvels method of calculating the
%           empirical copula, see: https://en.wikipedia.org/wiki/Copula_(probability_theory)#Empirical_copulas
%
% Outputs
%  C - the copula density

if(strcmpi(method, 'deheuvels'))
    [M,D] = size(U);
    sz = ones(1,D)*K;
    C = zeros(1,K^D);       % we create it as a vector, then reshape at the
                            % end to the proper dimensionality
    u = cell(1,D);
    [u{:}] = ndgrid(linspace(0.01,.99,K));     % we go from .1 to .99 to avoid boundary effects
    
    linearIdx = 1;
    for ii=1:K^D
        u_t = zeros(1,D);
        for jj=1:D
            u_t(jj) = u{jj}(ii);
        end
        mask = bsxfun(@le,U,u_t);
        maskSum = sum(mask,2);
        val = sum(maskSum==D)/M;
        
        C(linearIdx) = val;
        linearIdx = linearIdx + 1;
    end
                            
    C = reshape(C,sz);
else
    error('Only Dehuvels method supported currently!');
end

end

