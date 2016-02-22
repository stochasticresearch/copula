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

function [ c ] = empcopulapdf( U, h, K, method )
%EMPCOPULAPDF Calculates the empirical copula density over the grid
%points with K points in each dimension, using the specified method
% 
% Inputs
%  U - the pseudosamples of the copula (i.e.) F(X) = [F(X_1) ... F(X_D)],
%      should be a M x D input, where M is the number of samples and D is
%      the dimensionality
%  K - the spacing of the gridpoints over which to calculate the empirical
%      copula density
%  method - options are:
%           betak - Use Beta Kernels to estimate copula density, as
%                   specified by Copula Density Estimation Book - by Arthur
%                   Charpentier
%
% Outputs
%  c - the copula density

if(strcmpi(method, 'betak'))
    M = size(U,1);
    D = size(U,2);
    sz = ones(1,D)*K;
    c = zeros(1,K^D);       % we create it as a vector, then reshape at the
                            % end to the proper dimensionality
    linearIdx = 1;
    
    for ii=0:K^D-1
        % make the beta pdf param's vector
        gridPoints = zeros(1,D);
        
        value = ii;
        xIdx = D;
        while(value > 0)
            d = mod(value,K);

            gridPoints(xIdx) = d;
            xIdx = xIdx - 1;

            value = floor(value/K);
        end
        gridPoints = gridPoints/(K-1);      % this calculates 0 ... 1 (spacing will be 1/(K+1))
                                            % to maintain K points
        gridPoints = fliplr(gridPoints);    % flip this to calculate column wise so that reshape
                                            % will restore c to expected dimensions
        Kernel_vec = ones(M,1);
        for jj=1:D 
            gridPoint = gridPoints(jj);
            Kernel_vec = Kernel_vec .* betapdf(U(:,jj), gridPoint/h + 1, (1-gridPoint)/h + 1);            
        end

        c(linearIdx) = sum(Kernel_vec)/M;
        linearIdx = linearIdx + 1;
    end
                            
    c = reshape(c,sz);
else
    error('Only beta kernel method supported currently!');
end

end

