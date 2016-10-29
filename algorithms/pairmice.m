function [R] = pairmice( X )
%PAIRRDC - computes pairwise dependency metrics of a given vector
% Inputs:
%  X - a matrix of observations from which pairwise RSDM metrics are
%  computed.
% Outputs:
%  R - a matrix which contains pairwise correlations of all columns in X
%**************************************************************************
%*                                                                        *
%* Copyright (C) 2016  Kiran Karra <kiran.karra@gmail.com>                *
%*                                                                        *
%* This program is free software: you can redistribute it and/or modify   *
%* it under the terms of the GNU General Public License as published by   *
%* the Free Software Foundation, either version 3 of the License, or      *
%* (at your option) any later version.                                    *
%*                                                                        *
%* This program is distributed in the hope that it will be useful,        *
%* but WITHOUT ANY WARRANTY; without even the implied warranty of         *
%* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          *
%* GNU General Public License for more details.                           *
%*                                                                        *
%* You should have received a copy of the GNU General Public License      *
%* along with this program.  If not, see <http://www.gnu.org/licenses/>.  *
%*                                                                        *
%**************************************************************************

n = size(X,2);      % the dimensionality of the dataset
R = zeros(n,n);

for ii=1:n
    x = X(:,ii);
    parfor jj=ii+1:n
        y = X(:,jj);        
        miceVal = mine(x',y');
        R(ii,jj) = miceVal.mic;
    end
end

% assign the lower triangle by copying the upper-triangle of the R matrix
% make matrix symmetric.  The below works because R is initialized to
% zeros.
R=R+R';
R(1:n+1:n*n) = 1;   % set the diagonal to 1

end