function [ h, mapping, prob, combos ] = calcpmf( X )
%CALCPMF - computes the empirical discrete probability
%of occurance of each of the possible outcomes in a multivariate discrete
%probability distribution
% Inputs:
%  X - input M x N matrix, with N being the dimensionality of the discrete
%      distribution, and M being the number of samples.  It is assumed that
%      the data in each colum (i.e. each dimension) is a subset of the
%      natural numbers (i.e. 1, 2, ...) starting with 1, and ending with
%      some value.  The maximum value of each column will be used to
%      determine the the domain of the random variable .  The domain for
%      each marginal distribution will be 1:max(X(:,i)) for the ith
%      marginal distribution, represented by the i^th column of data in X
% Outputs
%  h - the contingency table generated from the data.  The dimensions of
%      this matrix are determined by the number of unique values in each
%      dimension
%  mapping - a cell array which contains the mapping between each element
%            in the contingency table and the associated values of the
%            discrete random variable for each dimension
%  prob - a vector containing the probabilities for each of the
%         combinations, referenced in the combo vector
%  combo - the combo for which the probability was calculated
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
%**************************************************************************

[M,N] = size(X);

if(N==1)
    numUniqueVals = length(unique(X));
    combos = 1:numUniqueVals;
    hSize = [numUniqueVals 1];
else
    % store the number of unique values for each dimension of the discrete
    % distribution
    hSize = zeros(1,N);
    for nn=1:N
        hSize(nn) = max(X(:,nn));
    end

    combos = combvec(1:hSize(1),1:hSize(2));
    for nn=3:N
        combos = combvec(combos, 1:hSize(nn));
    end
end
h = zeros(hSize);
mapping = cell(hSize);

% for each combo, count the number of times they occur in the input data
% and compute the associated probability
combos = combos';
prob = zeros(size(combos,1),1);

for ii=1:size(combos,1)
    combo = combos(ii,:);

    % count the number of times this combo appears in the data set
    cnt = sum(ismember(X,combo,'rows'));
    
    % compute probability
    probVal = cnt/M;
    prob(ii) = probVal;
    
    % store into the h matrix
    cellJ = num2cell(combo);
    h(cellJ{:}) = probVal;
    mapping{cellJ{:}} = combo;
end

end