function [ h, mapping, prob, combos ] = calcpmf( X )
%CALCPMF - computes the empirical discrete probability of occurance of each 
% of the possible outcomes in a multivariate discrete probability 
% distribution.  If the distribution is continuous, then it gets
% discretized according to the number of bins 
% Inputs:
%  X - input M x D matrix, with D being the dimensionality of the discrete
%      distribution, and M being the number of samples.  It is assumed that
%      each dimesion of X is discrete.  If this is not the case, use
%      discretizeRv to discretize the appropriate dimensions before using
%      this function!!!
%  Optional:
%   MAX_NUM_BINS - the maximum number of bins for which the contingency
%                  table should be generated
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

[M,D] = size(X);
if(D==1)
    error('D must be >= 2!');
end

% store the number of unique values for each dimension of the discrete
% distribution
hSize = zeros(1,D);
uniqueCell = cell(1,D);
for dd=1:D
    uniqueVals = unique(X(:,dd));
    hSize(dd) = length(uniqueVals);
    uniqueCell{dd} = uniqueVals;
end

combos = combvec(1:hSize(1),1:hSize(2));
for dd=3:D
    combos = combvec(combos, 1:hSize(dd));
end

h = zeros(hSize);
mapping = cell(hSize);

% for each combo, count the number of times they occur in the input data
% and compute the associated probability
combos = combos';
prob = zeros(size(combos,1),1);

for ii=1:size(combos,1)
    combo = combos(ii,:);
    mapVec = zeros(1,D);
    for dd=1:D
        unq = uniqueCell{dd};
        mapVec(dd) = unq(combo(dd));
    end

    % count the number of times this combo appears in the data set
    % WARNING: there may be some issues here b/c we are potentially
    % comparing double array's ..., not sure how ismember handles this?
    % the following test case works, but that doesn't mean all will ...
    %  ismember(1.111111111,1.111111111) --> 1
    cnt = sum(ismember(X,mapVec,'rows'));
    
    % compute probability
    probVal = cnt/M;
    prob(ii) = probVal;
    
    % store into the h matrix
    cellJ = num2cell(combo);
    h(cellJ{:}) = probVal;
    
    % store the mapping
    mapping{cellJ{:}} = mapVec;
end

end