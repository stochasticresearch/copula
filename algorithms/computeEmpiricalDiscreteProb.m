%******************************************************************************
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

function [ prob, combos ] = computeEmpiricalDiscreteProb( X )
%COMPUTEEMPIRICALDISCRETEPROB - computes the empirical discrete probability
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
%  prob - a vector containing the probabilities for each of the
%         combinations, referenced in the combo vector
%  combo - the combo for which the probability was calculated

[M,N] = size(X);

if(N==1)
    numUniqueVals = length(unique(X));
    combos = 1:numUniqueVals;
else
    % store the number of unique values for each dimension of the discrete
    % distribution
    uniqueVals = zeros(1,N);
    for nn=1:N
        uniqueVals(nn) = max(X(:,nn));
    end

    combos = combvec(1:uniqueVals(1),1:uniqueVals(2));
    for nn=3:N
        combos = combvec(combos, 1:uniqueVals(nn));
    end
end

% for each combo, count the number of times they occur in the input data
% and compute the associated probability
combos = combos';
prob = zeros(size(combos,1),1);

for ii=1:size(combos,1)
    combo = combos(ii,:);

    % count the number of times this combo appears in the data set
    cnt = 0;
    for mm=1:M
        if(isequal(combo,X(mm,:)))
            cnt = cnt + 1;
        end
    end

    % compute probability
    prob(ii) = cnt/M;
end



end