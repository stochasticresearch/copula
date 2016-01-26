function [ prob, combos ] = computeEmpiricalDiscreteProb( X )
%COMPUTEEMPIRICALDISCRETEPROB - computes the empirical discrete probability
%of occurance of each of the possible outcomes in a multivariate discrete
%probability distribution
% Inputs:
%  X - input M x N matrix, with N being the dimensionality of the discrete
%      distribution, and M being the number of samples
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
        uniqueVals(nn) = length(unique(X(:,nn)));
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