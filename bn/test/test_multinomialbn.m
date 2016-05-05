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

% A test script for multinomialbn
clear;
clc;
bntPath = '../bnt'; addpath(genpath(bntPath));

% setup the graphical model
aa = 1; bb = 2; cc = 3;
D = 3;
dag = zeros(D,D);
dag(aa,cc) = 1;
dag(bb,cc) = 1;
discreteNodes = [aa bb];

% generate the data
M = 2000;
Rho = [1 .4 .2; .4 1 -.8; .2 -.8 1];
U = copularnd('Gaussian', Rho, M);

probsA = [0.25 0.25 0.25 0.25];
probsB = [0.25 0.25 0.25 0.25];
a_dist = makedist('Multinomial','Probabilities', probsA);
b_dist = makedist('Multinomial','Probabilities', probsB);
X = zeros(size(U));
X(:,1) = a_dist.icdf(U(:,1));
X(:,2) = b_dist.icdf(U(:,2));
X(:,3) = norminv(U(:,3));       % standard normal distribution for now ...

X_train = X(1:M/2,:);
X_test = X(M/2+1:end,:);

multinomialbnObj = multinomialbn(X_train, discreteNodes, dag);
llVal = multinomialbnObj.dataLogLikelihood(X_test)

% for each combination of values for nodes A and B, compute what the
% "actual" distribution looks like, and what the discretized distribution
% looks like, hopefully they are similar :D
for aVal=1:length(a_dist.Probabilities)
    for bVal=1:length(b_dist.Probabilities)
        X_subset = [];
        combo = [aVal bVal];
        for mm=1:M/2
            if(isequal(X_test(mm,1:2),combo))
                X_subset = [X_subset X_test(mm,3)];
            end
        end
        
        % estimate the empirical pdf of X_subset
        [actual_epdf, actual_domain] = ksdensity(X_subset);
        estBnParams = multinomialbnObj.bnParams{3};
        numCombos = length(estBnParams);
        for comboIdx=1:numCombos
            comboInfo = estBnParams{comboIdx};
            if(isequal(comboInfo.combo, combo))
                estimated_info = comboInfo.empInfo;
                break;
            end
        end
        
        f = figure(1);
        plot(actual_domain, actual_epdf);
        hold on;
        stem(estimated_info.domain, estimated_info.density);
        grid on;
        title(sprintf('A=%d B=%d continuousINT=%f discreteINT=%f', ...
            aVal, bVal,trapz(actual_domain,actual_epdf), ...
            trapz(estimated_info.domain, estimated_info.density) ));
%         estimated_info.distribution
%         trapz(actual_domain,actual_epdf)
        pause;
        clf(f);
        
    end
end