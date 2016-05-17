function [corrcoef, tieCorrectedCoefMat] = ndcorr(X, type)
%NDCORR - compute the n-dimensional rank based correlation coefficient
%based on the standard generalization given by:
%      Friedrich Schmid and Rafael Schmidt. 
%      Multivariate extensions of Spearman's Rho and Related Statistics. 
%      Statistics and Probability Letters, 77(4):407 ? 416, 2007.
% Inputs:
%  X - an M x D data matrix for which the generalized m-dimensional rank
%      based correlation coefficient is desired.  M is the # of samples, D
%      is the dimensionality
%  type - should be either 'spearman' or 'kendall'
% Outputs:
%  corrcoef - the average of all the pair-wise correlations
%  pairwiseCoeffsVec - a vector of all the pairwise correlations
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

[~,D] = size(X);
numUniquePossibilities = nchoosek(D,2);
pairwiseCoeffsVec = zeros(1,numUniquePossibilities);

sRhoTieCorrectedVecIdx = 1;
tieCorrectedCoefMat = eye(D,D);
if(strcmpi(type,'spearman'))
    for jj=1:D
        for ii=1:D
            if(ii<jj)
                matrixSlice = [X(:,ii) X(:,jj)];
                sRho = corr(X(:,ii),X(:,jj),'type','spearman');
                sRhoTieCorrected = sRho/(12*sqrt(prod(var(matrixSlice))));
                pairwiseCoeffsVec(sRhoTieCorrectedVecIdx) = sRhoTieCorrected;
                
                % make a symmetric "correlation" matrix
                tieCorrectedCoefMat(ii,jj) = sRhoTieCorrected;
                tieCorrectedCoefMat(jj,ii) = sRhoTieCorrected;
                
                sRhoTieCorrectedVecIdx = sRhoTieCorrectedVecIdx + 1;
            end
        end
    end
elseif(strcmpi(type,'kendall'))
    error('Currently unimplemented!');
else
    error('Unknown type argument specified!');
end

corrcoef = mean(pairwiseCoeffsVec);

end