function [corrcoef] = ndcorr(X, type)
%NDCORR - compute the n-dimensional rank based correlation coefficient
%based on the standard generalization given by:
%      Friedrich Schmid and Rafael Schmidt. 
%      Multivariate exten- sions of spearman?s rho and related statistics. 
%      Statistics and Probability Letters, 77(4):407 ? 416, 2007.
% Inputs:
%  X - an M x D data matrix for which the generalized m-dimensional rank
%      based correlation coefficient is desired.  M is the # of samples, D
%      is the dimensionality
%  type - should be either 'spearman' or 'kendall'
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
% weight = 1/nchoosek(D,2);
sRhoTieCorrectedVec = zeros(1,nchoosek(D,2));

if(strcmpi(type,'spearman'))
    corrcoef = 0;
    for jj=1:D
        for ii=1:D
            if(ii<jj)
                matrixSlice = [X(:,ii) X(:,jj)];
                sRho = corr(X(:,ii),X(:,jj),'type','spearman');
                sRhoTieCorrected = sRho/(12*sqrt(prod(var(matrixSlice))));
%                 corrcoef = corrcoef + weight*sRhoTieCorrected;
                sRhoTieCorrectedVec(sRhoTieCorrectedVecIdx) = sRhoTieCorrected;
            end
        end
    end
elseif(strcmpi(type,'kendall'))
    error('Currently unimplemented!');
else
    error('Unknown type argument specified!');
end

% if all sRho's are similar in magnitude, we assume that this is a
% one-parameter archimedean copula so we just average all the
% sRhoTieCorrected values.  
% if there are sRho's that are significant then we assume that we have a 
% multiparameter copula, so we just return the vector as-is in the form of 
% a correlation matrix. Whether that correlation matrix is then used for a 
% Gaussian copula, or something else is upto the user!

% TODO: we need a statistically significant test for testing whether the
% sRho's magnitude differences are significant enough to warrant deeming
% this a multi-parameter copula rather than a single parameter copula.

end