function [ U ] = pobs( X, varargin )
%PSEUDOOBS Generates Pseudo-Observations from a given multivariate vector
%          See https://en.wikipedia.org/wiki/Copula_(probability_theory)#Empirical_copulas
%          for details
% Inputs:
%  X - an M x D multivariate input matrix, where M is the # of samples of
%      the multivariate joint distribution, and D is the dimensionality
% Optional Inputs:
%  varargin{1} - a string indicating which method to use.  By default, we
%                use ranks to generate the pseudo-observations, but an
%                empirical CDF method could also be used.  Valid values of
%                varargin{1} are: a.) rank
%                                 b.) ecdf
%  varargin{2} - if varargin{1} is ecdf, then varargin{2} must specify the
%                number of points to use in the ECDF
% Outputs:
%  U - an M x D matrix of the pseudo-observations
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
%* 
%**************************************************************************

[M,D] = size(X);

nVarargin = length(varargin);
method = 'rank';
correctionFlag = 0;
if(nVarargin>0)
    for ii=1:nVarargin
        if(strcmpi(varargin{ii},'method'))
            if(nVarargin>=ii+2)
                method = varargin{ii+1};
                if(strcmpi(method,'rank'))
                elseif(strcmpi(method,'ecdf'))
                    if(isnumeric(methodTypeCfg))
                        numEcdfPts = varargin{ii+2};
                    else
                        error('2nd argument for ECDF argmuent must be an integer (# of points)');
                    end
                end
                ii = ii + 2;        % we gobbled those arguments, so fast-foward
            else
                error('Not enough arguments provided for selecting pseudo-observations calculation method!');
            end
        elseif(strcmpi(varargin{ii},'correction'))
            if(nVarargin>=ii+3)
                correctionFlag = 1;
                sRhoMat = varargin{ii+1};
                sRhoMat
                alpha = varargin{ii+2};
                discreteIdxs = varargin{ii+3};
                if(~isnumeric(alpha))
                    error('alpha must be numeric values!');
                end
                ii = ii + 3;        % we gobbled those arguments, so fast-foward
            else
                error('Not enough arguments provided for selecting correction factor!');
            end
        end
    end
end

if(strcmpi(method,'rank'))
    U = tiedrank(X)/(M+1);      % scale by M+1 to mitigate boundary errors
elseif(strcmpi(method,'ecdf'))
    U = zeros(size(X));
    for nn=1:D
        isdiscrete = 0;     % TODO: this is suspect?
        domain = linspace(min(X(:,nn)),max(X(:,nn)),numEcdfPts);
        FX = ksdensity(X(:,nn), domain, 'function', 'cdf')';
        empInfoObj = rvEmpiricalInfo(domain, [], FX, isdiscrete);
        for mm=1:M
            U(mm,nn) = empInfoObj.cdf(X(mm,nn));
        end
    end
end

if(correctionFlag)
    
%     for jj=1:D
%         for ii=1:D
%             if(ii<jj)
%                 % check which of these indices is discrete (if any)
%                 pairwiseDiscreteIdxs = intersect( [ii jj], discreteIdxs );
%                 if(~isempty(pairwiseDiscreteIdxs))
%                     % get the pairwise correlation for this
%                     sRho = sRhoMat(ii,jj);
%                     % appropriately squeeze points
%                     if(length(pairwiseDiscreteIdxs)==2)
%                         % for now, we don't perform any squeezing if both
%                         % the indices are discrete
%                     else
%                         
%                     end
%                 end
%                 
%             end
%         end
%     end
    
    % TODO: generalize below for n-dimensional, right now only 2-d
    if(sRhoMat(1,2)>=0)
        projectionPoints = repmat(U(:,1), 1, 2);
    else
        projectionPoints = [U(:,1) 1-U(:,1)];
    end
    U_corrected = U + abs(sRhoMat(1,2))*alpha*(projectionPoints-U);
    U = U_corrected;
end

end



