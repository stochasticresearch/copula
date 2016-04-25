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

function [ U ] = pseudoobs( X, varargin )
%PSEUDOOBS Generates Pseudo-Observations from a given multivariate vector
%          See https://en.wikipedia.org/wiki/Copula_(probability_theory)#Empirical_copulas
%          for details
% Inputs:
%  X - an M x N multivariate input matrix, where M is the # of samples of
%      the multivariate joint distribution, and N is the dimensionality
% Optional Inputs:
%  varargin{1} - a string indicating which method to use.  By default, we
%                use ranks to generate the pseudo-observations, but an
%                empirical CDF method could also be used.  Valid values of
%                varargin{1} are: a.) rank
%                                 b.) ecdf
%  varargin{2} - if varargin{1} is ecdf, then varargin{2} must specify the
%                number of points to use in the ECDF
% Outputs:
%  U - an M x N matrix of the pseudo-observations

[M,N] = size(X);

nVarargin = length(varargin);
if(nVarargin==0)
    method = 'rank';
else
    method = varargin{1};
    if(strcmpi(method,'rank'))
    elseif(strcmpi(method,'ecdf'))
        if(nVarargin~=2)
            error('ecdf method requires a second argument of the number of ecdf points to generate!');
        end
        numEcdfPts = varargin{2};
    else
        error('varargin{1} must be either rank or ecdf!');
    end
end

if(strcmpi(method,'rank'))
    U = tiedrank(X)/(M+1);      % scale by M+1 to mitigate boundary errors
elseif(strcmpi(method,'ecdf'))
    U = zeros(size(X));
    for nn=1:N
        domain = linspace(min(X(:,nn)),max(X(:,nn)),numEcdfPts);
        FX = ksdensity(X(:,nn), domain, 'function', 'cdf')';
        empInfoObj = rvEmpiricalInfo(domain, [], FX);
        for mm=1:M
            U(mm,nn) = empInfoObj.cdf(X(mm,nn));
        end
    end
end