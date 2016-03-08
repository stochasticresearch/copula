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

function [ U, X_i ] = frankcopularnd( M, N, alpha, varargin )
%FRANKCOPULARND Generates M samples from a Frank copula of dimensionality
%N, with parameter alpha
% Inputs:
%  M - the number of samples to generate
%  N - the dimensionality of the data
%  alpha - the dependency parameter of the Gumbel copula
%
% Outputs:
%  U - an M x N matrix of generated samples
%  X_i - an M x N matrix of intermediary random variables generated in the
%        creation of U

if(N<2)
    error('N must be atleast 2');
end

% Algorithm 1 described in both the SAS Copula Procedure, as well as the
% paper: "High Dimensional Archimedean Copula Generation Algorithm"
if(alpha<=0)
    error('For N>=3, alpha > 0 for the Frank Copula')
end

U = zeros(M,N);
X_i = rand(M,N);
for ii=1:M
    p = -1.0*expm1(-1*alpha);
    if (abs(1 - p) <= eps(p))
        % boundary protection
        p = 1 - eps;
    end
    vv = logserrnd(p, 1);

    % sample N independent uniform random variables
    x_i = X_i(ii,:);
    t = -1*log(x_i)./vv;
    U(ii,:) = -1.0*log1p( exp(-t)*expm1(-1.0*alpha))/alpha;
end

end % function