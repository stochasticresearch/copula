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

function [ U ] = empcopularnd( c, M, varargin )
%EMPCOPULARND Sample from any given copula density.
%
% Inputs:
%  c - The copula density, provided as a [K x K x ... x K] matrix, where
%      the dimensions are [u_1,u_2, ..., u_D]
%  M - The number of samples to generate from this copula
%
% Outputs
%  U - an M x D vector of samples from the copula, given by the integral
%      of the provided copula density.

K = size(c,1);
uu = linspace(0,1,K);
D = length(size(c));

c_d = cell(1,D);    % c_d{i} stores c_i(u_1,...,u_i) = integral(c,du_{i+1} ... du_{D})
                    % For more information, refer to: Analysis and
                    % Generation of Random Vectors with Copulas, by Johann
                    % Strelen and Feras Nassaj
c_d{D} = c;
% integrate out each dimension successively
for ii=D-1:-1:1
    c_d{ii} = squeeze(sum(c_d{ii+1},ii+1));
end

U = rand(M,D);
genStartIdx = 2;
D_parents = 1;
if(~isempty(varargin))
    U_parents = varargin{1};
    if(size(U_parents,1)~=M)
        error('U_parents must be of length M');
    end
    D_parents = size(U_parents,2);
    if(D_parents>=D)
        error('U_parents must be of less dimensions than D!');
    end
    U(:,1:D_parents) = U_parents;
    genStartIdx = D_parents+1;
end

for ii=1:M
    idxVec = zeros(1,D); 
    for kk=1:D_parents
        u_j = U(ii,kk);
        idxVec(kk) = findClosest(uu,u_j);
    end
    for jj=genStartIdx:D
        C_d_num = getMarginalIntegral(c_d{jj},idxVec(1:jj-1),K);
        C_d_den = c_d{jj-1}(getLinearIdx(idxVec(1:jj-1),K));

        C_d = C_d_num./C_d_den;
        % perform numerical inverse
        u_j = uu(findClosest(C_d,U(ii,jj)));
        U(ii,jj) = u_j;

        idxVec(jj) = findClosest(uu,u_j);
    end
end

end

function [colLinearIdx] = getLinearIdx(arrIdx,K)
colLinearIdx = arrIdx(:,1);
multiplyFactor = K;
for ii=2:size(arrIdx,2)
    colLinearIdx = colLinearIdx + (arrIdx(:,ii)-1)*multiplyFactor;
    multiplyFactor = multiplyFactor*K;
end
end

function [y] = getMarginalIntegral(c, valVec, K)
linearIdxs = zeros(K, length(valVec)+1);
linearIdxs(:,length(valVec)+1) = 1:K;
linearIdxs(:,1:length(valVec)) = repmat(valVec,K,1);
linearIdxs = getLinearIdx(linearIdxs,K);
y = c(linearIdxs);
y = cumsum(y);      % integrate
end

function [idx] = findClosest(vec, val)
tmp = abs(val-vec);
[~,idx] = min(tmp);
end