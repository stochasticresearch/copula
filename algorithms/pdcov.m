function [pdcov_val, pdcor_val, cov_pval, cor_pval] = pdcov( X, Y, Z, replicates )
%PDCORR computes the partial distance covariance between two random 
% variables X and Y, given Z. Rows represent the examples, and columns 
% the variables.  This implementation somewhat based on of dcorr.m here:
%  https://www.mathworks.com/matlabcentral/fileexchange/49968-dcorr--x--y--
% Reference: https://arxiv.org/pdf/1310.2926.pdf
% Inputs:
%  X - 
%  Y -
%  Z - 
% Outputs:
%  metric - the partial distance correlation of {X & Y}|Z
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

a = pdist2(X, X); % matrix euclidean distances
b = pdist2(Y, Y);
c = pdist2(Z, Z);

A = a - bsxfun(@plus,mean(a),mean(a,2))+mean(a(:));
B = b - bsxfun(@plus,mean(b),mean(b,2))+mean(b(:));
C = c - bsxfun(@plus,mean(c),mean(c,2))+mean(c(:));

AC = uproduct(A,C);
BC = uproduct(B,C);
CC = uproduct(C,C);

% TODO: divide by 0 protection here
c1 = AC./CC;
c2 = BC./CC;

P_xz = A-c1.*C;
P_yz = B-c2.*C;

pdcov_val = uproduct(P_xz,P_yz);
n = size(X,1);
teststat = n*pdcov_val;
den = sqrt(uproduct(P_xz, P_xz) * uproduct(P_yz, P_yz));
if(den>0)
    pdcor_val = teststat/(n*den);
else
    pdcor_val = 0;
end
    
if(nargout>2)
    
    if(nargin<4)
        replicates = 199;
    end

    % compute the p-value for the covariance value
    Tk_vec = zeros(1,replicates);
    parfor kk=1:replicates
        % permute Y
        I = randperm(size(P_xz,1));
        P_xz_perm = P_xz(I,:);
        P_xz_perm = P_xz_perm(:,I);
        % compute the statistic
        Tk_vec(kk) = uproduct(P_xz_perm,P_yz);
    end
    Tk_vec = n*Tk_vec;
    cov_pval = (1 + sum(Tk_vec>teststat))/(1+replicates);
    
    if(pdcor_val>0)
        nRootV = teststat/pdcor_val;
    else
        nRootV = 1;
    end
    teststat = pdcor_val;
    pdcor_reps = Tk_vec/nRootV;
    cor_pval   = (1 + sum(pdcor_reps > teststat))/(1 + replicates);
end

end

function [prod] = uproduct(X, Y)
n = size(X,1);
X(1:n+1:n*n) = 0;
Y(1:n+1:n*n) = 0;
zz = tril(X).*tril(Y);

sumval = sum(zz(:));
prod = 2*sumval/(n*(n-3));
end