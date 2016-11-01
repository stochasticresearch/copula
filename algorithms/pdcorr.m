function [T0, pval] = pdcorr( X, Y, Z, replicates )
%PDCORR computes the partial distance correlation between two random 
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
%  pval - the associated p-value
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

T0 = pdcorr_impl(X,Y,Z);

if(nargout>1)
    if(nargin<4)
        replicates = 199;
    end

    % compute the p-value also
    Tk_vec = zeros(1,replicates);
    parfor kk=1:replicates
        % permute Y
        I = randperm(size(X,1));
        X_perm = X(I,:);
        % compute the statistic
        Tk_vec(kk) = pdcorr_impl(X_perm,Y,Z);
    end
    pval = (1 + sum(Tk_vec>=T0))/(1+replicates);
end

end

function [metric] = pdcorr_impl(X, Y, Z)

a = pdist2(X, X); % matrix euclidean distances
b = pdist2(Y, Y);
c = pdist2(Z, Z);

A = a - bsxfun(@plus,mean(a),mean(a,2))+mean(a(:));
B = b - bsxfun(@plus,mean(b),mean(b,2))+mean(b(:));
C = c - bsxfun(@plus,mean(c),mean(c,2))+mean(c(:));

n = size(X,1);

AA = sqrt( (sum(sum(A.*A))-sum(diag(A).*diag(A))) / (n*(n-3)) );
BB = sqrt( (sum(sum(B.*B))-sum(diag(B).*diag(B))) / (n*(n-3)) );
CC = sqrt( (sum(sum(C.*C))-sum(diag(C).*diag(C))) / (n*(n-3)) );

R_xy = ( (sum(sum(A.*B))-sum(diag(A).*diag(B))) / (n*(n-3)) ) / (AA*BB);
R_xz = ( (sum(sum(A.*C))-sum(diag(A).*diag(C))) / (n*(n-3)) ) / (AA*CC);
R_yz = ( (sum(sum(B.*C))-sum(diag(B).*diag(C))) / (n*(n-3)) ) / (BB*CC);

if(R_xz^2==0 || R_yz^2==0)
    metric = 0;
else
    metric = (R_xy - R_xz*R_yz) / (sqrt(1-R_xz^2)*sqrt(1-R_yz^2));
end


end