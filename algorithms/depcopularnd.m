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

function [ U ] = depcopularnd( Z, D, family, depParam )
%DEPCOPULARND Aids in the generation of generating random numbers from
%graphical models with common parents.
% Inputs:
%  Z - an [M x 1] vector of the parent's common random variable.  For the
%  Gaussian family, this should be U_i(:,1), and for the 
%  D - the dimensionality of the data to generate.  D can only be 2 for
%      Clayton/Frank/Gumbel families currently
%  Clayton/Frank/Gumbel families, this should be X_i(:,1)
%  family - the desired family, can be either:
%           Gaussian
%           Clayton
%           Frank
%           Gumbel
%  depParam - for the Gaussian family, this must be a correlation matrix,
%             for Clayton/Frank/Gumbel, this should be the alpha value
% Outputs:
%  U - the copula random variables

M = size(Z,1);
family_LC = lower(family);
switch family_LC
    case 'gaussian'
        UU = norminv(Z);
        VV = mvnrnd(zeros(1,D-1),eye(D-1),M);
        U_uncorrelated = [UU VV];
        A = chol(depParam);
        U = normcdf(U_uncorrelated*A);
    case 'gumbel'
        error('Not yet implemented :(');
        if(D>2)
            error('D>2 unsupported for Gumbel copula!');
        end
    case 'frank'
        if(D>2)
            error('D>2 unsupported for Frank copula!');
        end
        
        p = rand(M,1);
        alpha = depParam;
        u1 = Z;
        if abs(alpha) > log(realmax)
            u2 = (u1 < 0) + sign(alpha).*u1; % u1 or 1-u1
        elseif abs(alpha) > sqrt(eps)
            u2 = -log((exp(-alpha.*u1).*(1-p)./p + exp(-alpha))./(1 + exp(-alpha.*u1).*(1-p)./p))./alpha;
        else
            u2 = p;
        end
        U = [u1 u2];
        
    case 'clayton'
        if(D>2)
            error('D>2 unsupported for Clayton copula!');
        end
        
        p = rand(M,1);
        alpha = depParam;
        u1 = Z;
        if alpha < sqrt(eps)
            u2 = p;
        else
            u2 = u1.*(p.^(-alpha./(1+alpha)) - 1 + u1.^alpha).^(-1./alpha);
        end
        U = [u1 u2];
        
    otherwise
        error('Unsupported Copula Family Type!');
end

end