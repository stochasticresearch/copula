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

function [ X, llvec ] = genSynthData2( discreteType, continuousType, M )
%GENSYNTHDATA Generates synthetic data for testing HCBN, MTE, and CLG from
%the following BN (all arrows point downward):
%       A   B      
%        \ /
%         C 
% A, B, are multinomial discrete distributions, D is continuous
% random variables
%
% Inputs:
%  discreteType   - a cell array of 2 elements, where each element is a
%                   vector representing the probabilities in the
%                   multinomial distribution
%  continuousType - a string, options are: 
%                   'Gaussian' -> D will be Gaussian
%                   'Random'   -> D will be randomely picked 
%                                 continuous distributions
%  M - the number of realizations from this BN to generate
% Outputs:
%  X - the data matrix, col1->A, col2->B, col3->C

% let us define the following copula's:
%   F(A,B,C) -> C1   -> Gaussian Copula

a_dist = makedist('Multinomial','Probabilities',discreteType{1});
b_dist = makedist('Multinomial','Probabilities',discreteType{2});

% generate the copula random variables
Rho_C1 = [1 .4 .2; .4 1 -.8; .2 -.8 1];
U = copularnd('Gaussian', Rho_C1, M);

X12 =   [a_dist.icdf(U(:,1)) ...
          b_dist.icdf(U(:,2))];
% perform the inverse transform to generate the X data matrix
if(strcmpi(continuousType,'Gaussian'))
    rhoD = 0.6;
    X3 = norminv(U(:,3),0,rhoD);
else
%     X3 = betainv(U(:,3),2,5);
    % make it bimodal
    x = [normrnd(-2,0.3,M/2,1); normrnd(2,0.5,M/2,1)];
    [f,xi] = emppdf(x,0);
    F = empcdf(x,0);
    myObj = rvEmpiricalInfo(xi,f,F,0);
    X3 = zeros(M,1);
    for ii=1:M
        X3(ii) = myObj.icdf(U(ii,3));
    end
end

X = [X12 X3];
X = X(randperm(M),:);

% compute the likelihood of each point to the generative model
if(nargout>1)
    llvec = zeros(M,1);
    for ii=1:M
        xi = X(ii,:);
        u = [a_dist.cdf(xi(1)) b_dist.cdf(xi(2)) myObj.cdf(xi(3))];
        
        uuvec = zeros(1,length(u));
        for jj=1:length(u)
            if(abs(u(jj)-1)<=.001)
                uu = 0.999;
            elseif(abs(u(jj)-.001)<=0.001)
                uu = 0.001;
            else
                uu = u(jj);
            end
            uuvec(jj) = uu;
        end
        
        % here we are just computing f(x1,x2,x3) = f(x1)*f(x2)*f(x3)*c(F(x1),F(x2),F(x3)) and 
        % then taking the LOG
        llval = log( a_dist.pdf(xi(1))*b_dist.pdf(xi(2))*myObj.pdf(xi(3))* copulapdf('Gaussian', uuvec, Rho_C1) );
        fprintf('%s || %s copulapdf=%f\n', sprintf('%f,', xi), sprintf('%f,', uuvec), copulapdf('Gaussian', u, Rho_C1));
        llvec(ii) = llval;
    end
end

end