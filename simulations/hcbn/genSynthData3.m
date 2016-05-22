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

function [ X, llvec ] = genSynthData3( discreteType, continuousType, M )
%GENSYNTHDATA Generates synthetic data for testing HCBN, MTE, and CLG from
%the following BN (all arrows point downward):
%       A-->B        
% A, is a multinomial discrete distributions, B is continuous multimodal
%
% Inputs:
%  discreteType   - a cell array of 1 elements, where each element is a
%                   vector representing the probabilities in the
%                   multinomial distribution
%  continuousType - a string, options are: 
%                   'Gaussian' -> B will be Gaussian
%                   'Random'   -> B will be a multimodal gaussian
%  M - the number of realizations from this BN to generate
% Outputs:
%  X - the data matrix, col1->A, col2->B

% let us define the following copula's:
%   F(A,B,C) -> C1   -> Gaussian Copula

a_dist = makedist('Multinomial','Probabilities',discreteType{1});

% generate the copula random variables
D = 2; alpha = 4;
% U = gumbelcopularnd(M, D, alpha);
U = copularnd('Gumbel', alpha, M);
X = zeros(M,2);

X(:,1) = a_dist.icdf(U(:,1));
% perform the inverse transform to generate the X data matrix
if(strcmpi(continuousType,'Gaussian'))
    rhoD = 0.6;
    X(:,2) = norminv(U(:,2),0,rhoD);
else
    % make it bimodal
    x = [normrnd(-2,0.3,M/2,1); normrnd(2,0.5,M/2,1)];
    [f,xi] = emppdf(x,0);
    F = empcdf(x,0);
    myObj = rvEmpiricalInfo(xi,f,F,0);
    for ii=1:M
        X(ii,2) = myObj.icdf(U(ii,2));
    end
%     X(:,2) = unifinv(U(:,2),-2,2);
end

X = X(randperm(M),:);

% compute the likelihood of each point to the generative model
if(nargout>1)
    llvec = zeros(M,1);
    for ii=1:M
        xi = X(ii,:);
%         u = [a_dist.cdf(xi(1)) myObj.cdf(xi(2))];
        u = [a_dist.cdf(xi(1)) unifcdf(xi(2),-2,2)];
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
        
        % here we are just computing f(x1,x2) = f(x1)*f(x2)*c(F(x1),F(x2)) and then taking
        % the LOG
%         llval = log( a_dist.pdf(xi(1))*myObj.pdf(xi(2))* frankcopulapdf(uuvec, alpha) );
        llval = log( a_dist.pdf(xi(1))*unifpdf(xi(2),-2,2)* frankcopulapdf(uuvec, alpha) );
%         fprintf('%s || %s >> copulapdf=%0.02f >> f(x1)=%0.02f f(x2)=%0.02f ll=%0.02f\n', ...
%             sprintf('%f,', xi), sprintf('%f,', uuvec), frankcopulapdf(uuvec, alpha), ... 
%             a_dist.pdf(xi(1)), myObj.pdf(xi(2)), llval );
        llvec(ii) = llval;
    end
end

end