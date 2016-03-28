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

function [ X ] = genSynthData2( discreteType, continuousType, M )
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
    X3 = betainv(U(:,3),2,5);
end

X = [X12 X3];

end