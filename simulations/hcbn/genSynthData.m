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

function [ X, llvec ] = genSynthData( discreteType, continuousType, M )
%GENSYNTHDATA Generates synthetic data for testing HCBN, MTE, and CLG from
%the following BN (all arrows point downward):
%       A   B      
%      / \ / \
%     C   D   E
% A, B, are multinomial discrete distributions, C, D, E are continuous
% random variables
%
% Inputs:
%  discreteType   - a cell array of 2 elements, where each element is a
%                   vector representing the probabilities in the
%                   multinomial distribution
%  continuousType - a string, options are: 
%                   'Gaussian' -> C,D,E will be Gaussian
%                   'Random'   -> C,D,E will be randomely picked 
%                                 continuous distributions
%  M - the number of realizations from this BN to generate
% Outputs:
%  X - the data matrix, col1->A, col2->B, col3->C, col4 -> D, col5-> E

% let us define the following copula's:
%   F(A,C)   -> C1   -> Frank Copula
%   F(A,D,B) -> C2   -> Gaussian Copula
%   F(B,E)   -> C3   -> Clayton Copula

a_dist = makedist('Multinomial','Probabilities',discreteType{1});
b_dist = makedist('Multinomial','Probabilities',discreteType{2});

% generate the copula random variables
Rho_C2 = [1 .4 .2; .4 1 -.8; .2 -.8 1];
U_C2 = copularnd('Gaussian', Rho_C2, M);

c1_alpha = 2; p = rand(M,1);
U_C1_1 = U_C2(:,1);
U_C1_2 = -log((exp(-c1_alpha.*U_C1_1).*(1-p)./p + exp(-c1_alpha))./(1 + exp(-c1_alpha.*U_C1_1).*(1-p)./p))./c1_alpha;
U_C1 = [U_C1_1 U_C1_2];

c3_alpha = 6; p = rand(M,1);
U_C3_1 = U_C2(:,3);
U_C3_2 = U_C3_1.*(p.^(-c3_alpha./(1+c3_alpha)) - 1 + U_C3_1.^c3_alpha).^(-1./c3_alpha);
U_C3 = [U_C3_1 U_C3_2];

U = [U_C1(:,2) U_C2 U_C3(:,2)];

X12 =   [a_dist.icdf(U(:,1)) ...
          b_dist.icdf(U(:,2))];
% perform the inverse transform to generate the X data matrix
if(strcmpi(continuousType,'Gaussian'))
    % assume 0 mean, pick random std deviations 0.1 and 3
    rhoC = 0.6;
    rhoD = 2.3;
    rhoE = 0.1;
    X345 = [norminv(U(:,3),0,rhoC) ...
            norminv(U(:,4),0,rhoD) ...
            norminv(U(:,5),0,rhoE)];
else
    % pick random one parameter distributions
%     distChoices = randsample(5,3,true)';
    % choice to rv type mapping:
    %  1 -> Gaussian
    %  2 -> Exponential
    %  3 -> Gamma
    %  4 -> Beta
    %  5 -> Uniform
    
    % create X3, X4, X5
    X345 = zeros(M,3);
    
    X345(:,1) = betainv(U(:,3),2,5);
    X345(:,2) = gaminv(U(:,4),2,2);
    X345(:,3) = unifinv(U(:,5),2,6);
    
end

X = [X12 X345];
X = X(randperm(M),:);       % shuffle the data

for ii=1:M
    llvec = zeros(M,1);
	xi = X(ii,:);
	
	u1 = [a_dist.cdf(xi(1)) betacdf(xi(3),2,5)];
	u2 = [a_dist.cdf(xi(1)) b_dist.cdf(xi(2)) gamcdf(xi(4),2,2)];
	u2_parents = u2(1:2);
	u3 = [b_dist.cdf(xi(2)) unifcdf(xi(5),2,6)];

	% compute f(x1,x2,x3,x4,x5) = f(x1)*f(x2)*f(x3)*f(x4)*f(x5)*c1(F(x1),F(x3))*c2(F(x1),F(x2),F(x4))*c3(F(x2),F(x5))/c2_parents(F(x1),F(x2)) 
    % and then taking the LOG.  This equation is valid for the BN structure defined above.
	
	indepProbs = a_dist.pdf(xi(1))*b_dist.pdf(xi(2))*betapdf(xi(3),2,5)*gampdf(xi(4),2,2)*unifpdf(xi(5),2,6);
	numeratorCopulas = frankcopulapdf(u1,c1_alpha)*copulapdf('Gaussian',u2,Rho_C1)*claytoncopulapdf(u3,c3_alpha);
	denominatorCopula = copulapdf('Gaussian',u2_parents,Rho_C2(1:2,1:2));
	
	llval = log( indepProbs*numeratorCopulas/denominatorCopula );
	llvec(ii) = llval;
end

end