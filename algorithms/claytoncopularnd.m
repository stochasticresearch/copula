function [ U ] = claytoncopularnd( M, D, alpha )
%CLAYTONCOPULARND Generates M samples from a Clayton copula of 
% dimensionality D, with parameter alpha
% Inputs:
%  M - the number of samples to generate
%  N - the dimensionality of the data
%  alpha - the dependency parameter of the Clayton copula
%
% Outputs:
%  U - an M x N matrix of generated samples
%
% Acknowledgements:
% This code modeled after the function definitions in the paper:
% Estimators for Archimedean copulas in high dimensions, 
% by Marius Hofert1, Martin Mächler, and Alexander J. McNeil AND
% the R code in cop_objects.R in the "copula" package in R, found at:
% https://cran.r-project.org/web/packages/copula/
%
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

if(D<2)
    error('N must be atleast 2');
end
if(alpha<0)
    error('Clayton copula parameter must be between [0, inf)')
end

% Algorithm 1 described in both the SAS Copula Procedure, as well as the
% paper: "High Dimensional Archimedean Copula Generation Algorithm"
U = zeros(M,D);
for ii=1:M
    shape = 1.0/alpha;
    scale = 1;
    v = gamrnd(shape, scale);

    % sample N independent uniform random variables
    x_i = rand(1,D);
    t = -1*log(x_i)./v;

    tmp = 1.0 + t;
    U(ii,:) = tmp.^(-1.0/alpha);
end

end % function