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

function [ U ] = gumbelcopularnd( M, D, alpha )
%GUMBELCOPULARND Generates M samples from a Gumbel copula of dimensionality
%D, with parameter alpha
% Inputs:
%  M - the number of samples to generate
%  N - the dimensionality of the data
%  alpha - the dependency parameter of the Gumbel copula
%
% Outputs:
%  U - an M x N matrix of generated samples
%  X_i - an M x N matrix of intermediary random variables generated in the
%        creation of U

if(D<2)
    error('N must be atleast 2');
end
if alpha < 1
    error('Gumbel copula parameter must be between [1, inf)');
end

% Algorithm 1 described in both the SAS Copula Procedure, as well as the
% paper: "High Dimensional Archimedean Copula Generation Algorithm"
U = zeros(M,D);
for ii=1:M
    a  = 1.0/alpha;
    b  = 1;
    g  = cos(pi/(2.0*alpha)).^alpha;
    d  = 0;
    pm = 1;
    vv = rstable1(1,a,b,g,d,pm);

    % sample N independent uniform random variables
    x_i = rand(1,D);
    t = -1*log(x_i)./vv;

    U(ii,:) = exp(-1*(t.^(1.0/alpha)));
end

end % function