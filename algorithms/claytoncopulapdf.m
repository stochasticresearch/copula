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

function [ y ] = claytoncopulapdf( u, alpha, varargin )
%GUMBELCOPULAPDF Computes the PDF of the Clayton Copula for N>=2
% Inputs:
%  u - an M x N matrix of all the points over which to compute the Clayton
%      copula PDF, M is the number of points in a unit-hypercube of
%      dimension N
%  alpha - the dependency parameter of the Clayton copula.  
%      Please NOTE!! There is a difference in terminology between Mathworks
%      and R, in R (and seemingly in academic literature), theta is used as
%      the dependency parameter, and alpha = 1/theta.  Mathworks seems to
%      like to define the dependency parameter as alpha.  In this code, we
%      stick to the (unfortunate) Mathworks convention of using alpha as
%      the dependency parameter!!, 
%      Consequently we assign 1_div_alpha = 1/alpha
% Optional Inputs:
%  varargin{1} - if 0, then compute pdf directly, 
%                else, compute log of pdf at specified value
% Outputs:
%  y - the value of the Clayton copula density at the specified value in
%      the unit hypercube
%
% Acknowledgements:
% This code modeled after the function definitions in the paper:
% Estimators for Archimedean copulas in high dimensions, 
% by Marius Hofert1, Martin Mächler, and Alexander J. McNeil AND
% the R code in cop_objects.R in the "copula" package in R, found at:
% https://cran.r-project.org/web/packages/copula/

nVarargs = length(varargin);
if(nVarargs==0)
    wantLog = 0;
else
    if(isnumeric(varargin{1}))
        wantLog = varargin{1};
    else
        warning('Invalid varargin{1} for claytoncopulapdf, defaulting to LOG=FALSE');
        wantLog = 0;
    end
end


[~,N] = size(u);

% compute log of the pdf
lu = sum(log(u),2);
t = sum(iPsi_clayton(u, alpha), 2);
y = sum(log1p(alpha*(0:N-1))) - (1+alpha)*lu - (N + 1.0/alpha)*log1p(t);
% exponentiate the result if the desired value is not the log version
if(~wantLog)
    y = exp(y);
end

% respect the grounded property of the copula manually, in case the user
% requested boundary copula pdf values
idxs = ~all(u,2);
y(idxs) = 0;

end

function [ y ] = iPsi_clayton( u, alpha, varargin )

nVarargs = length(varargin);
if(nVarargs==0)
    useLog = 0;
else
    if(isnumeric(varargin{1}))
        useLog = varargin{1};
    else
        warning('Invalid varargin{1} for iPsi_clayton, defaulting to LOG=FALSE');
        useLog = 0;
    end
end

y = u.^(-alpha) - 1;
if(useLog)
    y = log(y);
end

end