function [ y ] = gumbelcopulapdf( u, alpha, varargin )
%GUMBELCOPULAPDF Computes the PDF of the Gumbel Copula for N>=2
% Inputs:
%  u - an M x N matrix of all the points over which to compute the Gumbel
%      copula PDF, M is the number of points in a unit-hypercube of
%      dimension N
%  alpha - the dependency parameter of the Gumbel copula.  
%      Please NOTE!! There is a difference in terminology between Mathworks
%      and R, in R (and seemingly in academic literature), theta is used as
%      the dependency parameter, and alpha = 1/theta.  Mathworks seems to
%      like to define the dependency parameter as alpha.  In this code, we
%      stick to the (unfortunate) Mathworks convention of using alpha as
%      the dependency parameter!!, 
%      Consequently we assign one_div_alpha = 1/alpha
% Optional Inputs:
%  varargin{1} - if 0, then compute pdf directly, 
%                else, compute log of pdf at specified value
%
% Outputs:
%  y - the value of the Gumbel copula density at the specified value in
%      the unit hypercube
%
% Acknowledgements:
% This code modeled after the function definitions in the paper:
% Estimators for Archimedean copulas in high dimensions, 
% by Marius Hofert1, Martin Mächler, and Alexander J. McNeil AND
% the R code in cop_objects.R in the "copula" package in R, found at:
% https://cran.r-project.org/web/packages/copula/
%
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

nVarargs = length(varargin);
if(nVarargs==0)
    wantLog = 0;
else
    if(isnumeric(varargin{1}))
        wantLog = varargin{1};
    else
        warning('Invalid varargin{1} for frankcopulapdf, defaulting to LOG=FALSE');
        wantLog = 0;
    end
end

[~,N] = size(u);

mlu = -log(u); % -log(u)
lmlu = log(mlu); % log(-log(u))

one_div_alpha = 1/alpha;
lx = one_div_alpha*log( sum(iPsi_gumbel(u, alpha),2) );     % might have to turn LOG for iPsi_gumbel call
                                                            % WARNING!! This computation of lx may not be numerically stable!! :(
                                                            % TODO: need to investigate                            
ls = polyG(lx, one_div_alpha, N, 1)-N*lx/one_div_alpha;

lnC = -exp(lx); % - t^alpha
y = lnC + N*log(alpha) + sum((alpha-1)*lmlu + mlu,2) + ls;

if(~wantLog)
    y = exp(y);
end

% respect the grounded property of the copula manually, in case the user
% requested boundary copula pdf values
idxs = ~all(u,2);
y(idxs) = 0;
% TODO: figure out better way to deal w/ this below.  I think according to
% Matlab, all edges of any copula density are 0 ... I don't think R
% follows this, so we need to decide which is more right.
y(isnan(y)) = 0;

end

function [ y ] = iPsi_gumbel( u, alpha, varargin )

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

if(useLog) 
    y = alpha*log(-log(u));
else
    y = (-log(u+0)).^alpha;
end

end