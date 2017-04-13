function [ F, xi ] = empcdf( x, isdiscrete, nbin )
%EMPCDF - Computes the empirical cdf
% Inputs:
%  x - the data from which the empirical cdf should be estimated
%  isdiscrete - 1 if x is discrete data, 0 if x is continuous data
% Outputs:
%  F - the empirical CDF
%  xi - the domain over which F is defined
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

TOL = 1e-3;

if(isdiscrete)
    [F,x] = ecdf(x);
    diffOutput = diff(x);
    
    if(length(diffOutput)<2)
        subtractAmt = 1;
    else
        subtractAmt = diffOutput(2);
    end
    
    xi = x; xi(1) = xi(2)-subtractAmt;
else
    [F, xi] = ksdensity(x,'function','cdf');
    if(nargin>2)
        % resample so that we get the desired number of bins
        if(nbin>100 || mod(100,nbin))
            error('nbin unsupported -- must be divisible by 100 and <=100');
        end
        F = downsample(F,100/nbin);
        xi = downsample(xi,100/nbin);
    end
end

F = F(:);
xi = xi(:);
F = F';
xi = xi';

F(F<0) = TOL;   % ensure distribution is bounded below by 0
F(F>1) = 1-TOL; % ensure distribution is bounded above by 1

end