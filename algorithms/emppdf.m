function [ f, xi ] = emppdf( x, isdiscrete, nbin )
%EMPPDF - Computes the empirical pdf
% Inputs:
%  x - the data from which the empirical cdf should be estimated
%  isdiscrete - 1 if x is discrete data, 0 if x is continuous data
% Outputs:
%  f - the empirical PDF
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

M = size(x,1);

if(isdiscrete)
    [~,xi] = empcdf( x, isdiscrete );
    diffOutput = diff(xi); binWidth = diffOutput(end);
    f = zeros(1,length(xi));
    idx = 1;
    for jj=1:length(xi)
        f(idx) = sum(x==xi(jj))/(M*binWidth);
        idx = idx + 1;
    end
else
    % make sure ksdensity only spreads mass over the support of the data
    pmVal = min(TOL,var(x)*TOL);
    supportVec = [min(x)-pmVal,max(x)+pmVal];
    
    if(nargin>2)
        [f,xi] = ksdensity(x,'Support',supportVec,'NumPoints',nbin);
    else
        [f,xi] = ksdensity(x,'Support',supportVec);
    end
end

f = f(:);
xi = xi(:);
f = f';
xi = xi';

f(f<0) = TOL;     % account for errors due to resample call
f(f>1) = 1-TOL;     % account for errors due to resample call

end