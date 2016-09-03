function [ h_XY ] = gencontingencytable( X, Y, numXbin, numYbin )
%GENCONTINGENCYTABLE - Generates a discrete contingency table
% Inputs:
%  X - the first variable data
%  Y - the second variable data
%  Optional:
%    numXbin - the number of bins to divide X into.  If numXbin is 0, we
%              assume that X is a discrete RV and automatically calculate
%              the number of bins.  Default is numXbin to be 0.
%    numYbin - the number of bins to divide Y into.  If numYbin is 0, we
%              assume that Y is a discrete RV and automatically calculate
%              the number of bins.  Default is numYbin to be 0.
%    NOTE -- the numXbin and numYbin should only be specified for continuous
%            random variables
% Outputs:
%  h_XY - the calculated contingency table
%
% TODO: expand to more than 2 variables -- merge with hist_discrete
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

% TODO: use computeEmpiricalDiscreteProb to compute the combinations and
% the contingency table, the code below currently is repeated code which
% should be consolidated.

if(nargin<4)
    numYbin = 0;
end
if(nargin<3)
    numXbin = 0;
end

% some simple error checking
if(length(X)~=length(Y))
    error('X and Y must be the same dimension!');
end
M = length(X);

numUniqueX = numXbin;
numUniqueY = numYbin;

if(numUniqueX==0)
	numUniqueX = unique(X);
else
    % bin the X data according to the desired number of points
    X = discretizeRv( X, numUniqueX );
end

if(numUniqueY==0)
	numUniqueY = unique(Y);
else
    % bin the Y data according to the desired number of points
    Y = discretizeRv( Y, numUniqueY );
end

h_XY = zeros(numUniqueX, numUniqueY);
% for each X/Y combo, we compute an empirical probability
for xx=1:numUniqueX
    for yy=1:numUniqueY
        xMatchIdx = find(xx==X);
        yMatchIdx = find(yy==Y);
        matchIdxs = intersect(xMatchIdx,yMatchIdx);
        h_XY(xx,yy) = length(matchIdxs)/M;
    end
end