function [ taucj ] = ktaucj( X, Y )
%ktaucj - computes a rescaled version of Kendall's tau that preserves
%         the definition of Kendall's tau, but assures that in the 
%         scenario of perfect concordance or discordance for discrete
%         or hybrid datatypes, taucj achieves +/- 1 respectively
% Inputs:
%  X - first variable input.
%  Y - second variable input.
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
%  taucj - the rescaled version of Kendall's tau
%  
% NOTE: See "A Primer on Copulas for Count Data" by Genest and Neslehova
%       for more information on taucj
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

MAX_NUM_BINS = 25;
numUniqueX = length(unique(X));
numUniqueY = length(unique(Y));
if(numUniqueX>MAX_NUM_BINS)
    % discretize X to that number of bins for calculation of carley bounds
    X_discretized = discretizeRv(X, MAX_NUM_BINS)';
else
    X_discretized = X;
end
if(numUniqueY>MAX_NUM_BINS)
    % discretize X to that number of bins for calculation of carley bounds
    Y_discretized = discretizeRv(Y, MAX_NUM_BINS)';
else
    Y_discretized = Y;
end

% compute the normal kendall's tau
% ktau = corr(X,Y,'type','kendall');
ktau = ktaub([X Y],0.05,0);      % calculate tau w/out correcting for ties
                                 % the last 2 arguments are dont cares

% compute the pmf
h_XY = calcpmf( [X_discretized Y_discretized] );

% compute the Carley bounds
[ tau_C_plus, tau_C_minus ] = carleybounds( h_XY );
fprintf('tau_C_plus=%0.02f tau_C_minus=%0.02f\n', tau_C_plus, tau_C_minus);

% rescale the kendall's tau
if(ktau<0)
    taucj = -1*ktau/tau_C_minus;
elseif(ktau>0)
    taucj = ktau/tau_C_plus;
else
    taucj = 0;
end