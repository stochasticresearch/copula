function [ taucj ] = ktaucj( X, Y, numXbin, numYbin )
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


if(nargin<4)
    numYbin = 0;
end
if(nargin<3)
    numXbin = 0;
end

% compute the normal kendall's tau
ktau = corr(X,Y,'type','kendall');

% compute the contingency table
h_XY = gencontingencytable( X, Y, numXbin, numYbin );

% compute the Carley bounds
[ tau_C_plus, tau_C_minus ] = carleybounds( h_XY );

% rescale the kendall's tau
if(ktau<0)
    taucj = -1*ktau/tau_C_minus;
elseif(ktau>0)
    taucj = ktau/tau_C_plus;
else
    taucj = 0;
end