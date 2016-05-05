function [u] = fixU(u, tol)
%FIXU - fixes the U query vector for querying copula densities
% u - a [M x D] matrix of the values in the unit hypercube to be fixed
% tol - (optional), how close the u values are allowed to get to the edge
%       of the hypercube, default is 0.01
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

if(nargin<=1)
    tol = 0.01;
end

u(abs(u-1)<=tol) = 1-tol;
u(abs(u-tol)<=tol) = tol;


end