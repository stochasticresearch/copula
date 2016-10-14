function [metric] = lnn_CD(x, y, z)
%MUTUALINFOCD - a function which uses novel mutual information estimators
%developed by other entities in order to compute a metric of conditional
%dependence by conditioning x and y onto z.
% WARNING!! - This function offloads processing onto Python!
% Inputs:
%   x - [1 x n] input vector x
%   y - [1 x n] input vector y
%   z - [1 x n] input vector z, which x and y are conditioned upon
% Outputs
%   cd - conditional dependence metric
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

% concatenate x,y,z and create split sizes so that we can disambiguate in
% python

I_XZ = py.lnn.mi_matlab(x, z);
I_X_YZ = py.lnn.mi_matlab(x, [y z]);
I_X_Y_givenZ = I_X_YZ - I_XZ;

if(I_X_Y_givenZ<0) I_X_Y_givenZ = 0;

% convert this into a informational coefficient of correlation
metric = sqrt(1-exp(-2*I_X_Y_givenZ));

end