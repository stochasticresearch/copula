function [projcor_val] = projcorr_R(x, y, z)
%PDCORR_R computes the projection distance correlation of {x,y}|z.
%Underneath, it calls the R function pdcor from Matlab using a CSV file
%interface.
% Inputs:
%  x - [n x 1] input vector
%  y - [n x 1] input vector
%  z - [n x 1] input vector
% Outputs:
%  projcor_val - the partial distance correlation value
% Notes:
%  Check the startup.m script at the root of this repository to see how R
%  is called and ensure that you have the appropriate dependencies.  For
%  this script, the required R packages are:
%   'SAM', 'energy', 'glasso', 'glmnet', 'pracma', 'pgraph'
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

if(ismac || isunix)
    tmpDir = '/tmp';
else
    tmpDir = 'C:\Windows\Temp';
end

csvwrite(fullfile(tmpDir, 'x_projcor.csv'), x(:));
csvwrite(fullfile(tmpDir, 'y_projcor.csv'), y(:));
csvwrite(fullfile(tmpDir, 'z_projcor.csv'), z(:));

% call the R function
pgraphRdir = getPGraphPackagePath();
retCode = system(['R CMD BATCH ' pgraphRdir '/projcor_matlab.R' ]);
if(retCode==0)
    projcor_val = csvread(fullfile(tmpDir, 'projcor_matlab_output.csv'));
else
    warning('Calling R projcor function failed!');
    projcor_val = -999;
end

end