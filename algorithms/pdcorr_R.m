function [pdcor_val] = pdcorr_R(x, y, z)
%PDCORR_R computes the partial distance correlation of {x,y}|z.
%Underneath, it calls the R function pdcor from Matlab using a CSV file
%interface.
% Inputs:
%  x - [n x 1] input vector
%  y - [n x 1] input vector
%  z - [n x 1] input vector
% Outputs:
%  pdcor_val - the partial distance correlation value
% Notes:
%  Check the startup.m script at the root of this repository to see how R
%  is called and ensure that you have the appropriate dependencies.  For
%  this script, the only R package required is the "energy" package, and
%  the "pracma" package
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

csvwrite(fullfile(tmpDir, 'x.csv'), x(:));
csvwrite(fullfile(tmpDir, 'y.csv'), y(:));
csvwrite(fullfile(tmpDir, 'z.csv'), z(:));

% call the R function
energyRdir = getEnergyPackagePath();
retCode = system(['R CMD BATCH ' energyRdir '/pdcor_matlab.R' ]);
if(retCode==0)
    pdcor_val = csvread(fullfile(tmpDir, 'pdcor_matlab_output.csv'));
else
    warning('Calling R pdcor function failed!');
    pdcor_val = -999;
end

% clean up temporary files
delete('pdcor_matlab.R', 'pdcor_matlab.Rout');

end