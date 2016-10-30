function [pdcor_val, pval] = pdcorr_R(x, y, z, R)
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

if(nargout>1)
    function_call = '/pdcor_pval_matlab.R';
else
    % we only need the pdcorr value here
    function_call = '/pdcor_matlab.R';
end

if(nargin<4)
    R = 199;    % the number of replicates for the permutation test to generate
                % a p-value
end

if(ismac || isunix)
    tmpDir = '/tmp';
else
    tmpDir = 'C:\Windows\Temp';
end

csvwrite(fullfile(tmpDir, 'x_pdcor.csv'), x(:));
csvwrite(fullfile(tmpDir, 'y_pdcor.csv'), y(:));
csvwrite(fullfile(tmpDir, 'z_pdcor.csv'), z(:));
csvwrite(fullfile(tmpDir, 'pdcorr_replicates.csv'), R(:));

% call the R function
energyRdir = getEnergyPackagePath();
retCode = system(['R CMD BATCH ' energyRdir function_call ]);
if(retCode==0)
    oo = csvread(fullfile(tmpDir, 'pdcor_matlab_output.csv'));
    pdcor_val = oo(1); 
    if(nargout>1)
        pval = oo(2);
    end
else
    warning('Calling R pdcor function failed!');
    pdcor_val = -999; pval = -999;
end

% clean up temporary files
if(nargout>1)
    delete('pdcor_pval_matlab.Rout');
else
    delete('pdcor_matlab.Rout');
end

end