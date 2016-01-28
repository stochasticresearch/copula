%******************************************************************************
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

function [ e ] = hyperFunctionError( f1, f2 )
%COMPUTEERRORBETWEENSURFACES Computes the error between two hyperfunctions
% as follows: sqrt( integral( (f1-f2).^2 ) )
%where the integral is computed over the entire domain of the functions.
%The two surfaces are expected to have the same domain.  I envision the 
% usage of this function to be a sort of a test of how good of an estimator 
% you have, where surf1 can be a surface (of any dimension) of the actual, 
% and surf2 can be an estimate.
% Inputs:
%  f1 - hyperfunction (of any dimension) one
%  f2 - hyperfunction (of the same dimension as f1) two
% Outputs:
%  e - the error

e = sqrt( nansum( (f1(:)-f2(:)).^2 ) ) / length(f1(:));    % normalize

end