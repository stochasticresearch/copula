function [ X_out ] = continueRv( X_in )
%CONTINUERV Continues a random variable according to Michel and Denuit
%(2005), and Neslehova (2007).
% Inputs:
%  X_in - Input matrix of dimension M x D, where D is the dimensionality of
%         the data, and M is the number of samples.
%
% Outputs:
%  X_out - the continued samples of the random variables
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

X_out = X_in + (rand(size(X_in))-1);

end

