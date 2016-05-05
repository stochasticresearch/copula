function [ y, edges ] = discretizeRv( x, N, varargin )
%DISCRETIZERV Discretizes a random variable using a specified method
% Inputs:
%  x - a vector of samples from a continuous random variable to be
%  discretized
%  N - the number of discrete intervals
%
% Optional Inputs:
%  varargin{1} - if 0, then compute pdf directly, 
%                else, compute log of pdf at specified value
% Outputs:
%  y - the discretized binned data
%  empInfo - the empirical information (i.e. pdf, cdf, domain)
%
% Acknowledgements:
% This file was inspired by the following research paper:
%  METHODS FOR DISCRETIZING CONTINUOUS VARIABLES WITHIN THE FRAMEWORK 
%  OF BAYESIAN NETWORKS
%  Mihaela-Daciana Craciun, Violeta Chis¸ and Cristina Bala 
%  ACTA UNIVERSITATIS APULENSIS - Special Issue
%
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

[~,edges,y] = histcounts(x, N, 'Normalization', 'pdf');
y = edges(y);