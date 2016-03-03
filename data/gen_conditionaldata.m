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
%**************************************************************************

function [ X, Y, Z ] = gen_conditionaldata( num_samps, model, m, n, ...
    depType, combineOperator, distX, distY, distZ )
%GEN_CONDITIONALDATA Generates either conditionally independent or not
%conditional independent data based on input parameters.  
% Conditionally independent data is generated with the following model:
%
%  +-------------------------------------------+
%  |                                           v
%+---------+     +-------+     +-------+     +---+
%|  Z_{1}  | --> |   X   | <-- | Z_{m} | --> | Y |
%+---------+     +-------+     +-------+     +---+
%                  ^                           ^
%                  |                           |
%                  |                           |
%+---------+     +-------+                     |
%| Z_{m+1} |     | Z_{2} | --------------------+
%+---------+     +-------+
%    .
%    .
%    .
%+---------+
%|  Z_{n}  |
%+---------+
% In other words, Z_{1} ... Z_{m| are parents of X and Y.  Z_{m+1} .. Z_{n}
% are additional nodes in the network but do not influence X and Y.
%
% Conditionally NOT independent data is generated with the following model:
%  +-------------------------------------------+
%  |                                           v
%+---------+     +-------+     +-------+     +---+
%|    E    | --> |       | <-- | Z_{m} | --> |   |
%+---------+     |   X   |     +-------+     | Y |
%+---------+     |       |     +-------+     |   |
%| Z_{m+1} |     |       | <-- | Z_{1} | --> |   |
%+---------+     +-------+     +-------+     +---+
%    .             ^                           ^
%    .             |                           |
%    .             |                           |
%+---------+     +-------+                     |
%|  Z_{n}  |     | Z_{2} | --------------------+
%+---------+     +-------+
% In other words, Z_{1} ... Z_{m| are parents of X and Y.  Z_{m+1} .. Z_{n}
% are additional nodes in the network but do not influence X and Y.  X and
% Y are also influenced by a common "error" term E.
%
% Inputs:
%  num_samps - the number of samples to generate
%  model - should be either 'CI' meaning the conditionally independent
%          model above, or 'nCI', which is the NOT conditionally
%          independent model (the 2nd diagram)
%  m - the number of nodes in Z which influence X and Y
%  n - the total number of nodes in Z
%  depType - a cell array of dimension [n x 2], where the ij^th element is
%            a function handle to the operation that should be applied to
%            induce the desired dependency.  The index i-1 represents 
%            Z_i -> X and i-2 represents Z_i -> Y
%  combineOperator - a cell array of dimension [2 x 1], where the first
%                    element is a function handle to combine the
%                    dependencies from Z_1 ... Z_n into X, and the second
%                    element is a function handle to combine the
%                    dependencies from Z_1 ... Z_n into Y.
%   distX - a probability distribution object representing the distribution
%           of X
%   distY - a probability distribution object representing the distribution
%           of Y
%   distZ - a cell array of size [n x 1] of probability distribution
%           objects representing the distribution of Z_i


end

