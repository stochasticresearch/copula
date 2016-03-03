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

function [ X, Y, Z_internal_X ] = gen_conditionaldata( num_samps, model, n, ...
                                        distX, distY, distZ, ...
                                        fpre_XY, fcmb_XY, fpost_XY, ...
                                        fpre_Z, fcmb_Z, fpost_Z )
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
%  distX - a probability distribution object representing the distribution
%          of X
%  distY - a probability distribution object representing the distribution
%          of Y
%  distZ - a cell array of size [n x 1] of probability distribution
%          objects representing the distribution of Z_i
%  fpre_XY - a cell array of dimension [1 x 2] with function handles
%            defining the pre-operator applied to X.  If the cell array
%            element is empty, then no operator is applied and X is kept as
%            is.  fpre_XY{1} = f(X);  fpre_XY = f(Y);
%  fcmb_XY - a cell array of dimension [1 x 2] with function handles.
%            fcmb_XY{1} = f(fpre_XY{1},Z_cmb); 
%            fcmb_XY{2} = f(fpre_XY{2},Z_cmb);
%  fpost_XY - a cell array of dimension [1 x 2] with function handles
%            fpost_XY{1} = f(fcmb_XY{1})
%            fpost_XY{2} = f(fcmb_XY{2})
%  fpre_Z - a cell array of dimension [n x 1] where the i1^th element is a
%           function handle which defines the transform that will happen to
%           Z_i as applied to X, and the i2^th element is applied to Y.
%  fcmb_Z - a cell array of dimension [1 x 2] which contains the
%           function handles for a function which accepts a vector of
%           column vectors each representing Z_i and has a weight matrix on
%           how they should be combined.  0 weighted items are in effect
%           "independent" of the rest of the data, but they could be part
%           of the data set.  As before:
%            fcmb_Z{1} applies to X
%            fcmb_Z{2} applies to Y
%  fpost_Z - a cell array of function handles that are applied to the 
%            combined Z. As before:
%            fpost_Z{1} applies to X
%            fpost_Z{2} applies to Y
%            
%
% TODO: describe how the function handles can be used to create nice unique
% and flexible distributional data that is CI or nCI.  Explain m and n from
% diagrams w/ the weights for fcmb_Z

% generate Z_1 .. Z_n
Z_internal_X = zeros(num_samps, n);
Z_internal_Y = zeros(num_samps, n);

Z = zeros(num_samps, n);
for nn=1:n
    % generate data according to desired distribution
    Zi = random(distZ{nn}, num_samps, 1);
    Z(:,nn) = Zi;
    
    % apply initial operators to Z_i
    fZ_i_X = fpre_Z{nn, 1}; Zi_X = fZ_i_X(Zi);
    fZ_i_Y = fpre_Z{nn, 2}; Zi_Y = fZ_i_Y(Zi);
    
    % store
    Z_internal_X(:,nn) = Zi_X;
    Z_internal_Y(:,nn) = Zi_Y;
end

% Combine Z
Z_cmb_X = fcmb_Z{1}(Z_internal_X);
Z_cmb_Y = fcmb_Z{2}(Z_internal_Y);

% Apply the post operator for Z
Z_post_X = fpost_Z{1}(Z_cmb_X);
Z_post_Y = fpost_Z{2}(Z_cmb_Y);

% Generate X and Y
X = random(distX, num_samps, 1);
Y = random(distY, num_samps, 1);

% apply the pre operators
X = fpre_XY{1}(X);
Y = fpre_XY{2}(Y);

% apply the combine operator to X and Y
X_cmb = fcmb_XY{1}([X Z_post_X]);
Y_cmb = fcmb_XY{2}([Y Z_post_Y]);

% apply post operators to X and Y
X = fpost_XY{1}(X_cmb);
Y = fpost_XY{2}(Y_cmb);

% induce dependency between X and Y if desired
% TODO: the type of dependency should be configurable
if(strcmpi(model, 'nCI'))
    ff = normrnd(0,1,num_samps,1) * 0.5;
    X = X + ff; Y = X + ff;
end

end

