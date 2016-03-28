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

function [ C_u ] = empcopula_val( C, U )
%EMPCOPULA_VAL Returns the value of C (the copula function) and 
% c (the copula density) at the specified U (vector).
% Inputs:
%  C - the empirical copula (or copula density) function, defined over an 
%      evenly spaced grid of points
%  U - A vector [U_1, U_2, ... U_d], where d is the dimensionality
%      of the copula (distribution and density functions)
%
% Outputs:
%  C_u - the value of C(U), if C is empty, C_u will be 0
%
% TODO:
%  [ ] - U should have the ability to be a matrix as well


d = length(size(C));

C_u = 0;

% calculate which grid points land in the U values of interest
linearIdx = 0;
for ii=1:d
    if(~isempty(C))
        K_ii = size(C,ii);
    elseif(~isempty(c))
        K_ii = size(c,ii);
    else
        error('Need to specify atleast one of the two: C or c');
    end
    
    idx = ceil(U(ii)*K_ii); % WARNING -- used to be :: idx = round(U(ii)*K_ii);
                            % but w/ testing, I think ceil works better,
                            % this is due to the way linspace generates
                            % data
    if(idx>K_ii)
        idx = K_ii;
    elseif(idx<1)
        idx = 1;
    end
    
    if(ii==1)
        linearIdx = idx;
    else
        linearIdx = linearIdx + (idx-1)*(K_ii.^(ii-1));
    end
end

% fprintf('u1=%f u2=%f linearIdx=%d\n', U(1), U(2), linearIdx);

% extract value of C_u from those
if(~isempty(C))
    C_u = C(linearIdx);
end

end

