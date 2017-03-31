function [model] = copmodelsel_td(X)
%COPMODELSEL_TD - performs copula model selection, based on tail-dependence
%and magnitude of correlation.  Only bivariate data is supported currently!
%focus is on speed, at the expense of accuracy.
%
% Inputs:
%  X - the input data matrix, of dimensions [N x 2]
% Outputs:
%  model - a string indicating the copula type, either:
%          Gaussian, Gumbel, Clayton
%
%**************************************************************************
%* 
%* Copyright (C) 2017  Kiran Karra <kiran.karra@gmail.com>
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
%* along with this program.  If not, see <http://www.gnu.org/licenses/>
%* 
%**************************************************************************

% compute srho, if it is negative, we can only model w/ Gaussian (we don't
% currently support rotated copulas, so Gaussian is our only choice for
% negative correlation modeling)
tauVal = taukl(X(:,1),X(:,2));
if(tauVal<0)
    model = 'Gaussian';
else
    % if there is significant lower or upper tail dependence, go w/ the
    % appropriate copula, else, use the Gaussian copula
    lt_dep = ltdepest(X);
    ut_dep = utdepest(X);
    lt_ut_ratio = lt_dep/ut_dep;
    ut_lt_ratio = 1/lt_ut_ratio;
    
    if(lt_ut_ratio > 1.25)
        % if we see strong tail lower-tail dependence, this is indicitave
        % of Clayton
        model = 'Clayton';
    elseif(ut_lt_ratio > 1.25)
        % if we see strong tail dependence on upper and lower tails, but
        % less dependence in the middle sections
        model = 'Gumbel';
    else
        model = 'Gaussian';
    end
end
end