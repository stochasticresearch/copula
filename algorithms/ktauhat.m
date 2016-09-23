function [ tau ] = ktauhat( X, Y )
%ktaucj - computes a rescaled version of Kendall's tau that preserves
%         the definition of Kendall's tau, but assures that in the 
%         scenario of perfect concordance or discordance for discrete
%         or hybrid datatypes, taucj achieves +/- 1 respectively
% Inputs:
%  X - first variable input.
%  Y - second variable input.
% Outputs:
%  tau - the rescaled version of Kendall's tau
%  
% NOTE: See Theorem 12 of the paper "On Rank Correlation Measures for
%       Non-Continuous Random Variables" - Neslehova (2007) for more info.
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

% TODO: some error checking that X and Y are the same length, NaN, Inf ...
% etc...


% TOOD: compare the 2 ways to generate ranks ... wonder if that affects the
% results?
% rank the data, do not account for ties
data_sorted = sort(X);
[~, U] = ismember(X,data_sorted);

data_sorted = sort(Y);
[~, V] = ismember(Y,data_sorted);

% U = tiedrank(X);
% V = tiedrank(Y);

% compute the numerator the tau_hat
K = 0;
len = length(X);
for k = 1:len-1
    K = K + sum( sign(U(k)-U(k+1:len)) .* sign(V(k)-V(k+1:len)) );
end

% compute the denominator ... compute the # of unique values of U and V and
% how many times each of those unique values occur
uniqueU = unique(U);
uniqueV = unique(V);

uniqueUCounts = zeros(1,length(uniqueU));
uniqueVCounts = zeros(1,length(uniqueV));

% TODO: we can combine the loops below after verification of functionality
for ii=1:length(uniqueU)
    uniqueUCounts(ii) = sum(U==uniqueU(ii));
end

for ii=1:length(uniqueV)
    uniqueVCounts(ii) = sum(V==uniqueV(ii));
end

u = 0;
k = 2;
for ii=1:length(uniqueU)
    n = uniqueUCounts(ii);
    if(k<=n)
        addVal = nchoosek(n,k);
    else
        addVal = 0;
    end
    u = u + addVal;
end
v = 0;
for ii=1:length(uniqueV)
    n = uniqueVCounts(ii);
    if(k<=n)
        addVal = nchoosek(n,k);
    else
        addVal = 0;
    end
    v = v + addVal;
end

if( (closeToZero(u, len) && v>0) || (u>0 && closeToZero(v, len)) )
    % special case of hybrid data
    t = max(u,v);
    tau = K/( sqrt(nchoosek(len,2)-t)*sqrt(nchoosek(len,2)-t) );
else
    % case of either all continuous or all discrete data
    tau = K/( sqrt(nchoosek(len,2)-u)*sqrt(nchoosek(len,2)-v) );
end

end

function [out] = closeToZero(in, len)
out = 1;
thresh = 0.02;      % if we are > 2% of length in terms of combinations;
lenFloor = floor(len*thresh);
if(in>nchoosek(lenFloor,2))
    out = 0;
end

end