function [ tau ] = ktauhat( X, Y, correctionFlagOpt )
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

X = X(:);
Y = Y(:);

if(nargin<3)
    correctionFlagSelect = 4;       % default correction factor processing
else
    correctionFlagSelect = correctionFlagOpt;
end

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
    if(closeToZero(u,len))
        continuousRvIndicator = 0;
    else
        continuousRvIndicator = 1;
    end
    numOverlapPtsVec = countOverlaps(U, V, continuousRvIndicator);
    
    switch(correctionFlagSelect)
        case 1
            correctionFactor = correctionFactor1(numOverlapPtsVec);
        case 2
            correctionFactor = correctionFactor2(numOverlapPtsVec);
        case 3
            correctionFactor = correctionFactor3(numOverlapPtsVec);
        case 4
            correctionFactor = correctionFactor4(numOverlapPtsVec);
        case 5
            correctionFactor = correctionFactor5(numOverlapPtsVec);
        otherwise
            error('Unknown Correction Factor Option!');
    end
    
    t = max(u,v)-correctionFactor;
    tau = K/( sqrt(nchoosek(len,2)-t)*sqrt(nchoosek(len,2)-t) );
    
%     fprintf('K=%0.02f max(u,v)=%0.02f [%0.02f %0.02f %0.02f], correctionFactor=%0.02f nchoosek(M,2)=%d C-t=%d\n', ...
%         K, max(u,v), numOverlapPtsVec(1), numOverlapPtsVec(2), numOverlapPtsVec(3), ...
%         correctionFactor, nchoosek(len,2), nchoosek(len,2)-t);
else
    % case of either all continuous or all discrete data
%     fprintf('len=%d\n', len);
    tau = K/( sqrt(nchoosek(len,2)-u)*sqrt(nchoosek(len,2)-v) );
    
%     fprintf('K=%0.02f u=%0.02f v=%0.02f nchoosek(M,2)=%d C-u=%d C-v=%d\n', ...
%         K, u, v, nchoosek(len,2), nchoosek(len,2)-u, nchoosek(len,2)-v);
end

end

function [cf] = correctionFactor1(numOverlapPtsVec)
    if(min(numOverlapPtsVec)<2)
        cf = 0;
    else
        cf = nchoosek(floor(min(numOverlapPtsVec)),2)*length(numOverlapPtsVec);
    end

end

function [cf] = correctionFactor2(numOverlapPtsVec)
    if(min(numOverlapPtsVec)<2)
        cf = 0;
    else
        cf = nchoosek(floor(max(numOverlapPtsVec)),2)*length(numOverlapPtsVec);
    end
end

function [cf] = correctionFactor3(numOverlapPtsVec)

    cf = 0;
    for ii=1:length(numOverlapPtsVec)
        nop = floor(numOverlapPtsVec(ii));
        if(nop>=2)
            cf = cf + nchoosek(floor(numOverlapPtsVec(ii)),2);
        end
    end

end

function [cf] = correctionFactor4(numOverlapPtsVec)
    meanVal = floor(mean(numOverlapPtsVec));
    if(meanVal==0)
        cf = 0;
    else
        cf = nchoosek(meanVal,2)*length(numOverlapPtsVec);
    end
end

function [cf] = correctionFactor5(numOverlapPtsVec)
    cf = 0;
end

function [out] = closeToZero(in, len)
out = 1;
thresh = 0.02;      % if we are > 2% of length in terms of combinations;
lenFloor = floor(len*thresh);

if(lenFloor>=2)
    cmpVal = nchoosek(lenFloor,2);
else
    cmpVal = 0;
end

if(in>cmpVal)
    out = 0;
end

end

function [numOverlapPtsVec] = countOverlaps(U, V, continuousRvIndicator)
% this function is only called for hybrid data, attempts to correct for
% overestimation of the number of ties in hybrid data

M = length(U);

if(continuousRvIndicator==0)
    % U is the continuous RV
    continuousOutcomes = U;
    discreteOutcomes = V;
    % get the number of unique discrete outcomes
    uniqueDiscreteOutcomes = unique(V);
else
    % V is the continuous RV
    continuousOutcomes = V;
    discreteOutcomes = U;
    % get the number of unique discrete outcomes
    uniqueDiscreteOutcomes = unique(U);
end

% for each unique outcome .. count the overlapping elements.
numOverlapPtsVec = zeros(1,length(uniqueDiscreteOutcomes)-1);
for discreteOutcomesIdx=1:length(uniqueDiscreteOutcomes)-1
    % find the min/max values of the continuous values for this idx and the
    % next
    I = discreteOutcomes==uniqueDiscreteOutcomes(discreteOutcomesIdx);
    relevantContinuousOutcomes_curIdx = continuousOutcomes(I);
    I = discreteOutcomes==uniqueDiscreteOutcomes(discreteOutcomesIdx+1);
    relevantContinuousOutcomes_nextIdx = continuousOutcomes(I);
    
    % compute the number of points which are overlapping
    minCur = min(relevantContinuousOutcomes_curIdx);
    maxCur = max(relevantContinuousOutcomes_curIdx);
    
    numOverlapPoints = length(find(relevantContinuousOutcomes_nextIdx>=minCur & ...
                                   relevantContinuousOutcomes_nextIdx<=maxCur));
%     numOverlapPtsVec(discreteOutcomesIdx) = numOverlapPoints;
    numOverlapPtsVec(discreteOutcomesIdx) = numOverlapPoints/length(relevantContinuousOutcomes_nextIdx)*(M/length(uniqueDiscreteOutcomes));

end

end