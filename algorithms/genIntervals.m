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

function [intervalPrune, numSubsets] = genIntervals(f, xi)

% determine inflection points
[peaks, valleys] = findinflections(0,f);
allPts = sort([peaks(:); valleys(:)]');
numSubsets = length(allPts)+1;
allPtsAugment = [0 allPts length(xi)];

% make a matrix of intervals over which to calculate the MTE parameters
intervals = zeros(numSubsets,2);
intervals(1,1) = allPtsAugment(1)+1; intervals(1,2) = allPtsAugment(2);    
minIntervalSize = 4;
for ii=2:numSubsets
    % if any intervals are less than a configurable amount of points,
    % extend the intervals so that curve fitting can work
    prevIntervalSize = intervals(ii-1,2)-intervals(ii-1,1);
    if(prevIntervalSize<minIntervalSize)
        endPoint = intervals(ii-1,2) + (minIntervalSize-prevIntervalSize);
        intervals(ii-1,2) = endPoint;
        % fast-forward ii to the end-point to avoid overlap points
        tmp = find(allPtsAugment<endPoint);
        ii = tmp(end)+1;
    end

    intervals(ii,1) = intervals(ii-1,2)+1;
    intervals(ii,2) = allPtsAugment(ii+1);
end
% if the last interval is less than the minimum interval size, merge it
% into the second to last interval
lastIntervalSize = intervals(end,2)-intervals(end,1);
if(lastIntervalSize<minIntervalSize)
    intervals(end-1,2) = length(xi);
    intervals(end,:) = [];
    numSubsets = numSubsets - 1;
end

% find any other intervals that have 0's and delete them
intervalPrune = [];
for ii=1:size(intervals,1)
    if(intervals(ii,1)~=0 && intervals(ii,2)~=0)
        intervalPrune = [intervalPrune; [intervals(ii,1) intervals(ii,2)]];
    end
end

end