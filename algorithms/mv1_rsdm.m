function [metric] = mv1_rsdm(x, y, varargin)
%MV1_RSDM - Multivariate 1 - Rank Statistics Dependency metric, 
%computes statistical dependency using rank statistics.  See associated paper.
% Inputs:
%  x - the x vector -- dimension [M x d1]
%  y - the y vector -- dimension [M x d2]
%  varargin{1} - minscanincr - the minimum scanning increment.  Large
%                              values will filter out high frequency
%                              dependencies, small values decrease the
%                              statistical power of the dependency metric
%  varargin{2} - diffthresh  - the threshold at which a change in
%                              concordance amount is detected.  Larger
%                              values are more robust to noise, but tend to
%                              miss high frequency changes.
%  varargin{3} - alpha       - the value used to determine significance
%                              level of a box's concordance level
% Outputs:
%  metric - the calculated dependency metric between x and y
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

nVarargin = length(varargin);
switch nVarargin
    case 0
    case 1
        minscanincr = varargin{1};
    case 2
        minscanincr = varargin{1};
        diffthresh = varargin{2};
    otherwise
        % the >=3 case
        minscanincr = varargin{1};
        diffthresh = varargin{2};
        alpha = varargin{3};
end

d1 = size(x,2);
d2 = size(y,2);
if(size(x,1)~=size(y,1))
    error('The number of examples provided for each multivariate distribution must be the same!');
end

% compute pairwise RSDM
combos = combvec(1:d1, 1:d2)';
metric = 0;
for ii=1:length(combos)
    combo = combos(ii,:);
    xx = x(:,combo(1));
    yy = y(:,combo(2));
    
    switch nVarargin
        case 0
        case 1
            metricTmp = rsdm(xx,yy,minscanincr);
        case 2
            metricTmp = rsdm(xx,yy,minscanincr,diffthresh);
        otherwise
            % the >=3 case
            metricTmp = rsdm(xx,yy,minscanincr,diffthresh, alpha);
    end
    metric = metric + metricTmp;
end

% scale the RSDM
metric = metric/size(combos,1);

end