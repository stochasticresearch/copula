function [metric] = mv2_rsdm(x, y)
%MV2_RSDM - Multivariate 2 - Rank Statistics Dependency metric, 
%computes statistical dependency using rank statistics.  See associated paper.
% Inputs:
%  x - the x vector -- dimension [M x d1]
%  y - the y vector -- dimension [M x d2]
% Outputs:
%  metric - the calculated dependency metric between x and y
%
% See paper: "Measuring Association and Dependence Between Random Vectors"
% for theoretical basis behind this multivariate measure of concordance
%  http://www.sciencedirect.com/science/article/pii/S0047259X13001802
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

% compute the multivariate Kendall function for x, y, and [x y]
intAdA = mv_kendall(x);
intBdB = mv_kendall(y);
intCdC = mv_kendall([x y]);
metric = (intCdC - intAdA * intBdB) / sqrt( intAdA * (1 - intAdA) * intBdB * (1 - intBdB) );

end

% TODO: vectorize the function below :)
function [res] = mv_kendall(x)
% we keep the variables the same 
n = size(x,1);
d = size(x,2);

res = 0;
for ii=1:n
    for jj=1:n
        if(ii~=jj)
            v = x(jj,:)<=x(ii,:);
            if(sum(v)==d)
                res = res + 1;
            end
        end
    end
end

res = res/(n*n);

end