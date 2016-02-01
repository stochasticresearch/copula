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

clear;
clc;

x = [10 10];
y = zeros(prod(x), length(x));
for ii=1:prod(x)
    remVal = ii-1;
    divVal = prod(x(1:end-1));
    for jj=1:length(x)
        y(ii,jj) = floor( remVal/divVal );
%         fprintf('ii=%d remVal=%d divVal=%d y(%d,%d)=%d \n', ...
%             ii, remVal, divVal, ii, jj, y(ii,jj));
        
        remVal = remVal-y(ii,jj)*divVal;
        divVal = divVal / x(1);
    end
%     fprintf('*****************************\n');
end
y = y+1;
y = y / x(end);
% fliplr(y)