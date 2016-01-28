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

clear;
clc;

D = 2;
K = 10;
x = zeros(K^D,D);

for ii=0:K^D-1
    
    value = ii;
    xIdx = D;
    while(value > 0)
        d = mod(value,K);
        
        x(ii+1,xIdx) = d;
        xIdx = xIdx - 1;
        
        value = floor(value/K);
    end
%     x(ii+1,:) = x(ii+1,:) + 1;
    x(ii+1,:) = x(ii+1,:)/(K-1);
end
