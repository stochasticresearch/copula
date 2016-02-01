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

function [Z_sorted, nU, F_emp] = empVF_v3(n,D,Z)

%Empirical distribution:

step=[1/n:1/n:1]';  % [1/n, 2/n, ... ,1]

for d=1:D
   [Z_sorted(d,:), index_alt(d,:)] = sort(Z(d,:));
end 


nU=zeros(D,n);     % to speed up processing

%Find  empirical CDF  F_emp
F_emp = zeros(4,n);
for d=1:D
    F_emp(d,:) = step';
end

% Sampling values repeatedly have different values F i/n, (i+1)/n, ...
for i=1:n
    for d=1:D        
        nU(d,index_alt(d,i)) =  i; 
        % values n*U(d,index_alt) in order to omit rounding errors
    end                
end

end