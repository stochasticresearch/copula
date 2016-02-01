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

function [ c, a_vec, b_vec ] = empcopdens_betak_v2( u,v,h,K )

K = K + 1;
c = zeros(K,K);
a_vec = zeros(K,K);
b_vec = zeros(K,K);
for ii=1:K
    for jj=1:K
        a = ii/K;
        b = jj/K;
        
        K1 = betapdf(u,a/h + 1,(1-a)/h + 1);
        K2 = betapdf(v,b/h + 1,(1-b)/h + 1);
        
        c(ii,jj) = sum(K1.*K2)/length(u);
        a_vec(ii,jj) = a;
        b_vec(ii,jj) = b;
    end
end

c = c(1:K-1,1:K-1);

end

