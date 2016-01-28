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

function [ c ] = empcopdens_betak( u,v,h,K )

K = K + 1;
c = zeros(K,K);
for ii=1:K
    for jj=1:K
        a = ii/K;
        b = jj/K;
        c(ii,jj) = sum(betapdf(a,u/h,(1-u)/h).*betapdf(b,v/h,(1-v)/h))/length(u);
    end
end
c = c(1:K-1,1:K-1);

end

