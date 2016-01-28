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

function [ c ] = empcopdens_betak_3d_v2( X1,X2,X3,h,K )

K = K + 1;
c = zeros(K,K,K);

for ii=1:K
    for jj=1:K
        for kk=1:K
            uu = ii/K;
            vv = jj/K;
            ww = kk/K;
            
            K1 = betapdf(X1, uu/h + 1, (1-uu)/h + 1);
            K2 = betapdf(X2, vv/h + 1, (1-vv)/h + 1);
            K3 = betapdf(X3, ww/h + 1, (1-ww)/h + 1);

            val = sum(K1.*K2.*K3)/length(X1);
            c(ii,jj,kk) = val;
        end
    end
end

c = c(1:K-1,1:K-1,1:K-1);

end

