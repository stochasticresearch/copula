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

% script which tests the continueRv function

M = 1000;
D = 2;
X_in = unidrnd(5,M,D)-1;
X_out = continueRv(X_in);

subplot(1,2,1);
scatter(X_in(:,1), X_in(:,2));
grid on;
title('Original Discrete RVs')
subplot(1,2,2);
scatter(X_out(:,1), X_out(:,2));
grid on;
title('Continued Discrete RVs')