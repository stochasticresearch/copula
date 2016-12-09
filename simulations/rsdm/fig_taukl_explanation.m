%**********************************************************************
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
%* along with this program.  If not, see <http://www.gnu.org/licenses/>
%* 
%**********************************************************************

% Generate the figure which explains why we need additional corrections for
% hybrid data for tau_CJ

sz =80;
x = [1 2 3 4 5 6 7 8 9 10];
y = [1 1 1 1 1 2 2 2 2 2];
scatter(x,y,sz,'MarkerEdgeColor',[0 .5 .5],...
              'MarkerFaceColor',[0 .7 .7],...
              'LineWidth',1.5)
axis([0 11 0.85 2.15])
grid on;
xlabel('x', 'FontSize', 20);
ylabel('y', 'FontSize', 20);