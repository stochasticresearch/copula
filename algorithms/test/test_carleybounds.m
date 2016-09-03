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
%* 
%**************************************************************************

% Test the carley bounds calcuation.  The answers are compared to the 
% examples given in "A Primer on Copulas for Count Data" by 
% Genest and Neslehova

%% Replicate Example 3 and Figure 3 of above referenced paper
clear;
clc;

h_XY = zeros(2,2);
pVector = 0:0.01:1;
rVector = 0:0.01:1;

tau_C_minusMat = zeros(length(pVector),length(rVector));
tau_C_plusMat = zeros(length(pVector),length(rVector));

for pIdx=1:length(pVector)
    for rIdx=1:length(rVector)
        p = pVector(pIdx);
        q = p;
        r = rVector(rIdx);
        
        lowerBound = max(0, p+q-1);
        upperBound = min(p,q);
        
        if(r >= lowerBound && r <= upperBound)
            % generate h_XY
            h_XY(1,1) = r;					% X=0, Y=0
            h_XY(2,1) = q-r;				% X=1, Y=0
            h_XY(1,2) = p-r;				% X=0, Y=1
            h_XY(2,2) = 1-p-q+r;        	% X=1, Y=1
        
            % compute the carley bound and save result
            [ tau_C_plus, tau_C_minus ] = carleybounds( h_XY ); 
        else
            tau_C_plus = NaN;
            tau_C_minus = NaN;
        end
        
        tau_C_minusMat(pIdx,rIdx) = tau_C_minus;
        tau_C_plusMat(pIdx,rIdx) = tau_C_plus;
    end
end

% notice the conjugate operator on tau_C_plusMat and tau_C_minusMat, this
% is important because of the way surf works.  From Matlab's documentation:
% "vertices of the surface faces are (X(j), Y(i), Z(i,j)) triples"
% NOTE: you will need to rotate the 3d plots to the same view angle as
% Figure 3 in the above mentioned paper to notice the plots are the same.

subplot(1,2,1);
surf(pVector,rVector,tau_C_minusMat');
grid on;
xlabel('p'); ylabel('r'); zlabel('\tau_{-1}');

subplot(1,2,2);
surf(pVector,rVector,tau_C_plusMat');
grid on;
xlabel('p'); ylabel('r'); zlabel('\tau_{+1}');