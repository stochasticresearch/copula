function [W1, W2, W3, W4] = nonlinearDGP(M, W1_prev, W2_prev)
% LINEARDGP - Generates a nonlinear time-series described in the paper
% K. S. Narendra. Neural networks for intelligent control. 
%                 American control conference workshop, 1997
% Inputs:
%  M - the number of samples to generate
%  W1_prev - a vector with 1 element, describing W1(t-1)
%  W2_prev - a vector with 1 element, describing W2(t-1)
% Outputs
%  W1 - W1(t)
%  W2 - W2(t)
%  W3 - W3(t)
%  W4 - W4(t)
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

% input argument checking
if(length(W1_prev)<1)
    error('W1_prev needs to include W1(t-1)');
end
if(length(W2_prev)<1)
    error('W2_prev needs to include W2(t-1)');
end

W1 = zeros(M,1);
W2 = zeros(M,1);
W3 = zeros(M,1);
W4 = rand(M,1);

n1 = randn(M,1)*0.1;     % generate from N(0,0.1^2)
n2 = randn(M,1)*0.1;     % generate from N(0,0.1^2)
n3 = randn(M,1)*0.1;     % generate from N(0,0.1^2)

% seed the 1st element
W1(1) = (1 + W1_prev(1)/(1 + W1_prev(1)^2))*sin(W2_prev(1)) + n1(1);
W2(1) = W1_prev(1)*exp(-1*(W1_prev(1)^2 + W2_prev(1)^2)/8) + W2_prev(1)*cos(W2_prev(1)) + ...
    W4(1)^3/( (1+W4(1))^2 + 0.5*cos(W1_prev(1) + W2_prev(1)) ) + n2(1);
W3(1) = W1_prev(1)/(1+0.5*sin(W2_prev(1))) + W2_prev(1)/(1+0.5*sin(W1_prev(1))) + n3(1);

for ii=2:M
    W1(ii) = (1 + W1(ii-1)/(1 + W1(ii-1)^2))*sin(W2(ii-1)) + n1(ii);
    W2(ii) = W1(ii-1)*exp(-1*(W1(ii-1)^2 + W2(ii-1)^2)/8) + W2(ii-1)*cos(W2(ii-1)) + ...
        W4(ii-1)^3/( (1+W4(ii-1))^2 + 0.5*cos(W1(ii-1) + W2(ii-1)) ) + n2(ii);
    W3(ii) = W1(ii-1)/(1+0.5*sin(W2(ii-1))) + W2(ii-1)/(1+0.5*sin(W1(ii-1))) + n3(ii);

end

end