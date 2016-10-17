function [W1, W2, W3, W4, W5] = linearDGP(M, W1_prev, W4_prev, W5_prev)
% LINEARDGP - Generates a linear time-series described in the paper
% C. Zou and J. Feng. Granger causality vs. dynamic Bayesian network 
%    inference: a comparative study. BMC Bioinformatics, 10(1):122+, 2009.
% Inputs:
%  M - the number of samples to generate
%  W1_prev - a vector with 2 elements, describing W1(t-1) and W1(t-2)
%  W4_prev - a vector with 1 elements, describing W4(t-1)
%  W5_prev - a vector with 2 elements, describing W5(t-1) and W5(t-2)
% Outputs
%  W1 - W1(t)
%  W2 - W2(t)
%  W3 - W3(t)
%  W4 - W4(t)
%  W5 - W5(t)
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
if(length(W1_prev)<2)
    error('W1_prev needs to include W1(t-1), and W1(t-2)');
end
if(length(W4_prev)<1)
    error('W4_prev needs to include W4(t-1)');
end
if(length(W5_prev)<2)
    error('W5_prev needs to include W5(t-2)');
end

W1 = zeros(M,1);
W2 = zeros(M,1);
W3 = zeros(M,1);
W4 = zeros(M,1);
W5 = zeros(M,1);

n1 = randn(M,1)*0.6;     % generate from N(0,0.6^2)
n2 = randn(M,1)*0.5;     % generate from N(0,0.5^2)
n3 = randn(M,1)*0.3;     % generate from N(0,0.3^2)
n4 = randn(M,1)*0.3;     % generate from N(0,0.3^2)   
n5 = randn(M,1)*0.6;     % generate from N(0,0.6^2)

% seed the first 2 elements
W1(1) = 0.95*sqrt(2)*W1_prev(1) - 0.9025*W1_prev(2) + n1(1);
W2(1) = 0.5*W1_prev(2) + n2(1);
W3(1) = -0.4*W1_prev(2) + n3(1);
W4(1) = -0.5*W1_prev(1) + 0.25*sqrt(2)*(W4_prev(1) + W5_prev(1)) + n4(1);
W5(1) = -0.25*sqrt(2)*(W4_prev(1) - W5_prev(2)) + n5(1);

W2(2) = 0.5*W1_prev(1) + n2(2);
W1(2) = 0.95*sqrt(2)*W1(1) - 0.9025*W1_prev(1) + n1(2);
W3(2) = -0.4*W1_prev(1) + n3(2);
W4(2) = -0.5*W1(1) + 0.25*sqrt(2)*(W4(1) + W5(1)) + n4(2);
W5(2) = -0.25*sqrt(2)*(W4(1) - W5_prev(1)) + n5(2);

for ii=3:M
    W1(ii) = 0.95*sqrt(2)*W1(ii-1) - 0.9025*W1(ii-2) + n1(ii);
    W2(ii) = 0.5*W1(ii-2) + n2(ii);
    W3(ii) = -0.4*W1(ii-2) + n3(ii);
    W4(ii) = -0.5*W1(ii-1) + 0.25*sqrt(2)*(W4(ii-1) + W5(ii-1)) + n4(ii);
    W5(ii) = -0.25*sqrt(2)*(W4(ii-1) - W5(ii-2)) + n5(ii);
end

end