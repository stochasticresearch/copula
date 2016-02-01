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

function [ X ] = gen_samples_clg_2D( x1_domain, x1_multinomial_est, x2_mle_params, N )
%GEN_SAMPLES_CONDITIONAL_LINEAR_GAUSSIAN Generates samples of a conditional
%linear gaussian model
% Input:
%  x1_domain - the domain of the discrete random variable
%  x1_multinomial_est - the estimated frequency of each value in the
%                       x1_domain occuring
%  x2_mle_params - a [2 x length(x_domain)] array of the associated mu and
%                  sigma for each value of x1_domain
%  N - the number of samples to generate
% 
% Outputs:
%  X - the generated samples

X = zeros(N, 2);
x1 = mnrnd(1,x1_multinomial_est,N);

for ii=1:N
    % generate discrete RV
    discreteIdx = find(x1(ii,:)~=0);
    discreteRv = x1_domain(discreteIdx);
    
    % generate continuous RV based on generated discrete RV
    mu = x2_mle_params(1,discreteIdx);
    sigma = x2_mle_params(2,discreteIdx);
    continuousRv = normrnd(mu,sigma);
    
    X(ii,:) = [discreteRv continuousRv];
end


end

