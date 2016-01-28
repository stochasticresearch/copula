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

function [ U_gen, Z_sorted, U_emp, f ] = emp_copularnd_old( x, N, K )
% function [ U_gen, Z_sorted, U_emp, f, nU ] = emp_copularnd( x, N, K )
%EMP_COPULARND Generates samples of the empirical copula of an underlying
%d-dimensional joint distribution function
% Inputs:
%   x - the D-dimensional distribution function.  This should be in the
%       format of an [M x D] matrix, where M is the number of samples and D
%       is the dimensionality of the distribution function.
%   K - the number of subsets of [0,1].  Should be a minimum of 100.  The
%       generated samples become more accurate as K is increased, but the
%       tradeoff is time required to compute copula samples.
%   N - the number of copula samples to generate.  This will be internally
%       rounded to a multiple of K.
% Outputs:
%   u - an [N x D] matrix of copula samples representing the empirical
%       copula of provided joint distribution x.
%
% TODO:
%   [ ] - Extend to more than 4 dimensions
%
%
% Credits: The algorithm and prototype code from which this function
%          was created from the algorithms outlined in the paper:
%          "Analysis and Generation of Random Vectors with Copulas"
%          and the code on Johann Strelen's website:
%          http://web.cs.uni-bonn.de/IV/strelen/Algorithmen/

M = size(x,1);
D = size(x,2);

x = x';

n_by_K = floor(M/K);
n = n_by_K * K;         %sample size = floor(n_whole_sample/K)*K

% truncate the data to n
x = x(:,1:n);

[Z_sorted, nU, ~] = empVF_v3(n,D,x); 


delta = 1/K;
j=zeros(1,D);
f2=zeros(K,K);
summand = 1*realpow(K,D)/n;

switch D
    case 4
        f4 = sparse(K*K+K, K*K+K); % Access f4(i,j,k,l) (not sparse) as follows:
                                   % f4(i*K + j, k*K + l) in the sparse
                                   % matrix!
        f3 = sparse(K*K+K,K);      % Access f3(i,j,k) (not sparse) as follows:
                                   % f3(i*K + j, l) in the sparse
                                   % matrix!
        for i=1:n
            for d=1:D
                j(d) = ceil( (nU(d,i)-0.00001)/n_by_K );   % rounding errors omitted
            end;
            f4(j(1)*K+j(2),j(3)*K+j(4)) = f4(j(1)*K+j(2),j(3)*K+j(4)) + summand;
            f3(j(1)*K+j(2),j(3)) = f3(j(1)*K+j(2),j(3)) + summand;
            f2(j(1),j(2)) = f2(j(1),j(2)) + summand;
        end
    case 3
        f3 = sparse(K*K+K,K);      % Access f3(i,j,k) (not sparse) as follows:
                                   % f3(i*K + j, l) in the sparse
                                   % matrix!
        for i=1:n
            for d=1:D
                j(d) = ceil( (nU(d,i)-0.00001)/n_by_K );   % rounding errors omitted
            end;
            
            f3(j(1)*K+j(2),j(3)) = f3(j(1)*K+j(2),j(3)) + summand;
            f2(j(1),j(2)) = f2(j(1),j(2)) + summand;
        end
    case 2
        for i=1:n
            for d=1:D
                j(d) = ceil( (nU(d,i)-0.00001)/n_by_K );   % rounding errors omitted
            end;
            
            f2(j(1),j(2)) = f2(j(1),j(2)) + summand;
        end
end %switch

% testing outputs only
f = {}; f{1} = [];
switch D
    case 4
        f{4} = f4;
        f{3} = f3;
        f{2} = f2;
    case 3
        f{3} = f3;
        f{2} = f2;
    case 2
        f{2} = f2;
end

U_emp=nU/n;

n_gen = N;
U_gen = zeros(D,n_gen);         %copula
U_gen(1,:) = rand(1,n_gen);

power_K = realpow(K,D-1);

for i=1:n_gen
    u_gen_ = rand(1,1);
    u_power_K = u_gen_ * power_K;

    summ = 0;
    j2 = 0;
    u1_gen_up = ceil(U_gen(1,i)*K);   % delta = 1/K
    while (summ < u_power_K) 
        j2 = j2+1;
        summ_old = summ;
        summ = summ + f2(u1_gen_up, j2);
    end;
    u2_down = j2-1;
    u2_gen_up = j2;
    u_gen_rest = u_gen_ - summ_old/power_K;
    u2_gen_rest = u_gen_rest*realpow(K,D-2) / f2(u1_gen_up, u2_gen_up);
    U_gen(2,i) = u2_down*delta + u2_gen_rest;

    if D>2
    u_gen_ = rand(1,1);
    summ = 0;
    j3 = 0;
    u_gen_f2u1u2 = u_gen_*f2(u1_gen_up,u2_gen_up);

    while (summ < u_gen_f2u1u2) 
        j3 = j3+1;
        summ_old = summ;
        summ = summ + f3(u1_gen_up*K+u2_gen_up, j3);
    end;
    u3_down = j3-1;
    u3_gen_up = j3;
    u_gen_rest = delta*(u_gen_f2u1u2 - summ_old);
    u3_gen_rest = u_gen_rest / f3(u1_gen_up*K+u2_gen_up, u3_gen_up);
    U_gen(3,i) = u3_down*delta + u3_gen_rest;
    end; % D>2

    if D>3
    u_gen_ = rand(1,1);
    summ = 0;
    j4 = 0;
    u_gen_f3u1u2u3 = u_gen_*f3(u1_gen_up*K+u2_gen_up,u3_gen_up);

    while (summ < u_gen_f3u1u2u3) 
        j4 = j4+1;
        summ_old = summ;
        summ = summ + f4(u1_gen_up*K+u2_gen_up,u3_gen_up*K+ j4);
    end;
    u4_down = j4-1;
    u_gen_rest = delta*(u_gen_f3u1u2u3 - summ_old); 
    u4_gen_rest = u_gen_rest / f4(u1_gen_up*K+u2_gen_up,u3_gen_up*K+u4_down+1);
    U_gen(4,i) = u4_down*delta + u4_gen_rest;
    end; % D>3

end;

U_gen = U_gen';
U_emp = U_emp';

end

