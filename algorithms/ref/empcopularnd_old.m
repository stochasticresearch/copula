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

function [ U_gen ] = empcopularnd_old( c, N )
%EMPCOPULARND Generates samples of the empirical copula of an underlying
%d-dimensional joint distribution function
% Inputs:
%   c - The conditional densities of the copula function.  This should be a
%       cell array with the format:
%                  c{1} = []
%                  c{2} = C(u_2|u_1)
%                  c{3} = C(u_3|u_1,u_2)
%                           .
%                           .
%                           .
%                  c{D} = C(u_D|u_1,u_2,...,u_{d-1})
%       Resultingly, c{D} will be the "copula density."
%   N - the number of copula samples to generate.  This will be internally
%       rounded to a multiple of K.
% Outputs:
%   U_gen - an [N x D] matrix of generated copula samples representing the 
%           empirical copula of provided joint distribution x.
% Credits: The algorithm and prototype code from which this function
%          was created from the algorithms outlined in the paper:
%          "Analysis and Generation of Random Vectors with Copulas"
%          and the code on Johann Strelen's website:
%          http://web.cs.uni-bonn.de/IV/strelen/Algorithmen/

K = size(c{end},1);
D = length(c);

% rescale c
f = c;
for d=2:length(c)
    f{d} = c{d}*(K^d);
end

n_gen = N;
U_gen = zeros(D,n_gen);         % empirical copula samples
U_gen(1,:) = rand(1,n_gen);

power_K = realpow(K,D-1);

delta = 1/K;

for i=1:n_gen
    u_gen_ = rand(1,1);
    u_power_K = u_gen_ * power_K;

    summ = 0;
    j2 = 0;
    u1_gen_up = ceil(U_gen(1,i)*K);   % delta = 1/K
    while (summ < u_power_K) 
        j2 = j2+1;
        summ_old = summ;
        summ = summ + f{2}(u1_gen_up, j2);
    end;
    u2_down = j2-1;
    u2_gen_up = j2;
    u_gen_rest = u_gen_ - summ_old/power_K;
    u2_gen_rest = u_gen_rest*realpow(K,D-2) / f{2}(u1_gen_up, u2_gen_up);
    U_gen(2,i) = u2_down*delta + u2_gen_rest;
    
    un_gen_up = [u1_gen_up u2_gen_up zeros(1,D-2)];
    for d=3:D
        u_gen_ = rand(1,1);
        summ = 0;
        jn = 0;
        
        % calculate linear index
        access_vec = un_gen_up(1:d-1);
        linear_idx = access_vec(1);
        for kk=2:d-1
            linear_idx = linear_idx + (access_vec(kk)-1)*(K.^(kk-1));
        end
        
        u_gen_fdm1u1udm1 = u_gen_*f{d-1}(linear_idx);      % u_gen f_{d-1} u_1 ... u_{d-1}
        
        while (summ < u_gen_fdm1u1udm1) 
            jn = jn+1;
            summ_old = summ;
            
            access_vec = [un_gen_up(1:d-1) jn];
            linear_idx = access_vec(1);
            for kk=2:d
                linear_idx = linear_idx + (access_vec(kk)-1)*(K.^(kk-1));
            end
            summ = summ + f{d}(linear_idx);
        end
        ud_down = jn-1;
        ud_gen_up = jn;
        u_gen_rest = delta*(u_gen_fdm1u1udm1 - summ_old);
        
        access_vec = [un_gen_up(1:d-1) jn];
        linear_idx = access_vec(1);
        for kk=2:d
            linear_idx = linear_idx + (access_vec(kk)-1)*(K.^(kk-1));
        end
        ud_gen_rest = u_gen_rest / f{d}(linear_idx);
        
        U_gen(d,i) = ud_down*delta + ud_gen_rest;
        un_gen_up(d) = ud_gen_up;
        
    end
    
end

U_gen = U_gen';

end

