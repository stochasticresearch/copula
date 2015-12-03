function [ U ] = empcopularnd( c, M )
%EMPCOPULARND Generates samples from an empirical copula
%
% Inputs:
%  c - The copula density, provided as a [K x K x ... x K] matrix, where
%      the dimensions are [u_1,u_2, ..., u_D]
%  M - The number of samples to generate from this copula
%
% Outputs
%  U - an M x D vector of samples from the copula, given by the integral
%      of the provided copula density.

    K = size(c,1);
    uu = linspace(0,1,K);
    D = length(size(c));

    c_d = cell(1,D);    % c_d{i} stores c_i(u_1,...,u_i) = integral(c,du_{i+1} ... du_{D})
                        % For more information, refer to: Analysis and
                        % Generation of Random Vectors with Copulas, by Johann
                        % Strelen and Feras Nassaj
    c_d{D} = c;
    % integrate out each dimension successively
    for ii=D-1:-1:1
        c_d{ii} = squeeze(sum(c_d{ii+1},ii+1));
    end

    U = rand(M,D);

    for ii=1:M
        u_j = U(ii,1);
        for jj=2:D
            idx = findClosest(uu,u_j);
        end
    end

end

function [idx] = findClosest(vec, val)
    tmp = abs(val-vec);
    [~,idx] = min(tmp);
end