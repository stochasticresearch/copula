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
        idxVec = zeros(1,D); idxVec(1) = findClosest(uu,u_j);
        for jj=2:D
            C_d_num = getMarginalIntegral(c_d{jj},idxVec(1:jj-1),K);
            C_d_den = c_d{jj-1}(getLinearIdx(idxVec(1:jj-1),K));
            
            C_d = C_d_num./C_d_den;
            % perform numerical inverse
            u_j = uu(findClosest(C_d,U(ii,jj)));
            U(ii,jj) = u_j;
            
            idxVec(jj) = findClosest(uu,u_j);
        end
    end

end

function [linearIdx] = getLinearIdx(arrIdx,K)
    linearIdx = arrIdx(1);
    multiplyFactor = K;
    for ii=2:length(arrIdx)
        linearIdx = linearIdx + (arrIdx(ii)-1)*multiplyFactor;
        multiplyFactor = multiplyFactor*K;
    end
end

function [y] = getMarginalIntegral(c, idxVec, K)
    y = zeros(1,K);
    % extract the desired dimension
    for ii=1:K
        arrIdx = [idxVec ii];
        linearIdx = getLinearIdx(arrIdx,K);
        y(ii) = c(linearIdx);
    end
    y = cumsum(y);      % integrate
end

function [idx] = findClosest(vec, val)
    tmp = abs(val-vec);
    [~,idx] = min(tmp);
end