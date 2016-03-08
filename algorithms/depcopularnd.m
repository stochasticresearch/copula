function [ U ] = depcopularnd( Z, D, family, depParam )
%DEPCOPULARND Aids in the generation of generating random numbers from
%graphical models with common parents.
% Inputs:
%  Z - an [M x 1] vector of the parent's common random variable.  For the
%  Gaussian family, this should be U_i(:,1), and for the 
%  D - the dimensionality of the data to generate
%  Clayton/Frank/Gumbel families, this should be X_i(:,1)
%  family - the desired family, can be either:
%           Gaussian
%           Clayton
%           Frank
%           Gumbel
%  depParam - for the Gaussian family, this must be a correlation matrix,
%             for Clayton/Frank/Gumbel, this should be the alpha value
% Outputs:
%  U - the copula random variables

M = size(Z,1);
family_LC = lower(family);
switch family_LC
    case 'gaussian'
        UU = norminv(Z);
        VV = mvnrnd(zeros(1,D-1),eye(D-1),M);
        U_uncorrelated = [UU VV];
        A = chol(depParam);
        U = normcdf(U_uncorrelated*A);
    case 'gumbel'
    case 'frank'
    case 'clayton'
        U = claytoncopularnd(M, D, depParam, Z);
    otherwise
        error('Unsupported Copula Family Type!');
end

end

