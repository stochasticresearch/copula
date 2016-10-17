function val = hd(X,Y,Z,param)

% val = hd(X,Y,Z,param) computes the Hellinger distance between conditional PDFs
% Input: X, Y and Z are column vectors of length n, X is [nx1], number of columns can vary for Y and Z
%             param is a struct with fields,
%                         kernel, for kernels such as gaus, laps
%                         kernelSize, a value such as 1 or a relation such as
%                               med for median intersample distance
%                               nnxx for mean xx-th nearest neighbor distance
% Output: val is the value of Hellinger distance
% Default: kernel is gauf and kernelSize is n^(-1/8.5)
% Author: Sohan Seth (sohan@cnel.ufl.edu)

n = size(X,1);
if size(Y,1) ~= n | size(Z,1) ~= n
    error('Length mismatch')
end

if nargin < 3
    error('Too few arguments')
end

if nargin == 3
    param.kernel = 'gauf';
    param.kernelSize = 'n^(-1/8.5)';
end

X = zscore(X); Y = zscore(Y); Z = zscore(Z);

weight.kernel = 'tris';
weight.kernelSize = 2;
weightVector = grammat([X,Y,Z],zeros(1,size([X,Y,Z],2)),weight) / 2^(size([X,Y,Z],2));

Kxz = mean(grammat([X,Z],param),2); Kxz(Kxz <= 0) = 0.1/n;
Kyz = mean(grammat([Y,Z],param),2); Kyz(Kyz <= 0) = 0.1/n;
Kz = mean(grammat(Z,param),2); Kz(Kz <= 0) = 0.1/n;
Kxyz = mean(grammat([X,Y,Z],param),2); Kxyz(Kxyz <= 0) = 0.1/n;

val = mean( (1 - sqrt((Kxz .* Kyz) ./ (Kz .* Kxyz))) .^2  .* weightVector);