function val = mci(X,Y,Z,param)

% val = mci(X,Y,Z,param) computes the measure of conditional independence between X and Y given Z
% Input: X, Y and Z are column vectors of length n, X is [nx1], number of columns can vary for Y and Z
%             param is a struct with fields,
%                         kernel, for kernels such as gaus, laps
%                         kernelSize, a value such as 1 or a relation such as
%                               med for median intersample distance
%                               nnxx for mean xx-th nearest neighbor distance
%                         reg, regularization, a value such as 10^-3 or a relation such as
%                               n^(-1) for 1/(number of samples)
% Output: val is the value of mci
% Default: kernel is laps, kernelSize is med and ref is n^(-1)
% Author: Sohan Seth (sohan@cnel.ufl.edu)

n = size(X,1);
if size(Y,1) ~= n | size(Z,1) ~= n
    error('Length mismatch')
end

if nargin < 3
    error('Too few arguments')
end

if nargin == 3
    param.kernel = 'laps';
    param.kernelSize = 'med';
    param.reg = 'n^(-1)';
end

if ~isnumeric(param.reg)
    reg = eval(param.reg);
else
    reg = param.reg;
end

Kzz = grammat(Z,param);
Kuu = grammat([Y,Z],param);
Kxx = grammat(X,'idfs');

A1 = (Kzz' * Kzz) / n + reg*eye(n);
B1 = Kzz'*Kxx/n;
a1 = A1 \ B1;

% U ~ (Y,Z)
A2 = (Kuu' * Kuu) / n + reg*eye(n);
B2 = Kuu'*Kxx/n;
a2 = A2 \ B2;

val = sum(sum(((Kzz*a1 - Kuu*a2).^2)))/n^2;