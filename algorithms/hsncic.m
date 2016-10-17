function d = hsncic(X,Y,Z,param)

% val = hsncic(X,Y,Z,param) computes the Hilbert-Schmidt normalized conditional independence criterion between X and Y given Z
% Input: X, Y and Z are column vectors of length n, X is [nx1], number of columns can vary for Y and Z
%             param is a struct with fields,
%                         kernel, for kernels such as gaus, laps
%                         kernelSize, a value such as 1 or a relation such as
%                               med for median intersample distance
%                               nnxx for mean xx-th nearest neighbor distance
%                         reg, regularization, a value such as 10^-3 or a relation such as
%                               n^(-1) for 1/(number of samples)
% Output: val is the value of hsncic
% Default: kernel is gaus, kernelSize is med and ref is n^(-1)
% Author: Sohan Seth (sohan@cnel.ufl.edu)

n = size(X,1);
if size(Y,1) ~= n | size(Z,1) ~= n
    error('Length mismatch')
end

if nargin < 3
    error('Too few arguments')
end

if nargin == 3
    param.kernel = 'gaus';
    param.kernelSize = 'med';
    param.reg = 'n^(-1)';
end

C = eye(n) - ones(n,1) * ones(1,n) /n;

Kxz = C * grammat([X,Z],param) * C; % Symmetric version
Kyz = C * grammat([Y,Z],param) * C;
Kz = C * grammat(Z,param) * C;

if ~isnumeric(param.reg)
    reg = eval(param.reg);
else
    reg = param.reg;
end

Rxz = Kxz * inv(Kxz + n*reg*eye(n));
Ryz = Kyz * inv(Kyz + n*reg*eye(n));
Rz = Kz * inv(Kz + n*reg*eye(n));

d = trace(Ryz * Rxz - 2* Ryz * Rxz * Rz + Ryz * Rz * Rxz * Rz);