function K = grammat(X,Y,param,varargin)

% This routine computes the Gram matrix K such that 
% K_{ij} = \kappa(x_i,y_j) where \kappa is a (strictly) positive definite kernel.

% INPUTS
% X is a (n x d) dimensional matrix
% Y is a (m x d) dimensional matrix
% param is a structure with fields kernel and kernelSize 
%       i.e. param.kernel and param.kernelSize.
% Optional parameters : 'diag', 'center',

% OUTPUT
% K is  a (n x m) matrix of kernel evaluations (distances)
%     is a (n x 1) matrix of diagonal entries of K is 'diag' is set.

% Kernels can be one of the following -- 
%       GAUSsian, LAPlacian, TRIangular, HEAVYside, Kronecker DELTA,  POLYnomial, LINear 
% of the following type --
%       RADial, Shift-Invariant and PROduct.
% Naming convention : func_type e.g. gaus_rad

% Kernel Size can be a number or choices such as 'med' and 'nnXX'.
%       'med' sets the kernel size to the median inter sample distance
%       'nnXX' sets the kernel size to the average XX-th neighbor distance
%       'var' sets the kernel size to the variance of the data
% of the samples from X.
% XX must be an integer e.g. nn20.

% 'diag' computes only the diagonal of the matrix i.e. 
%       \kappa(x_i,y_i). X and Y must have same length.
% 'center' computes the centered Gram matrix i.e. 
%       HKH where H = I - 11' is the centering matrix.

% USAGE:
%       K = grammat(X) computes Gram matrix \kappa(x_i,x_j) 
%               where \kappa is Gaussian with kernel size 1
%       K = grammat(X,Y) cmputes Gram matrix \kappa(x_i,y_j) 
%               where \kappa is Gaussian with kernel size 1
%       K = grammat(X,'func_type') computes Gram matrix \kappa(x_i,x_j) 
%               where \kappa is func_type with kernel size 1
%       K = grammat(X, param) computes Gram matrix \kappa(x_i,x_j)
%               where \kappa and kernelSize set by param
%       K = grammat(X,Y,param,option) uses option

if nargin == 1
    Y = X;
    param.kernel = 'gaus';
    param.kernelSize = 1;
end

if nargin == 2 & isnumeric(Y)
    param.kernel = 'gaus';
    param.kernelSize = 1;
end

if nargin == 2 & isstruct(Y)
    param = Y;
    Y = X;
end

if nargin == 2 & ischar(Y)
    param.kernel = Y;
    param.kernelSize = 1;
    Y = X;
end

if (nargin == 3 | nargin == 4) & ischar(param)
    temp = param; clear param;
    param.kernel = temp;
    param.kernelSize = 1;
end

if ~isfield(param,'kernelSize')
    param.kernelSize = 1;
end

diagFlag = 0;
centering = 0;
distance = 0;
if length(varargin) ~= 0
    option = varargin{1};
    if strcmp(option,'diag')
        diagFlag = 1;
    end
    if strcmp(option,'center')
        centering = 1;
    end
    if strcmp(option,'distance')
        distance = 1;
    end
end

[nx,dx] = size(X); [ny,dy] = size(Y);
if dx ~= dy
    error('X and Y must same number of columns')
end

if diagFlag
    if nx ~= ny
        error('X and Y must same number of rows')
    end
else
    K = ones(nx,ny);
end

switch param.kernel
    case 'gaus' % k(x,y) = prod_i exp(-(x_i-y_i)^2/s^2)
        if isnumeric(param.kernelSize)
            if length(param.kernelSize) == dx
                kernelSize = param.kernelSize;
            else
                kernelSize = param.kernelSize * ones(1,dx);
            end
        else
            switch param.kernelSize(1:2)
                case 'me'
                    kernelSize = med(X) * ones(1,dx);
                case 'va'
                    kernelSize = myvar(X);
                case 'nn'
                    kernelSize = nn(X,str2num(param.kernelSize(3:end))) * ones(1,dx);
                case 'n^'
                    n = size(X,1);
                    kernelSize = eval(param.kernelSize) * ones(1,dx);
                otherwise
                    error('Invalid kernel size')
            end
        end
        if diagFlag
            K = ones(nx,1);
        else
            for count = 1:dx
                K = K .* exp(-( repmat(X(:,count),1,ny) - repmat((Y(:,count))',nx,1) ).^2/kernelSize(count)^2);
            end
        end

    case 'gauf' % k(x,y) = prod_i (3 - (x_i - y_i)^2) *exp(-(x_i-y_i)^2/s^2) / 2
        if isnumeric(param.kernelSize)
            if length(param.kernelSize) == dx
                kernelSize = param.kernelSize;
            else
                kernelSize = param.kernelSize * ones(1,dx);
            end
        else
            switch param.kernelSize(1:2)
                case 'me'
                    kernelSize = med(X) * ones(1,dx);
                case 'va'
                    kernelSize = myvar(X);
                case 'nn'
                    kernelSize = nn(X,str2num(param.kernelSize(3:end))) * ones(1,dx);
                case 'n^'
                    n = size(X,1);
                    kernelSize = eval(param.kernelSize) * ones(1,dx);
                otherwise
                    error('Invalid kernel size')
            end
        end
        if diagFlag
            K = 3*ones(nx,1)/2;
        else
            for count = 1:dx
                K = K .* exp(-( repmat(X(:,count),1,ny) - repmat((Y(:,count))',nx,1) ).^2/kernelSize(count)^2) ...
                    .* (3 - ( repmat(X(:,count),1,ny) - repmat((Y(:,count))',nx,1) ).^2)/2;
            end
        end

    case 'tris' % k(x,y) = prod_i (1 - |x_i-y_i|/s) I(|x_i - y_i|<s)
        if isnumeric(param.kernelSize)
            if length(param.kernelSize) == dx
                kernelSize = param.kernelSize;
            else
                kernelSize = param.kernelSize * ones(1,dx);
            end
        else
            switch param.kernelSize(1:2)
                case 'me'
                    kernelSize = med(X) * ones(1,dx);
                case 'va'
                    kernelSize = myvar(X);
                case 'nn'
                    kernelSize = nn(X,str2num(param.kernelSize(3:end))) * ones(1,dx);
                case 'n^'
                    n = size(X,1);
                    kernelSize = eval(param.kernelSize) * ones(1,dx);
                otherwise
                    error('Invalid kernel size')
            end
        end
        if diagFlag
            K = ones(nx,1);
        else
            for count = 1:dx
                temp = abs(repmat(X(:,count),1,ny) - repmat((Y(:,count))',nx,1)) / kernelSize(count);
                temp(temp > 1) = 1;
                eval(['K = K .* (1 - temp);']);
            end
        end

    case 'poly' % k(x,y) = (1 + <x,y>)^s
        if isnumeric(param.kernelSize)
            kernelSize = param.kernelSize;
        end
        if mod(kernelSize * 10,10) ~= 0
            error('Order of the polynomial kernel must be integer')
        end
        if diagFlag
            K = (1 + sum(X.*X,2)).^kernelSize;
        else
            K = (1 + X * Y').^kernelSize;
        end

    case 'laps' % k(x,y) = prod_i exp(-|x_i-y_i|/s)
        if isnumeric(param.kernelSize)
            if length(param.kernelSize) == dx
                kernelSize = param.kernelSize;
            else
                kernelSize = param.kernelSize * ones(1,dx);
            end
        else
            switch param.kernelSize(1:2)
                case 'me'
                    kernelSize = med(X) * ones(1,dx);
                case 'va'
                    kernelSize = myvar(X);
                case 'nn'
                    kernelSize = nn(X,str2num(param.kernelSize(3:end))) * ones(1,dx);
                case 'n^'
                    n = size(X,1);
                    kernelSize = eval(param.kernelSize) * ones(1,dx);
                otherwise
                    error('Invalid kernel size')
            end
        end
        if diagFlag
            K = ones(nx,1);
        else
            for count = 1:dx
                K = K .* exp(-abs( repmat(X(:,count),1,ny) - repmat((Y(:,count))',nx,1) )/kernelSize(count));
            end
        end

    case 'l2' % k(x,y) = sqrt (sum_i (x_i - y_i)^2) % This is not a kernel but distance (supporting function)
        if diagFlag
            K = zeros(nx,1);
        else
            K = zeros(nx,ny);
            for count = 1:dx
                K = K + ( repmat(X(:,count),1,ny) - repmat((Y(:,count))',nx,1)).^2;
            end
            K = sqrt(K);
        end

    case 'idfs' % k(x,y) = prod_i I(x_i<y_i), I(a) = 1 if a is true or 0 otherwise
        if diagFlag
            K = prod(X <= Y,2);
        else
            for count = 1:dx
                K = K .* (repmat(X(:,count),1,ny) <= repmat((Y(:,count))',nx,1));
            end
        end

    case 'dels' % k(x,y) = prod_i I(x_i=y_i), I(a) = 1 if a is true or 0 otherwise
        if diagFlag
            K = prod(X == Y,2);
        else
            for count = 1:dx
                K = K .* (repmat(X(:,count),1,ny) == repmat((Y(:,count))',nx,1));
            end
        end

    otherwise
        error('Invalid kernel')
end

if centering == 1
    if nx ~= ny
        error('centering not defined : not square matrix')
        H = eye(nx) - ones(nx,1) * ones(1,nx)/nx;
        K = H * K * H;
    end
end

function kernelSizeX = med(X)

param.kernel = 'l2';
Kxx = grammat(X,param);
Kxx = Kxx(:);
kernelSizeX = quantile(Kxx(Kxx>0),0.5);

function kernelSizeX = myvar(X)

kernelSizeX = var(X);

function kernelSizeX = nn(X,num)

param.kernel = 'l2';
Kxx = sort(grammat(X,param));
kernelSizeX = mean(Kxx(min(num + 1,size(X,1)),:));