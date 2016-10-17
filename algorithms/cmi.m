function [val info] = cmi(data,varargin)

% This function computes the conditional mutual information 
% Input: 
%   (data) is a struct with fields,
%       Matrices (X),(Y),(Z) with rows as observations
%       (OR)
%       Matrices (distMatX),(distMatY),(distMatZ) of pairwise distances
%           and (dx),(dy),(dz) for dimensionality of (X),(Y),(Z) and (nn)
%   (optional)
%       (nn) for the nearest neighbor value
%                   (default) is 'ceil(n^0.5)'
%       (note) nn can be expression such as 'n/2' or 'ceil(n^0.5)'
% Output:
%   (val) is the value of CMI
%   (info) (optional) is a struct containing fields 
%           (distMatX),(distMatY),(distMatZ) and 
%           (dx),(dy),(dz) the dimensions of X,Y,Z
% Author: Sohan Seth, sohan.seth@hiit.fi

if isfield(data,'X')
    distMatX = computel2matrix(data.X);
    [n dx] = size(data.X);
else
    if isfield(data,'distMatX')
        distMatX = data.distMatX;
        if isfield(data,'dx')
            dx = data.dx;
        else
            error('Dimension of X missing.')
        end
        error('X or distMatX must be provided.')
    end
end
if isfield(data,'Y')
    distMatY = computel2matrix(data.Y);
    [n dy] = size(data.Y);
else
    if isfield(data,'distMatY')
        distMatY = data.distMatY;
        if isfield(data,'dy')
            dy = data.dy;
        else
            error('Dimension of Y missing.')
        end
        error('Y or distMatY must be provided.')
    end
end
if isfield(data,'Z')
    distMatZ = computel2matrix(data.Z);
    [n dz] = size(data.Z);
else
    if isfield(data,'distMatZ')
        distMatZ = data.distMatZ;
        if isfield(data,'dz')
            dz = data.dz;
        else
            error('Dimension of Z missing.')
        end
        error('Z or distMatZ must be provided.')
    end
end

n = length(distMatX);

nn = 'ceil(n^0.5)';
for count = 1:2:length(varargin)
    switch varargin{count}
        case 'nn'
            nn = varargin{count+1};
    end
end

if ~isnumeric(nn)
    nn = eval(nn);
end

distMatXYZ = sort(sqrt(distMatX.^2 + distMatY.^2 + distMatZ.^2));
distMatXZ = sort(sqrt(distMatX.^2 + distMatZ.^2));
distMatYZ = sort(sqrt(distMatY.^2 + distMatZ.^2));
distMatZ = sort(distMatZ);

initXYZ = sum(distMatXYZ  == 0); initIndXYZ = sub2ind([n,n], initXYZ + nn, 1:n);
initXZ = sum(distMatXZ  == 0); initIndXZ = sub2ind([n,n], initXZ + nn, 1:n);
initYZ = sum(distMatYZ  == 0); initIndYZ = sub2ind([n,n], initYZ + nn, 1:n);
initZ = sum(distMatZ  == 0); initIndZ = sub2ind([n,n], initZ + nn, 1:n);

val = mean( (dy + dz)*log(distMatYZ(initIndYZ)) + (dx + dz)*log(distMatXZ(initIndXZ)) ...
    - (dx + dy + dz)*log(distMatXYZ(initIndXYZ)) - dz*log(distMatZ(initIndZ)) ) ...
    + log( gamma((dx+dy+dz)/2 + 1) .* gamma((dz)/2 + 1) ./ gamma((dy+dz)/2 + 1) ./ gamma((dx+dz)/2 + 1)) ...
    - mean(log( (n - initXYZ) .* (n - initZ) ./ (n - initXZ) ./ (n - initYZ)));

if nargout == 2
    info.distMatX = distanceX;
    info.distMatY = distanceY;
    info.distMatZ = distanceZ;
    info.dx = dx; info.dy = dy; info.dz = dz;
    info.nn = nn;
end

function K = computel2matrix(X)
% This function computes the pairwise Euclidean distances between the rows of X
% Input:
%   (X) is a matrix with rows as observations
%       (note) (X) can be a array of cells,
% Output:
%   (K) is the matrix of pairwise distances
if any(iscell(X(:)))
    X = cell2mat(X);
end
[n d] = size(X);
innerProduct = X*X';
K = repmat(diag(innerProduct),1,n) + repmat(diag(innerProduct)',n,1) - 2*innerProduct;
K = sqrt(K);