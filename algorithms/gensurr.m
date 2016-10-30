function [P info] = gensurr(data,varargin)

% This function generates surrogate for assessing conditional dependence
% Input: 
%   (data) is a struct with fields,
%       Matrices (X),(Y),(Z) with rows as observations
%       (OR)
%       Matrices (distMatX),(distMatY),(distMatZ) of pairwise distances
%           and (dx),(dy),(dz) for dimensionality of (X),(Y),(Z) and (nn)
%   (optional)
%       (numDraw) for number of draws i.e. independent permutations
%                   (default) is 1
%       (nn) for the nearest neighbor value
%                   (default) is 'ceil(n^0.5)'
%       (note) nn can be expression such as 'n/2' or 'ceil(n^0.5)'
%
% Output: 
%   (P) is a cell of structs where each struct has fields (X),(Y),(Z)
%       containing the locations of the resamples
%   (note) If numDraw is 1 then P is struct not cell
%   (info) (optional) is a struct containing fields 
%           (distMatX),(distMatY),(distMatZ) and 
%           (dx),(dy),(dz) the dimensions of X,Y,Z
% 
% Author: Sohan Seth, sohan.seth@hiit.fi

% Computing distance metrices
if isfield(data,'X')
    distanceX = computel2matrix(data.X);
    [n dx] = size(data.X);
else
    if isfield(data,'distMatX')
        distanceX = data.distMatX;
        if isfield(data,'dx')
            dx = data.dx;
        else
            error('Dimension of X missing.')
        end
        error('X or distMatX must be provided.')
    end
end
if isfield(data,'Y')
    distanceY = computel2matrix(data.Y);
    [n dy] = size(data.Y);
else
    if isfield(data,'distMatY')
        distanceY = data.distMatY;
        if isfield(data,'dy')
            dy = data.dy;
        else
            error('Dimension of Y missing.')
        end
        error('Y or distMatY must be provided.')
    end
end
if isfield(data,'Z')
    distanceZ = computel2matrix(data.Z);
    [n dz] = size(data.Z);
else
    if isfield(data,'distMatZ')
        distanceZ = data.distMatZ;
        if isfield(data,'dz')
            dz = data.dz;
        else
            error('Dimension of Z missing.')
        end
        error('Z or distMatZ must be provided.')
    end
end

n = length(distanceX);


numDraw = 1;
nn = 'ceil(n^0.5)';
for count = 1:2:length(varargin)
    switch varargin{count}
        case 'numDraw'
            numDraw = varargin{count+1};
        case 'nn'
            nn = varargin{count+1};
    end
end

if ~isnumeric(nn)
    nn = eval(nn);
end

% Initialization
if numDraw == 1
    P.X = reshape([1:n]',n,1); P.Y = reshape([1:n]',n,1); P.Z = reshape([1:n]',n,1);
else
    for countDraw = 1:numDraw
        P{countDraw}.X = reshape([1:n]',n,1); P{countDraw}.Y = reshape([1:n]',n,1); P{countDraw}.Z = reshape([1:n]',n,1);
    end
end

% Computing distribution matrices
distXZ = zeros(n); distYZ = zeros(n);
for count = 1:n
    distanceXZ = sort(sqrt(distanceX.^2 + repmat(distanceZ(:,count).^2,1,n)));
    distanceYZ = sort(sqrt(distanceY.^2 + repmat(distanceZ(:,count).^2,1,n)));
    
    initXZ = sum(distanceXZ  == 0); initIndXZ = sub2ind([n,n], initXZ + nn, 1:n);
    initYZ = sum(distanceYZ  == 0); initIndYZ = sub2ind([n,n], initYZ + nn, 1:n);
    
    distXZ(:,count) = (nn./(n - initXZ)./distanceXZ(initIndXZ).^(dx+dz))';
    distYZ(:,count) = (nn./(n - initYZ)./distanceYZ(initIndYZ).^(dy+dz))';
    
    %%% The following code prevents blowup of volume due to high
    %%%     dimensionality, but it increases the complexity of the method
    %distXZlog(:,count) = (log(nn) - log(n - initXZ) - (dx+dz) * log(distanceXZ(initIndXZ)))';
    %temp = repmat(distXZlog(:,count),1,n) - repmat(distXZlog(:,count)',n,1);
    %distXZlog(:,count) = (cumsum(1./sum(exp(temp))));
    %distYZlog(:,count) = (log(nn) - log(n - initYZ) - (dy+dz) * log(distanceYZ(initIndYZ)))';
    %temp = repmat(distYZlog(:,count),1,n) - repmat(distYZlog(:,count)',n,1);
    %distYZlog(:,count) = (cumsum(1./sum(exp(temp))));
    %distXZ = distXZlog; distYZ = distYZlog;
    
end
distXZ = cumsum(distXZ ./ repmat(sum(distXZ),n,1));
distYZ = cumsum(distYZ ./ repmat(sum(distYZ),n,1));

% Sampling from distribution matrices
if numDraw == 1
    P.Z = reshape([1:n]',n,1); P.X = reshape(samplefrom(distXZ),n,1); P.Y = reshape(samplefrom(distYZ),n,1);
else
    for countDraw = 1:numDraw
        P{countDraw}.Z = reshape([1:n]',n,1); P{countDraw}.X = samplefrom(distXZ); P{countDraw}.Y = samplefrom(distYZ);
    end
end

if nargout == 2
    info.distMatX = distanceX;
    info.distMatY = distanceY;
    info.distMatZ = distanceZ;
    info.dx = dx; info.dy = dy; info.dz = dz;
    info.nn = nn;
end

function loc = samplefrom(dist)
% This function retuns a sample from a given distribution
% Input:
%   (dist) is a matrix with each column a distribution function over integers
% Output:
%   (loc) is a vector of samples from each column
n = size(dist,2);
distShifted = [zeros(1,n); dist(1:end-1,:)];
temp = rand(1,n);
[loc temp] = ind2sub([n,n], find(repmat(temp,n,1) > distShifted & repmat(temp,n,1) < dist)); loc = loc';

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