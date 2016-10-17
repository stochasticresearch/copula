function val = cassor(data,varargin)

% This function computes the "conditional association" of Y with X given Z
% (note) conditional association is not symmetric,
%           i.e. (Y,X)|Z is not the same as (X,Y)|Z. 
%           So, conditional association matrix of all triplets will be asymmetric.
%
% Input: 
%   (data) is a structure with fields
%       Matrices (X),(Y), and (Z) with observations on each row
%           (note) if only (X) is provided, cassor computes the conditional
%                       association between all pairs given the rest
%       (and)
%       Structures (distFuncX), (distFuncY), and (distFuncZ) with fields,
%           Function handle name e.g. distFuncX.name = @distFunc;
%               (default) is the l2 distance metric, and
%           associated parameters as distFuncX.param = param;            
%               (note) parameters are simply passed to the @distFunc; 
%                       cassor does not use it explicitly.
%           (note) If parameters not needed, distFuncX can be the distance handle 
%       (OR)
%       Matrices (distMatX),(distMatY), and (distMatZ) of pairwise distances
%           (note) Observations and pairwise distances can be mixed
%           (note) In case both provided for a variable, observations have preference
%
% Output:
%   (val) is the value of the conditional association
%       (or) a matrix of pairwise conditional association given rest of variables
%   (info) (optional) is a struct with fields
%       Matrices (distMatX),(distMatY), and (distMatZ) of pairwise distances
%
% Optional:
%   (embed) is embedding dimension
%       (note) if this is provided then the observations are treated as time series 
%               and conditional association is evaluated between embedded vector.
%
% Example 1:
%   data.X = randn(100,1); data.Y = randn(100,1); data.Z = randn(100,1);
%   val = cassor(data);
% Example 2:
%   data.X = dgp7713(100,0.8);
%   val = cassp(data,'embed',10);
%
% Author - Sohan Seth (sohan.seth@hiit.fi)

embed = 0;
for count = 1:2:length(varargin)
    switch varargin{count}
        case 'embed'
            embed = varargin{count+1};
    end
end

pairwiseCondAsso = 0;
if ((~isfield(data,'Y') & ~isfield(data,'distMatY')) || (~isfield(data,'Z') & ~isfield(data,'distMatZ')))
    if  size(data.X,2)<3
        if embed
            pairwiseCondAsso = 1;
        else
            error('Either X should have three columns or Y and Z should be provided')
        end
    else
        pairwiseCondAsso = 1;
    end
end

if ~pairwiseCondAsso
    if isstruct(data)
        if isfield(data,'X')
            X = data.X;
            if ~isfield(data,'distFuncX')
                distMatX = computel2matrix(X);
            else
                if ~isstruct(data.distFuncX)
                    distMatX = data.distFuncX(X) + diag(Inf*ones(size(X,1),1));
                else
                    distMatX = data.distFuncX.name(X,data.distFuncX.param) + diag(Inf*ones(size(X,1),1));
                end
            end
        else
            if isfield(data, 'distMatX')
                distMatX = data.distMatX;
            else
                error('X field must be provided in the struct data')
            end
        end
        if isfield(data,'Y')
            Y = data.Y;
            if ~isfield(data,'distFuncY')
                distMatY = computel2matrix(Y);
            else
                if ~isstruct(data.distFuncY)
                    distMatY = data.distFuncY(Y) + diag(Inf*ones(size(Y,1),1));
                else
                    distMatY = data.distFuncY.name(Y,data.distFuncY.param) + diag(Inf*ones(size(Y,1),1));
                end
            end
        else
            if isfield(data, 'distMatY')
                distMatY = data.distMatY;
            else
                error('Y field must be provided in the struct data')
            end
        end
        if isfield(data,'Z')
            Z = data.Z;
            if ~isfield(data,'distFuncZ')
                distMatZ = computel2matrix(Z);
            else
                if ~isstruct(data.distFuncZ)
                    distMatZ = data.distFuncZ(Z) + diag(Inf*ones(size(Z,1),1));
                else
                    distMatZ = data.distFuncZ.name(Z,data.distFuncZ.param) + diag(Inf*ones(size(Z,1),1));
                end
            end
        else
            if isfield(data, 'distMatZ')
                distMatZ = data.distMatZ;
            else
                error('Z field must be provided in the struct data')
            end
        end
    else
        error('Input should be a structure. See help.')
    end
    
    n = size(distMatX,1);
    
    rankX = dist2rank(distMatX);
    rankY = dist2rank(distMatY);
    rankZ = dist2rank(distMatZ);
    
    clear data
    
    [valTemp nnZ] = min(rankZ);
    rankXZ = rankX + rankZ;
    [sortVal sortLoc] = sort(rankXZ);
    nnXZ = sortLoc(1,:);
    % Checking for idential values
    % for count = 1:n
    %     simLoc = sortLoc(sortVal(:,count) == sortVal(1,count),count);
    %     if length(simLoc) > 1
    %         temp = Inf*ones(n,1); temp(simLoc) = 0;
    %         [valTemp nnXZ(1,count)] = min(rankZ(:,count) + temp);
    %     end
    % end
    
    clear rankX rankZ rankXZ temp valTemp simLoc sortVal sortLoc
    
    rankAll = zeros(n,1);
    for sampleCount = 1:n
        indexTemp = nnZ(sampleCount);
        rankAll(sampleCount) = rankY(indexTemp,sampleCount);
    end
    
    empPDF = zeros(n-1,1);
    for sampleCount = 1:n-1
        empPDF(sampleCount) = mean(rankAll == sampleCount);
    end
    valZ = sum((n - [1:n-1]') .* empPDF) / length(empPDF);
    
    rankAll = zeros(n,1);
    for sampleCount = 1:n
        indexTemp = nnXZ(sampleCount);
        rankAll(sampleCount) = rankY(indexTemp,sampleCount);
    end
    
    empPDF = zeros(n-1,1);
    for sampleCount = 1:n-1
        empPDF(sampleCount) = mean(rankAll == sampleCount);
    end
    valXZ = sum((n - [1:n-1]') .* empPDF) / length(empPDF);
    
    val = valXZ - valZ;
    
    if nargout == 2
        info.distMatX = distMatX;
        info.distMatY = distMatY;
        info.distMatZ = distMatZ;
    end
else
    W = data.X;
    [n d] = size(W);
    val = zeros(d) + NaN;
    for count1 = 1:d
        for count2 = 1:d
            if count1 ~= count2
                if ~sum(embed)
                    X = W(:,count1);
                    Y = W(:,count2);
                    Z = W; Z(:,[count1,count2]) = [];
                else
                    Y = W(embed+1:n,count2);
                    if any(iscell(W(:)))
                        X = cell(n - embed, embed); Z = cell(n - embed,(d-1)*embed);
                    else
                        X = zeros(n - embed, embed); Z = zeros(n - embed,(d-1)*embed);
                    end
                    for countLag = 1:embed
                        X(:,countLag) = W(embed - countLag + 1:n - countLag,count1);
                        tempZ = W(embed - countLag + 1:n - countLag,:);
                        tempZ(:,count1) = [];
                        Z(:,(countLag -1) * (d-1)+ 1:countLag * (d-1)) = tempZ;
                    end
                    %X = zscore(X); Y = zscore(Y); Z = zscore(Z); 
                end
                if ~isfield(data,'distFuncX')
                    distMatX = computel2matrix(X);
                    distMatY = computel2matrix(Y);
                    distMatZ = computel2matrix(Z);
                else
                    if ~isstruct(data.distFuncX)
                        distMatX = data.distFuncX(X) + diag(Inf*ones(size(X,1),1));
                        distMatY = data.distFuncX(Y) + diag(Inf*ones(size(X,1),1));
                        distMatZ = data.distFuncX(Z) + diag(Inf*ones(size(X,1),1));
                    else
                        distMatX = data.distFuncX.name(X,data.distFuncX.param) + diag(Inf*ones(size(X,1),1));
                        distMatY = data.distFuncX.name(Y,data.distFuncX.param) + diag(Inf*ones(size(X,1),1));
                        distMatZ = data.distFuncX.name(Z,data.distFuncX.param) + diag(Inf*ones(size(X,1),1));
                    end
                end
                
                n = size(distMatX,1);
                
                rankX = dist2rank(distMatX);
                rankY = dist2rank(distMatY);
                rankZ = dist2rank(distMatZ);                
                
                [valTemp nnZ] = min(rankZ);
                rankXZ = rankX + rankZ;
                [sortVal sortLoc] = sort(rankXZ);
                nnXZ = sortLoc(1,:);
                % Checking for idential values
                % for count = 1:n
                %     simLoc = sortLoc(sortVal(:,count) == sortVal(1,count),count);
                %     if length(simLoc) > 1
                %         temp = Inf*ones(n,1); temp(simLoc) = 0;
                %         [valTemp nnXZ(1,count)] = min(rankZ(:,count) + temp);
                %     end
                % end
                
                clear rankX rankZ rankXZ temp valTemp simLoc sortVal sortLoc
                
                rankAll = zeros(n,1);
                for sampleCount = 1:n
                    indexTemp = nnZ(sampleCount);
                    rankAll(sampleCount) = rankY(indexTemp,sampleCount);
                end
                
                empPDF = zeros(n-1,1);
                for sampleCount = 1:n-1
                    empPDF(sampleCount) = mean(rankAll == sampleCount);
                end
                valZ = sum((n - [1:n-1]') .* empPDF) / length(empPDF);
                
                rankAll = zeros(n,1);
                for sampleCount = 1:n
                    indexTemp = nnXZ(sampleCount);
                    rankAll(sampleCount) = rankY(indexTemp,sampleCount);
                end
                
                empPDF = zeros(n-1,1);
                for sampleCount = 1:n-1
                    empPDF(sampleCount) = mean(rankAll == sampleCount);
                end
                valXZ = sum((n - [1:n-1]') .* empPDF) / length(empPDF);
                
                val(count1,count2) = valXZ - valZ;
            end
        end
    end
    
    if nargout == 2
        info = [];
        warning('Too many matrices');
    end
end

function rank = dist2rank(distMat)
% This function converts distance matrix to rank matrix
% Input: (distMat) is pairwise distance matrix
% Output: (rank) is pairwise rank matrix
distMat = distMat+ diag(Inf*ones(length(distMat),1));
[valTemp rank] = sort(distMat);
[valTemp rank] = sort(rank);

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