function [ X, Y, Z ] = genCItestData( numSamps, m, arch, depType, rvType, xyDep, varargin )
%GENCITESTDATA Generates test data for {X indep Y | Z} tests with various 
%options. See detailed documentation for more information.
% Inputs:
%  numSamps - the number of samples to generate
%  m - the number of Z nodes
%  arch - the architecture type.  See documentation below for more details.
%  depType - the dependency type to be induced in the architecture.
%            Currently, the supported options are:
%            - 'Gaussian' - all dependencies are modeled w/ Gaussian copula
%            - 'Mixed' - dependencies are modeled w/ different kinds of
%                        copulas
%  rvType - the type of RV's.  Currently, the supported options are:
%           - 'continuous' - all RV's are continuous
%           - 'discrete' - all RV's are discrete
%           - 'mixed' - a mix of continuous and discrete RVs
%  xyDep - if 1, we induce a direct dependency between X and Y, 
%          if 0, no direct dependency between X and Y is induced.
%  varargin -
%             varargin{1} - the specific case within the dependency type,
%                           if empty, this is defaulted to 1.
%             varargin{2} - the specific case within the random variable 
%                           type, if empty, this is defaulted to 1
% 
% Outputs:
%  X - the X vector, of dimension [numSamps x 1]
%  Y - the Y vector, of dimension [numSamps x 1]
%  Z - the Z vector, of dimension [numSamps x m]

depType_LC = lower(depType);
rvType_LC = lower(rvType);

depTypeErrorStr = 'Dependency Types other than Gaussian or Mixed are currently unsupported!';
rvTypeErrorStr = 'Random Variable Types other than continuous, discrete, or mixed are currently unsupported!';
depTypeCaseErrorStr = 'Dependency Type Case currently unsuppored!';

% Parse input arguments
nVarargs = length(varargin);
depTypeCase = 1;
rvTypeCase = 1;
if(nVarargs>0)
    depTypeCase = varargin{1};
end
if(nVarargs>1)
    rvTypeCase = varargin{2};
end

switch m
    case 1
        switch depType_LC
            case 'gaussian'
                
            case 'mixed'
            otherwise
                error(depTypeErrorStr);
        end
    case 2
        switch arch
            case 1
                switch depType_LC
                    case 'gaussian'
                    case 'mixed'
                    otherwise
                        error(depTypeErrorStr);
                end
            case 2
                switch depType_LC
                    case 'gaussian'
                    case 'mixed'
                    otherwise
                        error(depTypeErrorStr);
                end
            case 3
                switch depType_LC
                    case 'gaussian'
                    case 'mixed'
                    otherwise
                        error(depTypeErrorStr);
                end
            otherwise
                error('For m=2, arch>3 unsupported!\n');
        end
    case 3
        switch arch
            case 1
                switch depType_LC
                    case 'gaussian'
                    case 'mixed'
                    otherwise
                        error(depTypeErrorStr);
                end
            case 2
                switch depType_LC
                    case 'gaussian'
                    case 'mixed'
                    otherwise
                        error(depTypeErrorStr);
                end
            case 3
                switch depType_LC
                    case 'gaussian'
                    case 'mixed'
                    otherwise
                        error(depTypeErrorStr);
                end
            case 4
                switch depType_LC
                    case 'gaussian'
                    case 'mixed'
                    otherwise
                        error(depTypeErrorStr);
                end
            case 5
                switch depType_LC
                    case 'gaussian'
                    case 'mixed'
                    otherwise
                        error(depTypeErrorStr);
                end
            case 6
                switch depType_LC
                    case 'gaussian'
                    case 'mixed'
                    otherwise
                        error(depTypeErrorStr);
                end
            otherwise
                error('For m=3, arch>6 unsupported!\n');
        end
    otherwise
        error('m > 3 unsupported!\n');
end

% add RV info to the desired structure

% add dependency between X and Y to test for Type II errors if desired
if(xyDep)
    ff = normrnd(0,1,numSamps,1) * 0.5;
    X = X + ff; Y = Y + ff;
end

end