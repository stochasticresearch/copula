classdef hcbn < handle
    %HCBN Definition of a Hybrid Copula Bayesian Network
    %
    %**********************************************************************
    %* 
    %* Copyright (C) 2016  Kiran Karra <kiran.karra@gmail.com>
    %*
    %* This program is free software: you can redistribute it and/or modify
    %* it under the terms of the GNU General Public License as published by
    %* the Free Software Foundation, either version 3 of the License, or
    %* (at your option) any later version.
    %*
    %* This program is distributed in the hope that it will be useful,
    %* but WITHOUT ANY WARRANTY; without even the implied warranty of
    %* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    %* GNU General Public License for more details.
    %*
    %* You should have received a copy of the GNU General Public License
    %* along with this program.  If not, see <http://www.gnu.org/licenses/>
    %* 
    %**********************************************************************

    properties
        dag;        % the Directed Acyclic Graph structure.  The format of
                    % the DAG is an adjancency matrix.  Element (i,j)
                    % represents a connection from the ith random variable
                    % to the jth random variable
                    % TODO: consider changing this to a sparse matrix
                    % representation
        nodeNames;  % a cell array of the list of node names  
        nodeVals;   % an array of the topological order of the nodes w.r.t
                    % the DAG
        D;          % the number of nodes in the graph
        empInfo;    % a cell array of the empirical information for each
                    % node.  empInfo{i} is the empirical information for
                    % the node nodeNames{nodeVals{i}}.  Empirical
                    % information consits of 3 things, 1.) the domain of
                    % the random variable, the empirical marginal density 
                    % function for node i, and the empirical marginal 
                    % distribution function for node i
        copulaFamilies; % a cell array of the copula of the family associated 
                        % with nodeNames{nodeVals{i}}
        X;          % an M x D matrix of the observable data.  This can be
                    % considered the "training" data.
        X_xform;    % a transformed M x D matrix of the observable data,
                    % where we "continue" the discrete random variables, as
                    % described in Neslehova's paper, and Denuit and
                    % Lambert's paper
        K;          % the number of points over the unit-hypercube for
                    % which to calculate the empirical copula.  The higher
                    % the value of K, the greater the accuracy, but also
                    % the higher the memory and processing time
                    % requirements
        discreteNodes;  % the names of the nodes which should be modeled as
                        % discrete random variables
        discNodeIdxs;   % the indices of of the nodes of the discrete
                        % random variables
        
        h;          % the beta-kernel density estimation weighting
                        
        LOG_CUTOFF;
                                
        SIM_NUM;    % used for debugging
        TYPE_NAME;  % used for debugging
        
        PSEUDO_OBS_CALC_METHOD;     % should be either RANK or ECDF
    end
    
    methods
        function obj = hcbn(bntPath, X, nodeNames, discreteNodes, K, h, varargin)
            % HCBN - Constructs a HCBN object
            %  Inputs:
            %   X - a M x D matrix of the the observable data which the
            %       HCBN will model
            %   nodes - a cell array of names which represent the columns,
            %           where the ith 
            %  Optional Inputs:
            %   dag - an adjacency matrix representing the DAG structure.
            %         If this is not specified, it is assumed that the user
            %         will learn the structure through an available
            %         learning algorithm
            %
            %  TODO
            %   [ ] - 
            %
            
            obj.SIM_NUM = 0;
            obj.TYPE_NAME = '';
            obj.LOG_CUTOFF = 10^-5;
            
            obj.PSEUDO_OBS_CALC_METHOD = 'RANK';    % can be either RANK or ECDF
            
            % add BNT to the path
            addpath(genpath(bntPath));
            
            obj.D = size(X,2);      
            
            obj.nodeNames = nodeNames;
            obj.nodeVals = 1:obj.D;
            obj.copulaFamilies = cell(1,obj.D);
            obj.empInfo = cell(1,obj.D);
            
            obj.dag = zeros(obj.D,obj.D);
            
            obj.K = K;
            obj.h = h;
            
            obj.X = X;
            obj.X_xform = X;

            obj.discreteNodes = discreteNodes;
            obj.discNodeIdxs = zeros(1,length(obj.discreteNodes));
            for ii=1:length(obj.discNodeIdxs)
                nodeName = obj.discreteNodes{ii};
                for jj=1:length(obj.nodeNames)
                    if(isequal(nodeName,obj.nodeNames{jj}))
                        obj.discNodeIdxs(ii) = jj;
                        break;
                    end
                end
            end
            
            % continue the discrete random variables.  X* = X + (U - 1)
            for idx=obj.discNodeIdxs
                obj.X_xform(:,idx) = continueRv(obj.X(:,idx));
            end
            
            obj.calcEmpInfo();
            
            nVarargs = length(varargin);
            if(nVarargs>0)
                candidateDag = varargin{1};
                if(~acyclic(candidateDag))
                    error('Specified DAG is not acyclic!\n');
                end
                if(nVarargs>4)
                    opt1 = varargin{2};
                    if(strcmpi(opt1,'copulaFamilyInput'))
                        obj.setDag(candidateDag, 0);    % do NOT compute estimated copula families
                                                        % b/c they were provided as debugging
                                                        % input
                        copulaFamiliesInput = varargin{3};
                        obj.inflateCopulaInformation(copulaFamiliesInput);
                    elseif(strcmpi(opt1,'empInfoInput'))
                        % overwrite the empirical info w/ the provided
                        % debugging empirical info
                        empInfoInput = varargin{3};
                        obj.inflateEmpInfo(empInfoInput);
                    else
                        error('Unrecognized var-input!');
                    end
                    opt2 = varargin{4};
                    if(strcmpi(opt2,'copulaFamilyInput'))
                        obj.setDag(candidateDag, 0);    % do NOT compute estimated copula families
                                                        % b/c they were provided as debugging
                                                        % input
                        copulaFamiliesInput = varargin{5};
                        obj.inflateCopulaInformation(copulaFamiliesInput);
                    elseif(strcmpi(opt2,'empInfoInput'))
                        % overwrite the empirical info w/ the provided
                        % debugging empirical info
                        empInfoInput = varargin{5};
                        obj.inflateEmpInfo(empInfoInput);
                    else
                        error('Unrecognized var-input');
                    end
                elseif(nVarargs>2)
                    opt1 = varargin{2};
                    if(strcmpi(opt1,'copulaFamilyInput'))
                        obj.setDag(candidateDag, 0);    % do NOT compute estimated copula families
                                                        % b/c they were provided as debugging
                                                        % input
                        copulaFamiliesInput = varargin{3};
                        obj.inflateCopulaInformation(copulaFamiliesInput);
                    elseif(strcmpi(opt1,'empInfoInput'))
                        % overwrite the empirical info w/ the provided
                        % debugging empirical info
                        empInfoInput = varargin{3};
                        obj.inflateEmpInfo(empInfoInput);
                        obj.setDag(candidateDag);
                    else
                        error('Unrecognized var-input!');
                    end
                else
                    obj.setDag(candidateDag);       % estimate the copula families
                end
            end
        end
        
        function [] = inflateEmpInfo(obj, empInfoInput)
            for dd=1:obj.D
                obj.empInfo{dd} = empInfoInput{dd};
            end
        end
        
        function [] = inflateCopulaInformation(obj, copulaFamiliesInput)
            % inflate the copula families argument input into the
            % required format for processing w/ HCBN.  The input format
            % of copulaFamiliesInput is a cell array of dimension 
            % [1 x D].  copulaFamiliesInput{ii}{1} is a string of the
            % copula TYPE for node ii, and copulaFamiliesInput{ii}{2}
            % is the dependency parameter for that copula type.  If a
            % node is not dependent upon any other nodes, then
            % copulaFamiliesInput{ii} = [].  The dependency
            % parameter details are as follows -- if it is an
            % archimedean copula, then because all dimensions in an
            % archimedean copula are the same dependency, a scalar
            % is provided.  In the case of the Gaussian copula, the
            % input MUST be permuted such that the first column in
            % the correlation matrix refers to the child, the 2nd
            % to the first parent, the 3rd to the second parent etc
            % etc ... so in the 3-D copula case, it should be:
            %  |1 Rho(child,parent_1) Rho(child,parent_2)   |
            %  |                                            |
            %  |Rho(parent1, child)  1 Rho(parent1, parent2)|
            %  |                                            |
            %  |Rho(parent2, child) Rho(parent2, parent1) 1 |
            for dd=1:obj.D
                if(isempty(copulaFamiliesInput{dd}))
                    obj.copulaFamilies{dd} = [];
                else
                    node = obj.nodeNames{dd};
                    nodeIdx = obj.nodeVals(dd);
                    [parentIdxs, parentNames] = obj.getParents(nodeIdx);

                    u = linspace(0,1,obj.K);
                    % generate the points over which the copulas
                    % will be computed
                    familyD = length(parentIdxs) + 1;
                    sz = ones(1,familyD)*obj.K;
                    ndgridInput = cell(1,familyD);
                    for family_dd=1:familyD
                        ndgridInput{family_dd} = u;
                    end
                    ndgridOutput = cell(1,numel(ndgridInput));
                    [ndgridOutput{:}] = ndgrid(ndgridInput{:});
                    gridPoints_all = zeros(numel(ndgridOutput{1}), familyD);
                    for family_dd=1:familyD
                        gridPoints_all(:,family_dd) = reshape(ndgridOutput{family_dd},numel(ndgridOutput{family_dd}),1);
                    end

                    copulaType = copulaFamiliesInput{dd}{1};
                    copulaDep_all  = copulaFamiliesInput{dd}{2};
                    copulaDep_parents = copulaDep_all;
                    % compute PDF of family copula
                    if(strcmpi(copulaType, 'Gaussian'))
                        c = reshape(copulapdf('Gaussian', gridPoints_all, copulaDep_all), sz);
                        copulaDep_parents = copulaDep_all(2:end,2:end); % remove child correlation from corr matrix
                    elseif(strcmpi(copulaType, 'Frank'))
                        c = reshape(frankcopulapdf(gridPoints_all, copulaDep_all), sz);
                    elseif(strcmpi(copulaType, 'Gumbel'))
                        c = reshape(gumbelcopulapdf(gridPoints_all, copulaDep_all), sz);
                    elseif(strcmpi(copulaType, 'Clayton'))
                        c = reshape(claytoncopulapdf(gridPoints_all, copulaDep_all), sz);
                    else
                        error('Unsupported copula!');
                    end
                    % NOTE: for the archimedean copulas, since the
                    % copula dependency parameter remains the
                    % same,across all the dimensions,
                    % copulaDep_parents = copulaDep_all

                    % compute CDF of family copula
                    C = c;
                    for family_dd=1:familyD
                        C = cumtrapz(u, C, family_dd);
                    end

                    % integrate discrete dimensions out
                    allIdxs = [nodeIdx parentIdxs];
                    [~,discreteDimensions,~] = intersect(allIdxs,obj.discNodeIdxs); discreteDimensions = discreteDimensions';
                    C_discrete_integrate = c;
                    for discreteDimension=discreteDimensions
                        C_discrete_integrate = cumtrapz(u,C_discrete_integrate,discreteDimension);
                    end

                    if(length(parentIdxs)==1)
                        % the density will be used directly, so there
                        % is no need to calculate these copula's (they
                        % don't actually make sense because there is
                        % only one dimension, no concept of "joint")
                        C_parents = [];
                        c_parents = [];
                        C_parents_discrete_integrate = [];
                    else
                        parentsD = length(parentIdxs);
                        sz = ones(1,parentsD)*obj.K;
                        ndgridInput = cell(1,parentsD);
                        for parents_dd=1:parentsD
                            ndgridInput{parents_dd} = u;
                        end
                        ndgridOutput = cell(1,numel(ndgridInput));
                        [ndgridOutput{:}] = ndgrid(ndgridInput{:});
                        gridPoints_parents = zeros(numel(ndgridOutput{1}), parentsD);
                        for parents_dd=1:parentsD
                            gridPoints_parents(:,parents_dd) = reshape(ndgridOutput{parents_dd},numel(ndgridOutput{parents_dd}),1);
                        end

                        % compute PDF of parents copula
                        if(strcmpi(copulaType, 'Gaussian'))
                            c_parents = reshape(copulapdf('Gaussian', gridPoints_parents, copulaDep_parents), sz);
                        elseif(strcmpi(copulaType, 'Frank'))
                            c_parents = reshape(frankcopulapdf(gridPoints_parents, copulaDep_parents), sz);
                        elseif(strcmpi(copulaType, 'Gumbel'))
                            c_parents = reshape(gumbelcopulapdf(gridPoints_parents, copulaDep_parents), sz);
                        elseif(strcmpi(copulaType, 'Clayton'))
                            c_parents = reshape(claytoncopulapdf(gridPoints_parents, copulaDep_parents), sz);
                        else
                            error('Unsupported copula!');
                        end
                        % compute CDF of parents copula
                        C_parents = c_parents;
                        for parent_dd=1:parentsD
                            C_parents = cumtrapz(u, C_parents, parent_dd);
                        end

                        % integrate discrete dimensions out
                        [~,discreteDimensions,~] = intersect(parentIdxs,obj.discNodeIdxs); discreteDimensions = discreteDimensions';
                        C_parents_discrete_integrate = c_parents;
                        for discreteDimension=discreteDimensions
                            C_parents_discrete_integrate = cumtrapz(u,C_parents_discrete_integrate,discreteDimension);
                        end

                    end
                    % store into obj.copulaFamilies{dd}
                    copFam = hcbnfamily(node, nodeIdx, parentNames, parentIdxs, ...
                            C, c, C_discrete_integrate, C_parents, c_parents, C_parents_discrete_integrate);
                    obj.copulaFamilies{nodeIdx} = copFam;
                end
            end
        end
        
        function [] = setSimNum(obj, simnum)
            obj.SIM_NUM = simnum;
        end
        
        function [] = setTypeName(obj, typeName)
            obj.TYPE_NAME = typeName;
        end
        
        function [] = calcEmpInfo(obj)
            %CALCEMPINFO - calculates the empirical distribution function
            %              and the empirical density function via kernel
            %              based methods
            for ii=1:obj.D
                % check if this is a discrete or continuous node for
                % density estimation, we handle these separately
                isdiscrete = 0;
                if(any(obj.discNodeIdxs==ii))
                    isdiscrete = 1;
                end
                [F,x] = empcdf(obj.X(:,ii), isdiscrete);
                f = emppdf(obj.X(:,ii), isdiscrete);
                empInfoObj = rvEmpiricalInfo(x, f, F, isdiscrete);
                obj.empInfo{ii} = empInfoObj;
            end
        end
        
        function [parentIdxs, parentNames] = getParents(obj, node)
            %GETPARENTS - returns the indices and the names of a nodes
            %                 parents
            %
            % Inputs:
            %  node - the node index or name of the node for which the
            %         parents are desired
            %
            % Outputs:
            %  parentIdxs - a vector of the indices of all the parents of
            %               this node
            %  parentNames - a cell array of the names of the parents of
            %               this node
            
            if(ischar(node))
                % node name was provided
                nodeName = node;
                for ii=1:obj.D
                    if(isequal(nodeName, obj.nodeNames{ii}))
                        nodeIdx = ii;
                        break;
                    end
                end
            else
                % assume node index was provided
                nodeIdx = node;
            end
            
            % find the node's parents
            parentIdxs = find(obj.dag(:,nodeIdx))';
            parentNames = cell(1,length(parentIdxs));
            for jj=1:length(parentIdxs)
                parentNames{jj} = obj.nodeNames{parentIdxs(jj)};
            end
        end
        
        function [] = setDag(obj, candidateDag, varargin)
            % if you set varargin=1, then the copula families will not be
            % estimated automatically
            obj.dag = candidateDag;
            nVarargin = length(varargin);
            if(nVarargin==0 || varargin{1}~=0)
                obj.estFamilyCopula();
            end
        end
        
        function [] = estFamilyCopula(obj)
            %ESTFAMILYCOPULA - estimates the copulas for each family (child
            %                  and parents), based on the current definition
            %                  of the dag, which is stored in obj.dag
            
            % for each node, estimate the copula of that node and it's
            % parents
            for ii=1:obj.D
                node = obj.nodeNames{ii};
                nodeIdx = obj.nodeVals(ii);
                [parentIdxs, parentNames] = obj.getParents(nodeIdx);
                
                if(isempty(parentIdxs))
                    % no parents situation
                    obj.copulaFamilies{nodeIdx} = [];
                else
                    % grab the appropriate values 
                    X_in = zeros(size(obj.X_xform,1), 1+length(parentNames));
                    X_in(:,1) = obj.X_xform(:,nodeIdx);
                    kk = 2:2+length(parentIdxs)-1;
                    X_in(:,kk) = obj.X_xform(:,parentIdxs);
                    
                    % fit all continuous and hybrid models w/ an empirical
                    % copula
                    if(strcmpi(obj.PSEUDO_OBS_CALC_METHOD, 'ecdf'))
                        U_in = pobs(X_in, 'ecdf', 100);
                    else
                        U_in = pobs(X_in);
                    end

                    c = empcopulapdf(U_in, obj.h, obj.K, 'betak');
                    C = c; u = linspace(0,1,obj.K);
                    for dim=1:1+length(parentNames)
                        C = cumtrapz(u,C,dim);
                    end
                    
                    allIdxs = [nodeIdx parentIdxs];
                    % find which dimensions are discrete, and integrate
                    % those out of c.  This is the equivalent of takign the
                    % partial derivative of the copula function w.r.t. only
                    % the continuous variables
                    [~,discreteDimensions,~] = intersect(allIdxs,obj.discNodeIdxs); discreteDimensions = discreteDimensions';
                    C_discrete_integrate = c;
                    for discreteDimension=discreteDimensions
                        C_discrete_integrate = cumtrapz(u,C_discrete_integrate,discreteDimension);
                    end
                    
                    if(length(parentIdxs)==1)
                        % the density will be used directly, so there
                        % is no need to calculate these copula's (they
                        % don't actually make sense because there is
                        % only one dimension, no concept of "joint")
                        C_parents = [];
                        c_parents = [];
                        C_parents_discrete_integrate = [];
                    else
                        U_in_parents = U_in(:,2:end);
                        % TODO: probably we can just integrate OUT the
                        % child node dimension and get the same value?
                        % speed up some processing this way :D
                        c_parents = empcopulapdf(U_in_parents, obj.h, obj.K, 'betak');
                        C_parents = c_parents; u = linspace(0,1,obj.K);
                        for dim=1:length(parentNames)
                            C_parents = cumtrapz(u,C_parents,dim);
                        end
                        allIdxs = parentIdxs;
                        [~,discreteDimensionsParents,~] = intersect(allIdxs,obj.discNodeIdxs); discreteDimensionsParents = discreteDimensionsParents';
                        C_parents_discrete_integrate = c_parents;
                        for discreteDimension=discreteDimensionsParents
                            C_parents_discrete_integrate = cumtrapz(u,C_parents_discrete_integrate,discreteDimension);
                        end
                    end
                    copFam = hcbnfamily(node, nodeIdx, parentNames, parentIdxs, ...
                            C, c, C_discrete_integrate, C_parents, c_parents, C_parents_discrete_integrate);
                    obj.copulaFamilies{nodeIdx} = copFam;
                end
            end
        end
        
        function [X] = genFamilySamples(obj, node, M)
            %GENFAMILYSAMPLES - generates samples from a specified family
            % Inputs:
            %  node - the child node, whose family is sampled
            %  M - the number of samples to generate
            % Outputs:
            %  X - the samples from the family associated with the input
            %      node.  The columns of X are [node, parent1, parent2 ...]
            
            if(ischar(node))
                % node name was provided
                nodeName = node;
                for ii=1:obj.D
                    if(isequal(nodeName, obj.nodeNames{ii}))
                        nodeIdx = ii;
                        break;
                    end
                end
            else
                % assume node index was provided
                nodeIdx = node;
            end
            
            % generate U from the copula that was "learned"
            copFam = obj.copulaFamilies{nodeIdx};
            U = empcopularnd(copFam.c, M);
            D_family = size(U,2);
            X = zeros(size(U));
            % invert each U appropriately to generate X
            allIdxs = [nodeIdx copFam.parentNodeIdxs];
            for ii=1:M
                for dd=1:D_family
                    X(ii,dd) = obj.empInfo{allIdxs(dd)}.icdf(U(ii,dd));
                end
            end
        end
        
        function [mixedProbability] = computeMixedJointProbability_(obj, X, idxs, nodeNum, parentsFlag)
            %COMPUTEMIXEDPROBABILITY - computes the joint probability of a
            %mixed probability distribution for the indices given by idxs.
            % Inputs:
            %  X - the data point, should be the entire dimension obj.D,
            %      this function picks the appropriate points from that
            %      according to idxs
            %  idxs - the indices of the nodes for which to calculate the
            %         mixed joint probability.
            %  parentsFlag - 1 if we are calculating the mixed prob of the
            %                parents
            
            if(length(idxs)==1)
                mixedProbability = obj.empInfo{idxs}.pdf(X(idxs));
            else
                [~,discreteIdxs,~] = intersect(idxs,obj.discNodeIdxs); discreteIdxs = discreteIdxs';
                continuousIdxs = setdiff(1:length(idxs),discreteIdxs);

                u = zeros(1,length(idxs));

                % fill the u w/ the continuous values, since we
                % already performed the partial derivative w.r.t.
                % the continuous variables for the copula
                for continuousIdx=continuousIdxs
                    continuousNodeNum = idxs(continuousIdx);
                    % query that node's distribution and insert into u
                    u(continuousIdx) = obj.empInfo{continuousNodeNum}.cdf(X(continuousNodeNum));
                end

                % compute the coupla value for the discrete
                % portions through rectangle differencing
                vals = 0:2^length(discreteIdxs) - 1;
                rectangleDiffStates = dec2bin(vals)-'0';
                mixedProbability = 0;
                for ii=1:size(rectangleDiffStates,1)
                    rectangleDiffState = rectangleDiffStates(ii,:);
                    diffStateIdx = 1;
                    for diffState=rectangleDiffState
                        discreteIdx = discreteIdxs(diffStateIdx);
                        discreteNodeNum = idxs(discreteIdx);
                        u(discreteIdx) = obj.empInfo{discreteNodeNum}.cdf(X(discreteNodeNum)-diffState);
                        diffStateIdx = diffStateIdx + 1;
                    end
                    if(parentsFlag)
                        tmp = (-1)^(sum(rectangleDiffState))*empcopulaval(obj.copulaFamilies{nodeNum}.C_parents_discrete_integrate, u, 1/obj.K);
                    else
                        tmp = (-1)^(sum(rectangleDiffState))*empcopulaval(obj.copulaFamilies{nodeNum}.C_discrete_integrate, u, 1/obj.K);
                    end

                    mixedProbability = mixedProbability + tmp;
                end

                % multiply w/ the marginal distributions of the
                % continuous variables
                for continuousIdx=continuousIdxs
                    continuousNodeNum = idxs(continuousIdx);
                    mixedProbability = mixedProbability * obj.empInfo{continuousNodeNum}.pdf(X(continuousNodeNum));
                end
            end
        end
        
        function [conditionalProb] = computeMixedConditionalProbability_(obj, X, idxs, nodeNum)
            jointProbAllNodes = obj.computeMixedJointProbability_(X, idxs, nodeNum, 0);
            jointProbParentNodes = obj.computeMixedJointProbability_(X, idxs(2:end), nodeNum, 1);
            conditionalProb = jointProbAllNodes/jointProbParentNodes;
        end
        
        function [ ll_val, totalProbVec ] = dataLogLikelihood(obj, X)
            M = size(X,1);
            if(size(X,2)~=obj.D)
                error('Input data for LL calculation must be the same dimensions as the BN!');
            end      
            
            if(nargout>1)
                totalProbVec = zeros(1,M);
            end
            
            ll_val = 0;
            for mm=1:M
                % for each D-dimensional data point, compute the likelihood
                totalProb = 1;
                for dd=1:obj.D
                    % get the parents for this node
                    parentIdxs = obj.getParents(dd);
                    if(isempty(parentIdxs))
                        tmp = obj.empInfo{dd}.pdf(X(mm,dd));
                    else
                        allIdxs = [dd parentIdxs];
                        tmp = obj.computeMixedConditionalProbability_(X(mm,:),allIdxs, dd);
                    end
                    totalProb = totalProb * tmp;
                end
                
                if(isnan(totalProb))
                    1;      % for DEBUGGING :D
                end
                if(totalProb<=obj.LOG_CUTOFF)
                    totalProb = obj.LOG_CUTOFF;
                end
                ll_val = ll_val + log(totalProb);
                
                if(nargout>1)
                    totalProbVec(mm) = totalProb;
                end
            end
        end
        
        function [ll_val] = copulaLogLikelihood(obj, X)
            % WARNING - this method should *ONLY* be used if you are
            % dealing w/ all continuous nodes (i.e. using the hcbn as a
            % non-parametric CBN)
            M = size(X,1);
            if(size(X,2)~=obj.D)
                error('Input data for LL calculation must be the same dimensions as the BN!');
            end      
            ll_val = 0;
            for mm=1:M
                for dd=1:obj.D
                    R_ci = obj.computeCopulaRatio(dd, X(mm,:));
                    
                    if(isinf(R_ci) || isnan(R_ci))
                        error('R_ci is inf/nan!');
                    end
                    
                    if(R_ci < obj.LOG_CUTOFF)
                        R_ci = obj.LOG_CUTOFF;
                    end
                    
                    ll_val = ll_val + log(R_ci);
                    if(~isreal(ll_val))
                        error('LL Value imaginary!');
                    end
                end
            end
        end
        
        function [rcVal] = computeCopulaRatio(obj, nodeIdx, x)
            % WARNING - this method should *ONLY* be used if you are
            % dealing w/ all continuous nodes (i.e. using the hcbn as a
            % non-parametric CBN)
            copFam = obj.copulaFamilies{nodeIdx};
            if(isempty(copFam))
                rcVal = 1;
            else
                idxs_all = [nodeIdx copFam.parentNodeIdxs];
                x_all = x(idxs_all);                
                u_all = zeros(1,length(x_all));
                % convert to pseudo-observations via ECDF
                for ii=1:length(x_all)
                    u_all(ii) = obj.empInfo{idxs_all(ii)}.cdf(x_all(ii));
                end                
                u_parents = u_all(2:end);
                
                Rc_num = empcopulaval(copFam.c, u_all);
                if(length(u_parents)==1)
                    Rc_den = 1;
                else
                    Rc_den = empcopulaval(copFam.c_parents, u_parents);
                end
                rcVal = Rc_num/Rc_den;
            end
        end
        
    end
end
