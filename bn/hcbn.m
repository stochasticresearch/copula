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
        X;          % an N x D matrix of the observable data.  This can be
                    % considered the "training" data.
        X_xform;    % a transformed N x D matrix of the observable data,
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
                                
        DEBUG_MODE; % if turned on, will print out extra stuff to screen to
                    % monitor what is going on and store extra debugging
                    % variables (slower)
        SIM_NUM;    % used for debugging
        TYPE_NAME;  % used for debugging
        
        PSEUDO_OBS_CALC_METHOD;     % should be either RANK or ECDF
    end
    
    methods
        function obj = hcbn(bntPath, X, nodes, discreteNodes, K, h, varargin)
            % HCBN - Constructs a HCBN object
            %  Inputs:
            %   X - a N x D matrix of the the observable data which the
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
            
            obj.DEBUG_MODE = 0;
            obj.SIM_NUM = 0;
            obj.TYPE_NAME = '';
            obj.LOG_CUTOFF = 10^-5;
            
            obj.PSEUDO_OBS_CALC_METHOD = 'RANK';    % can be either RANK or ECDF
            
            % add BNT to the path
            addpath(genpath(bntPath));
            
            obj.D = size(X,2);      
            
            obj.nodeNames = nodes;
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
                obj.setDag(candidateDag);
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
                empInfoObj = rvEmpiricalInfo(x, f, F);
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
                
                if(obj.DEBUG_MODE)
                    fprintf('Estimating Copula for Node=%s:%d <-- ', node, nodeIdx);
                    for jj=1:length(parentNames)
                        fprintf('%d:%s ', parentIdxs(jj), parentNames{jj});
                    end
                    fprintf('\n');
                end
                
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
                        U_in = pseudoobs(X_in, 'ecdf', 100);
                    else
                        U_in = pseudoobs(X_in);
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
                    copFam = copulafamily(node, nodeIdx, parentNames, parentIdxs, ...
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
                        tmp = (-1)^(sum(rectangleDiffState))*empcopulaval(obj.copulaFamilies{nodeNum}.C_parents_discrete_integrate, u);
                    else
                        tmp = (-1)^(sum(rectangleDiffState))*empcopulaval(obj.copulaFamilies{nodeNum}.C_discrete_integrate, u);
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
        
        function [ ll_val ] = dataLogLikelihood(obj, X)
            M = size(X,1);
            if(size(X,2)~=obj.D)
                error('Input data for LL calculation must be the same dimensions as the BN!');
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
            end
        end
        
    end
end
