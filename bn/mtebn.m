classdef mtebn < handle
    %MTE Definition of a Mixture of Truncated Exponentials Bayesian Network
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
        N;
        D;
        discreteNodes;
        bnParams;
        X;
        uniqueVals; % the number of unique values for each discrete node
        MTE_DEFAULT_GAUSSIAN;
        
        DEBUG_MODE;
        LOG_CUTOFF;
    end
    
    methods
        function obj = mtebn(X, discreteNodes, varargin)
            % DAG - Constructs a HCBN object
            %  Inputs:
            %   X - a N x D matrix of the the observable data which the
            %       HCBN will model
            %   discreteNodes - an array w/ the indexes of the discrete
            %                   nodes
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
            obj.LOG_CUTOFF = 10^-5;
            
            obj.N = size(X,1);
            obj.D = size(X,2);
            obj.X = X;
            obj.discreteNodes = discreteNodes;
            obj.dag = zeros(obj.D,obj.D);
            obj.bnParams = cell(1,size(obj.dag,1));
            obj.uniqueVals = zeros(1,size(obj.dag,1));
            
            load('mte_default_gaussian.mat');
            obj.MTE_DEFAULT_GAUSSIAN = mte_default_gaussian;
            
            for ii=obj.discreteNodes
                obj.uniqueVals(ii) = length(unique(X(:,ii)));
            end
            
            nVarargs = length(varargin);
            if(nVarargs>0)
                candidateDag = varargin{1};
                if(~acyclic(candidateDag))
                    error('Specified DAG is not acyclic!\n');
                end
                obj.setDag(candidateDag);
            end
        end
        
        function [] = setDag(obj, candidateDag)
            obj.dag = candidateDag;
            obj.calcMteParams();
        end
        
        function [parentIdxs] = getParents(obj, node)
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
            
            nodeIdx = node;
            parentIdxs = find(obj.dag(:,nodeIdx))';
        end
        
        function [] = calcMteParams(obj)
            %CALCCLGPARAMS - computes the CLG parameters for the specified
            %DAG with the data the CLG object was initialized with
            
            numNodes = size(obj.dag,1);
            for node=1:numNodes
                
                if(obj.DEBUG_MODE)
                    fprintf('Processing node %d\n', node);
                end
                
                % get the node's parents
                parentNodes = obj.getParents(node);
                
                if(isempty(parentNodes))
                    X_univariate = obj.X(:,node);
                    % estimate the univariate distribution
                    if(any(node==obj.discreteNodes))
                        isdiscrete = 1;
                        [F,x] = empcdf(X_univariate, isdiscrete);
                        f = emppdf(X_univariate, isdiscrete);
                        
                        empInfoObj = rvEmpiricalInfo(x,f,F,isdiscrete);
                        obj.bnParams{node} = empInfoObj;
                    else
                        % node is continuous, estimate as Mixture of
                        % Truncated Exponentials (MTE) and store paramters
                        mte_params = estMteDensity(X_univariate);
                        obj.bnParams{node} = mte_params;
                    end
                else
                    % if both the current node and its parents are discrete,
                    % then we can create a joint multinomial distribution,
                    % otherwise the parent must be discrete and the child be
                    % continuous
                    
                    % make a list of all the discrete parents indices
                    discreteParents = intersect(parentNodes, obj.discreteNodes);
                    
                    % make a list of all the continuous parents indices
                    continuousParents = setdiff(parentNodes, discreteParents);
                    
                    % if continousParents is not empty and our current node
                    % is discrete, this violates the CLG model so we throw
                    % an error
                    if(~isempty(continuousParents) && any(node==obj.discreteNodes))
                        error('MTE implementation doesnt support continuous parents w/ discrete children at the moment!');
                    end
                    % make combinations for all the parents
                    if(length(discreteParents)==1)
                        combos = 1:obj.uniqueVals(discreteParents(1));
                    else
                        combos = 1:obj.uniqueVals(discreteParents(1));
                        for pNode=discreteParents(2:end)
                            combos = combvec(combos, 1:obj.uniqueVals(discreteParents(pNode)));
                        end
                    end
                    combos = combos';

                    % for each combination, estimate the MTE parameter
                    numCombos = size(combos,1);
                    nodeBnParams = cell(1,numCombos);
                    nodeBnParamsIdx = 1;
                    for comboNum=1:numCombos
                        combo = combos(comboNum,:);
                        if(any(node==obj.discreteNodes))
                            % for each unique value, create probability
                            error('Currently unsupported :/ - need to add this functionality!');
                        else
                            % find all the data rows where this combo occurs
                            X_subset = [];
                            for ii=1:obj.N
                                if(isequal(obj.X(ii,discreteParents),combo))
                                    X_subset = [X_subset; obj.X(ii,:)];
                                end
                            end
                            parentCombinationProbability = size(X_subset,1)/size(obj.X,1);
                            
                            continuousNodesIdxs = [node continuousParents];
                            if(~isempty(X_subset))
%                                 fprintf('Training MTE w/ subset size = %d for combo %s\n', size(X_subset,1), sprintf('%d,',combo'));
                                % get all the continuous data associated with this combo
                                X_subset_continuous = X_subset(:,continuousNodesIdxs);

                                if(length(unique(X_subset_continuous))==1)
                                    mte_info = obj.MTE_DEFAULT_GAUSSIAN;
                                elseif(size(X_subset_continuous,2)==1)
                                    % estimate univariate MTE parameters
                                    mte_info = estMteDensity(X_subset_continuous);
                                else
                                    % estimate the Multivariate MTE parameters
                                    error('MTE -- Currently unsupported!');
                                end
                            else
                                if(size(X_subset_continuous,2)==1)
                                    mte_info = obj.MTE_DEFAULT_GAUSSIAN;
                                else
                                    error('MTE -- Currently unsupported!');
                                end
                            end
                            
                            if(any(isnan(mte_info.density)))
                                1;
                            end
                            
                            nodeBnParams{nodeBnParamsIdx} = mteNodeBnParam(node, combo, parentCombinationProbability, mte_info);
                            nodeBnParamsIdx = nodeBnParamsIdx + 1;
                        end
                    end
                    obj.bnParams{node} = nodeBnParams;
                end
            end
        end
        
        function [llVal, llVec] = dataLogLikelihood(obj, X)
            %DATALOGLIKELIHOOD - calculates the log-likelihood of the given
            %dataset to the calculated model of the data
            % Inputs:
            %  X - the dataset for which to calculate the log-likelihood
            % Outputs:
            %  llVal - the log-likelihood value
            M = size(X,1);
            
            if(nargout>1)
                llVec = zeros(M,1);
            end
            
            llVal = 0;
            for mm=1:M
                familyProb = 1;
                for node=1:obj.D
                    parentNodes = obj.getParents(node);
                    if(isempty(parentNodes))
                        nodeBnParam = obj.bnParams{node};
                        familyProb = familyProb * nodeBnParam.pdf(X(mm,node));
                    else
                        parentNodes = obj.getParents(node);
                        X_parent = X(mm,parentNodes);
                        % find the correct combo
                        numCombos = length(obj.bnParams{node});
                        for combo=1:numCombos
                            if(isequal(X_parent,obj.bnParams{node}{combo}.parentCombination))
                                mte_info = obj.bnParams{node}{combo}.mte_info;
                                parentProb = obj.bnParams{node}{combo}.parentCombinationProbability;
                                familyProb = familyProb * mte_info.pdf(X(mm,node));
                                break;
                            end
                        end
                    end
                end
                if(familyProb<=obj.LOG_CUTOFF)
                    familyProb = obj.LOG_CUTOFF;
                end
                
                if(isnan(familyProb))
                    1;
                end
                
                log_familyProb = log(familyProb);
                llVal = llVal + log_familyProb;
                if(nargout>1)
                    llVec(mm) = log_familyProb;
                end
            end
        end
        
    end
end