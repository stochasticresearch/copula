classdef multinomialbn < handle
    %MULTINOMIALBN Definition of a Multinomial Discrete Bayesian Network
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
        bnParams;
        continuousNodesDiscretizedEdges;
        X;
        X_discretized;
        continuousNodes;
        uniqueVals; % the number of unique values for each discrete node
        
        NUM_DISCRETE_INTERVALS;
        
        DEBUG_MODE;
        LOG_CUTOFF;
    end
    
    methods
        function obj = multinomialbn(X, discreteNodes, varargin)
            % CLG - Constructs a multinomial BN object
            %  Inputs:
            %   X - a N x D matrix of the the observable data which the
            %       HCBN will model
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
            obj.NUM_DISCRETE_INTERVALS = 10;
            
            obj.N = size(X,1);
            obj.D = size(X,2);
            obj.X = X;
            
            % compute the set of continuous marginal distributions
            obj.continuousNodesDiscretizedEdges = cell(1,obj.D);
            obj.continuousNodes = setdiff(1:obj.D,discreteNodes);
            
            obj.X_discretized = obj.X;
            for continuousNode=obj.continuousNodes
                [binnedData, edgeInfo] = discretizeRv(obj.X(:,continuousNode), obj.NUM_DISCRETE_INTERVALS);
                obj.X_discretized(:,continuousNode) = binnedData;
                obj.continuousNodesDiscretizedEdges{continuousNode} = edgeInfo;
            end
            
            obj.uniqueVals = zeros(1,obj.D);
            obj.bnParams = cell(1,size(obj.dag,1));
            for ii=1:obj.D
                obj.uniqueVals(ii) = length(unique(obj.X_discretized(:,ii)));
            end
            
            nVarargs = length(varargin);
            if(nVarargs>0)
                candidateDag = varargin{1};
                if(~acyclic(candidateDag))
                    error('Specified DAG is not acyclic!\n');
                end
                obj.setDag(candidateDag);
            end
            if(nVarargs>1)
                if(isnumeric(varargin{2}))
                    obj.NUM_DISCRETE_INTERVALS = varargin{2};
                end
            end
        end
        
        function [] = setDag(obj, candidateDag)
            obj.dag = candidateDag;
            obj.calcMultinomialParams();
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
        
        function [] = calcMultinomialParams(obj)
            for node=1:obj.D
                parentNodes = obj.getParents(node);
                if(isempty(parentNodes))
                    X_univariate = obj.X_discretized(:,node);
                    % just perform a empirical pdf/cdf fit for this node
                    isdiscrete = 1;
                    [F,x] = empcdf(X_univariate, isdiscrete);
                    f = emppdf(X_univariate, isdiscrete);

                    empInfoObj = rvEmpiricalInfo(x,f,F,isdiscrete);
                    obj.bnParams{node} = empInfoObj;
                else
                    % compute all the combinations for the parent nodes
                    if(length(parentNodes)==1)
                        combos = 1:obj.uniqueVals(parentNodes(1));
                    else
                        combos = 1:obj.uniqueVals(parentNodes(1));
                        for pNode=parentNodes(2:end)
                            combos = combvec(combos, 1:obj.uniqueVals(parentNodes(pNode)));
                        end
                    end
                    combos = combos';
                    numCombos = size(combos,1);
                    nodeBnParams = cell(1,numCombos);
                    nodeBnParamsIdx = 1;
                    for comboNum=1:numCombos
                        combo = combos(comboNum,:);
                        % find this combo in the dataset and compute
                        % associated probability
                        X_discrete_subset = [];
                        for ii=1:obj.N
                            if(isequal(obj.X_discretized(ii,parentNodes),combo))
                                X_discrete_subset = [X_discrete_subset; obj.X_discretized(ii,node)];
                            end
                        end
                        
                        % estimate the discrete domain/pdf/cdf for
                        % X_discrete_subset
                        isdiscrete = 1;
                        if(~isempty(X_discrete_subset))
                            [f,xi] = emppdf(X_discrete_subset, isdiscrete);
                            F = empcdf(X_discrete_subset, isdiscrete);
                        else
                            % this is the default ...
                            xi = 1:10;
                            F = xi/10;
                            f = .1*ones(1,10);
                        end
                        empInfo = rvEmpiricalInfo(xi, f, F,isdiscrete);
                        
                        tmp = struct;
                        tmp.combo = combo;
                        tmp.empInfo = empInfo;
                        nodeBnParams{nodeBnParamsIdx} = tmp;
                        nodeBnParamsIdx = nodeBnParamsIdx + 1;
                    end
                    obj.bnParams{node} = nodeBnParams;
                end
            end
        end
        
        function [discretizedValues] = getDiscreteValueMapping(obj, continuousValues, associatedNode)
            edgeInfo = obj.continuousNodesDiscretizedEdges{associatedNode};
%             discretizedIdxs = discretize(continuousValues, edgeInfo);
            % using below because it is compatible w/ R2014b ... discretize
            % was introduced in R2015a :(
            [~,~,discretizedIdxs] = histcounts(continuousValues, edgeInfo);
            
            % find any NaN's and fix
            % NaN's occur due to out of range exceptions
            badIdxs = [find(isnan(discretizedIdxs)) find(discretizedIdxs==0)];
            belowMinIdxs = continuousValues(badIdxs)<min(edgeInfo);
            aboveMaxIdxs = continuousValues(badIdxs)>max(edgeInfo);
            discretizedIdxs(badIdxs(belowMinIdxs)) = 1;
            discretizedIdxs(badIdxs(aboveMaxIdxs)) = length(edgeInfo);
            discretizedIdxs = floor(discretizedIdxs);
            
            discretizedValues = edgeInfo(discretizedIdxs);
        end
        
        function [llVal] = dataLogLikelihood(obj, X_ll)
            M = size(X_ll,1);
            
            % discretize all the continuous node's in the same manner in
            % which the probability distributions were estimated
            X_discretized_ll = X_ll;
            for continuousNode=obj.continuousNodes
                X_discretized_ll(:,continuousNode) = ...
                    obj.getDiscreteValueMapping(X_ll(:,continuousNode), continuousNode);
            end
            
            llVal = 0;
            for mm=1:M
                familyProb = 1;
                for node=1:obj.D
                    dataPoint = X_discretized_ll(mm,node);
                    
                    nodeBnParam = obj.bnParams{node};
                    if(isa(nodeBnParam, 'rvEmpiricalInfo'))
                        familyProb = familyProb * obj.bnParams{node}.pdf(dataPoint);
                    else
                        % get the parents datapoints
                        parentNodes = obj.getParents(node);
                        parentDataPoints = X_discretized_ll(mm,parentNodes);
                        
                        numCombos = length(nodeBnParam);
                        for comboNum=1:numCombos
                            combo = nodeBnParam{comboNum}.combo;
                            empInfo = nodeBnParam{comboNum}.empInfo;
                            if(isequal(combo,parentDataPoints))
                                familyProb = familyProb * empInfo.pdf(dataPoint);
                                break;
                            end
                        end
                    end
                end
                
                if(familyProb < obj.LOG_CUTOFF)
                    familyProb = obj.LOG_CUTOFF;
                end
                log_familyProb = log(familyProb);
                llVal = llVal + log_familyProb;
            end
        end
        
        function [] = learnStruct_hc(obj, seeddag)
            %LEARNSTRUCT_HC learn the structure of the HCBN network using
            %               the hill climbing algorithm as described in 
            %               Koller and Friedman (2009).  The code for this 
            %               is based off the structure learning toolbox 
            %               within BNT.  See BNT/SLP/learn_struct_hc.m.
            %               This will use the X dataset that was passed to 
            %               the constructor to learn the structure.
            %
            % Inputs:
            %  seeddag - a DAG which can be used as a reference DAG as a
            %            starting point for the search process.  This is
            %            optional, and can be an empty array if no seed is
            %            desired.
            
            % ensure that the seeddag is acyclic
            if(~obj.acyclic(seeddag))
                obj.setDag(zeros(obj.D,obj.D));
            else
                obj.setDag(seeddag);
            end
            
            % get the baseline score
            bestScore = obj.dataLogLikelihood(obj.X);
            done = 0;
            while ~done
                % make dag's which are addition, reversal, and subtraction
                % of edges
                [candidateDags,~,~] = mk_nbrs_of_dag(obj.dag);
                
                % score all the dags
                scores = -Inf*ones(1,length(candidateDags));
                for ii=1:length(candidateDags)
                    obj.setDag(candidateDags{ii});
                    scores(ii) = obj.hcbnLogLikelihood(obj.X);
                end
                
                % find the maximum scoring DAG, and see if it is better
                % than the current best
                maxScore = max(scores);
                % see if multiple candidate dag's had the same maximum
                % score, and if so, choose randomely among those dag's
                new = find(scores == maxScore );
                % update best candidate dag as new dag and continue search
                if ~isempty(new) && (maxScore > bestScore)
                    p = sample_discrete(normalise(ones(1, length(new))));
                    bestScore = maxScore;
                    obj.setDag(candidateDags{new(p)});
                else
                    done = 1;
                end 
            end
            
            % TODO: topo-sort the DAG, and sort the names cell array to
            % match the topologically sorted DAG
            
            % TODO: print out DAG structure
        end
        
    end
end