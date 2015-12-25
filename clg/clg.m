classdef clg < handle
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
    end
    
    methods
        function obj = clg(X, discreteNodes, dag)
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
            obj.N = size(X,1);
            obj.D = size(X,2);      
            obj.discreteNodes = discreteNodes;
            obj.dag = dag;
            obj.bnParams = cell(1,size(obj.dag,1));
            obj.uniqueVals = zeros(1,size(obj.dag,1));
            
            for ii=obj.discreteNodes
                obj.uniqueVals(ii) = length(unique(X(:,ii)));
            end
        end
        
        function [] = calcClgParams(obj)
            %CALCCLGPARAMS - computes the CLG parameters for the specified
            %DAG with the data the CLG object was initialized with
            
            numNodes = size(obj.dag,1);
            for node=1:numNodes
                % get the node's parents
                parentNodes = obj.getParents(node);
                
                if(~isempty(parentNodes))
                    % if both the current node and its parents are discrete,
                    % then we can create a joint multinomial distribution,
                    % otherwise the parent must be discrete and the child be
                    % continuous
                    
                    % make a list of all the discrete parents indices
                    discreteParents = intersect(parentNodes, obj.discreteNodes);
                    
                    % make a list of all the continuous parents indices
                    continuousParents = setdiff(obj.discreteNodes, discreteParents);
                    
                    % if continousParents is not empty and our current node
                    % is discrete, this violates the CLG model so we throw
                    % an error
                    if(~isempty(continuousParents) && any(node==obj.discreteNodes))
                        error('CLG model must not have discrete children w/ continuous parents!');
                    else
                        % make combinations for all the parents
                        if(len(discreteParents)==1)
                            combos = 1:obj.uniqueVals(discreteParents(1));
                        elseif(len(discreteParents)==2)
                            combos = combvec(1:obj.uniqueVals(discreteParents(1)), ...
                                             1:obj.uniqueVals(discreteParents(2)));
                        elseif(len(discreteParents)==3)
                            combos = combvec(1:obj.uniqueVals(discreteParents(1)), ...
                                             1:obj.uniqueVals(discreteParents(2)), ...
                                             1:obj.uniqueVals(discreteParents(3)));
                        elseif(len(discreteParents)==4)
                            combos = combvec(1:obj.uniqueVals(discreteParents(1)), ...
                                             1:obj.uniqueVals(discreteParents(2)), ...
                                             1:obj.uniqueVals(discreteParents(3)), ...
                                             1:obj.uniqueVals(discreteParents(3)));
                        else
                            error('Figure out a better syntax to generalize this :)');
                        end
                        combos = combos';
                        
                        % for each combination, estimate the CLG parameter
                        numCombos = size(combos,1);
                        for comboNum=1:numCombos
                            combo = combos(comboNum,:);
                            if(any(node==obj.discreteNodes))
                                % for each unique value, create probability
                                error('Currently unsupported - need to add this functionality!');
                            else
                                % find all the data rows where this combo occurs
                                X_subset = [];
                                for ii=1:obj.N
                                    comboFound = 1;
                                    for jj=1:length(combo)
                                        if(obj.X(ii,discreteParents(jj))~=combo(jj))
                                            comboFound = 0;
                                        end
                                    end
                                    if(comboFound)
                                        X_subset = [X_subset; obj.X(ii,:)];
                                    end
                                end
                                
                                % get all the continuous data associated with this combo
                                continuousNodesIdxs = setdiff(1:D, discreteParents+1);  % +1 for the current node
                                X_subset_continuous = zeros(length(X_subset),length(continuousNodesIdxs));
                                X_subset_continuous(:,1) = X_subset(:,node);
                                idx = 2;
                                for ii=continuousNodesIdxs
                                    X_subset_continuous(:,idx) = X_subset(:,ii);
                                    idx = idx + 1;
                                end

                                if(size(X_subset_continuous,2)==1)
                                    % estimate univariate Gaussian parameters
                                    [muhat,rhohat] = normfit(X_subset_continuous);
                                else
                                    % estimate the Multivariate Gaussian parameters
                                    [Mean, Covariance] = ecmnmle(X_subset_continuous);
                                end
                                
                            end
                        end
                    end
                end
            end
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
        
    end
end